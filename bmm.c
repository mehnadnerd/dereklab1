#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>
#include "gen_matrix.h"
#include "my_malloc.h"

//         0 4 | 8 c
//         1 5 | 9 d
//         2 6 | a e
//         3 7 | b f
// 0 1 2 3           // ynumber 0
// 4 5 6 7
// -------
// 8 9 a b           // ynumber 1
// c d e f
//         xnumber
//         0      1
// xdimsize = 4
// ydimsize = 2

void debug_print_matrix(double *result, int xdim_size, int ydim_size, int rank) {
    int x, y;
    for (y = 0; y < ydim_size; ++y) {
        printf("r%i ", rank);
        for (x = 0; x < xdim_size; ++x) {
            printf("%f ", result[y * xdim_size + x]);
        }
        printf("\n");
    }
    printf("\n");
}

void mm_kernel_accum(
        double *__restrict result, // row major
        const double *const __restrict a, // row major
        const double *const __restrict b, // column major
        const int xdim_size,
        const int ydim_size,
        const int xnumber,
        const int ynumber
        // n.b. we don't care about ynumber/y because we would have to offset and then unoffset it
) {
    int x, y, k, affinex;
    for (y = 0; y < ydim_size; ++y) {
        for (x = xnumber * ydim_size; x < (xnumber * ydim_size) + ydim_size; ++x) {
            affinex = x - xnumber * ydim_size;
            double r = result[y * xdim_size + x];
            for (k = 0; k < xdim_size; ++k) {
                r += a[y * xdim_size + k] * b[affinex * xdim_size + k];
            }
            result[y * xdim_size + x] = r;
        }
    }
}

double matrix_sum(const double *const __restrict a, // row major
                  const int xdim_size, const int ydim_size) {
    double accum = 0.0;
    int x, y;
    for (y = 0; y < ydim_size; ++y) {
        for (x = 0; x < xdim_size; ++x) {
            accum += a[y * xdim_size + x];
        }
    }
    return accum;
}

void print_matrix(double *result, int xdim_size, int ydim_size) {
    int x, y;
    for (y = 0; y < ydim_size; ++y) {
        for (x = 0; x < xdim_size; ++x) {
            printf("%f ", result[y * xdim_size + x]);
        }
        printf("\n");
    }
    printf("\n");
}

#define SWP(u, i) do { tmp = (u); (u) = (i); (i) = tmp; } while (0);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    bool ammonarch;
    int rankme, rankleft, rankup, rankright, rankdown, rankmonarch;
    int num_arg_matrices;
    int c[1] = {0};

    if (argc != 4) {
        printf("usage: debug_perf test_set matrix_dimension_size\n");
        exit(1);
    }
    int debug_perf = atoi(argv[1]);
    int test_set = atoi(argv[2]);
    matrix_dimension_size = atoi(argv[3]);
    num_arg_matrices = init_gen_sub_matrix(test_set);

    //stolen from https://web.cels.anl.gov/~thakur/sc16-mpi-tutorial/slides.pdf
    int p;
    int dims[1] = {0};
    int coords[1] = {0};
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Dims_create(p, 1, dims);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
    int periods[1] = {1};

    MPI_Comm topocomm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &topocomm);
    if (topocomm == MPI_COMM_NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
        printf("Node at rank %i is unused\n", rankme);
    }

    MPI_Comm_rank(topocomm, &rankme);
    c[0] = rankme + 1;
    MPI_Cart_rank(topocomm, c, &rankright);
    c[0] = rankme - 1;
    MPI_Cart_rank(topocomm, c, &rankleft);
    c[0] = 0;
    MPI_Cart_rank(topocomm, c, &rankmonarch);
    MPI_Cart_coords(topocomm, rankme, 1, coords);
    printf("I am m %i l %i r %i m %i\n", rankme, rankleft, rankright, rankmonarch);
    ammonarch = rankme == rankmonarch;

    int xtag = 0x00ff;
    int endtag = 0xbeef;
    MPI_Request rightsend;
    MPI_Request leftrecv;

    int xdim_size = matrix_dimension_size;
    int ydim_size = matrix_dimension_size / dims[0];

    int mystartx = 0;
    int mystarty = ydim_size * coords[0];
    int matsize = xdim_size * ydim_size;
    printf("I am rank %i xdim %i ydim %i matsize %i starty %i\n", rankme, xdim_size, ydim_size, matsize, mystarty);

    double *matrices[5];
    // at start
    // 0/1 are used for a(horizontal)
    // 2/3 are used for b(vertical)
    // 4 is used for output
    double *au;
    double *ai; // not used
    double *bu;
    double *bi;
    double *o;
    double *tmp;

    // allocate arrays
    for (int i = 0; i < 5; ++i) {
        matrices[i] = (double *) my_calloc(matsize, sizeof(double));
    }

    au = matrices[0];
    ai = matrices[1];
    bu = matrices[2];
    bi = matrices[3];
    o = matrices[4];
    if (au == NULL || ai == NULL || bu == NULL || bi == NULL || o == NULL) {
        printf("mem failed\n");
        exit(1);
    }

    // get first a matrix
    if (gen_sub_matrix(rankme, test_set, 0, au,
                       mystartx, mystartx + xdim_size - 1, 1,
                       mystarty, mystarty + ydim_size - 1, 1, 1) ==
        NULL) {
        printf("inconsistency in gen_sub_matrix for first matrix\n");
        exit(1);
    }
    printf("Begin first print for %i\n", rankme);
    debug_print_matrix(au, xdim_size, ydim_size, rankme);
    for (int matrixnum = 1; matrixnum < num_arg_matrices; ++matrixnum) {
        if (gen_sub_matrix(rankme, test_set, matrixnum, bu,
                           mystarty, mystarty + ydim_size - 1, 1,
                           mystartx, mystartx + xdim_size - 1, 1, 0) ==
            NULL) {
            printf("inconsistency in gen_sub_matrix for %i matrix\n", matrixnum);
            exit(1);
        }
        printf("Begin print for %i\n", rankme);
        debug_print_matrix(bu, xdim_size, ydim_size, rankme);
        for (int iteration = 0; iteration < dims[0]; ++iteration) {
#ifdef DEBUG
            //            printf("matrixnum %i rank %i iteration %i\n", matrixnum, rankme, iteration);
            //            fflush(stdout);
#endif
            // start to send/receive
            MPI_Isend(bu, matsize, MPI_DOUBLE,
                      rankright, xtag, topocomm, &rightsend);
            MPI_Irecv(bi, matsize, MPI_DOUBLE,
                      rankleft, xtag, topocomm, &leftrecv);
            // do matrix multiply
            printf("Begin printa for %i\n", rankme);
            debug_print_matrix(au, xdim_size, ydim_size, rankme);
            printf("Begin printb for %i\n", rankme);
            debug_print_matrix(bu, xdim_size, ydim_size, rankme);
            mm_kernel_accum(o, au, bu, xdim_size, ydim_size, (iteration + coords[0]) % dims[0], coords[0]);
            printf("Begin printo for %i\n", rankme);
            debug_print_matrix(o, xdim_size, ydim_size, rankme);

            // finish send/receive
            MPI_Wait(&rightsend, MPI_STATUS_IGNORE);
            MPI_Wait(&leftrecv, MPI_STATUS_IGNORE);
            // shuffle matrices
            // swap bu and bi
            SWP(bu, bi)
        }
        // assign output to a (b will be overwritten at top of loop)
        SWP(au, o)
        memset(o, 0, matsize * sizeof(double));
        printf("Begin inter print for %i\n", rankme);
        debug_print_matrix(au, xdim_size, ydim_size, rankme);
    }
    o = au;
    printf("Begin final print for %i\n", rankme);
    debug_print_matrix(o, xdim_size, ydim_size, rankme);
    double accum = matrix_sum(o, xdim_size, ydim_size);
    printf("%i local sum %f\n", rankme, accum);
    MPI_Gather(&accum, 1, MPI_DOUBLE, ai, 1, MPI_DOUBLE, rankmonarch, MPI_COMM_WORLD);
    printf("I am %i and monarch is %i\n", rankme, rankmonarch);
    if (ammonarch) {
        printf("Begin ai print\n");
        debug_print_matrix(ai, dims[0], 1, 100+rankme);
        double newaccum = matrix_sum(ai, dims[0], 1);
        printf("Reduce sum %f\n", newaccum);
    } else {
        printf("I'm done %i\n", rankme);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}