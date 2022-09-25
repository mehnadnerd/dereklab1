#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
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
        for (x = xnumber * ydim_size; x < (xnumber * ydim_size); ++x) {
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
    double accum;
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
    // debug
    MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
    // end debug
    int periods[1] = {1};
//#ifdef DEBUG
//    printf("precarted %i\n", dims[0]);
//    fflush(stdout);
//#endif
    MPI_Comm topocomm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &topocomm);
    if (topocomm == MPI_COMM_NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
        printf("Node at rank %i is unused\n", rankme);
    }
//#ifdef DEBUG
//    printf("carted one\n");
//    fflush(stdout);
//#endif
    MPI_Comm_rank(topocomm, &rankme);
    MPI_Cart_shift(topocomm, 0, 1, &rankme, &rankright);
    MPI_Cart_shift(topocomm, 0, -1, &rankme, &rankleft);
    //MPI_Cart_shift(topocomm, 1, 1, &rankme, &rankup);
    //MPI_Cart_shift(topocomm, 1, -1, &rankme, &rankdown);

    MPI_Cart_coords(topocomm, rankme, 1, coords);
    ammonarch = coords[0] == 0;
//#ifdef DEBUG
//    printf("mpi carted\n");
//    fflush(stdout);
//#endif

    MPI_Datatype dbl = MPI_DOUBLE;
    int xtag = 0xff;
    //int ytag = 0xf00f;
    int endtag = 0xbeef;
    int deadtag = 0xdead;
    MPI_Request rightsend;
    //MPI_Request downsend;
    MPI_Request leftrecv;
    //MPI_Request uprecv;

    int xdim_size = matrix_dimension_size;
    int ydim_size = matrix_dimension_size / dims[0];

    int mystartx = 0;
    int mystarty = ydim_size * coords[0];
    int matsize = xdim_size * ydim_size;

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
//#ifdef DEBUG
//    printf("alloced\n");
//    fflush(stdout);
//#endif
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
    for (int matrixnum = 1; matrixnum < num_arg_matrices; ++matrixnum) {
        if (gen_sub_matrix(rankme, test_set, matrixnum, bu,
                           mystarty, mystarty + ydim_size - 1, 1,
                           mystartx, mystartx + xdim_size - 1, 1, 0) ==
            NULL) {
            printf("inconsistency in gen_sub_matrix for %i matrix\n", matrixnum);
            exit(1);
        }
        for (int iteration = 0; iteration < dims[0]; ++iteration) {
#ifdef DEBUG
//            printf("matrixnum %i rank %i iteration %i\n", matrixnum, rankme, iteration);
//            fflush(stdout);
#endif
            // start to send/receive
            MPI_Isend(bu, matsize, MPI_DOUBLE,
                       rankright, xtag, topocomm, &rightsend);
//            MPI_Isend(bu, matsize, MPI_DOUBLE,
//                       rankdown, ytag, topocomm, &downsend);
            MPI_Irecv(bi, matsize, MPI_DOUBLE,
                      rankleft, xtag, topocomm, &leftrecv);
//            MPI_Irecv(bi, matsize, MPI_DOUBLE,
//                      rankup, ytag, topocomm, &uprecv);
            // do matrix multiply
            mm_kernel_accum(o, au, bu, xdim_size, ydim_size, iteration, coords[0]);

            // finish send/receive
            MPI_Wait(&rightsend, MPI_STATUS_IGNORE);
            //MPI_Wait(downsend, MPI_STATUS_IGNORE);
            MPI_Wait(&leftrecv, MPI_STATUS_IGNORE);
            //MPI_Wait(uprecv, MPI_STATUS_IGNORE);
            // shuffle matrices
            // swap au and ai, bu and bi
            SWP(au, ai)
            SWP(bu, bi)
        }
        // assign output to a (b will be overwritten at top of loop)
        SWP(au, o)
    }
    o = au;
    if (ammonarch) {
        // i am monarch
        if (debug_perf == 0) {
            print_matrix(o, xdim_size, ydim_size);
            for (int i = 1; i < dims[0]; ++i) {
                int c[1] = {i};
                MPI_Cart_rank(topocomm, c, &rankmonarch);
                MPI_Recv(ai, matsize, MPI_DOUBLE,
                         rankmonarch, endtag, topocomm, MPI_STATUS_IGNORE);
                print_matrix(ai, xdim_size, ydim_size);
            }
        } else {
            double accum = matrix_sum(o, xdim_size, ydim_size);
            double tmpdouble;
            for (int i = 1; i < dims[0]; ++i) {
                int c[1] = {i};
                MPI_Cart_rank(topocomm, c, &rankmonarch);
                MPI_Recv(&tmpdouble, 1, MPI_DOUBLE,
                         rankmonarch, endtag, topocomm, MPI_STATUS_IGNORE);
                accum += tmpdouble;
            }
            printf("%f\n", accum);
        }
    } else {
        // i am non-monarch, need to send to monarch
        int c[1] = {0};
        MPI_Cart_rank(topocomm, c, &rankmonarch);
        if (debug_perf == 0) {
            MPI_Send(o, matsize, MPI_DOUBLE,
                      rankmonarch, endtag, topocomm);
        } else {
            double accum = matrix_sum(o, xdim_size, ydim_size);
            MPI_Send(&accum, 1, MPI_DOUBLE,
                     rankmonarch, endtag, topocomm);
#ifdef DEBUG
            printf("%i %f\n", rankme, accum);
#endif
        }
        // dummy send to stop it from finishing
        MPI_Send(o, 1, MPI_DOUBLE, rankmonarch, deadtag, topocomm);
    }
    return 0;
}