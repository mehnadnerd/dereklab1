#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include "gen_matrix.h"
#include "my_malloc.h"

//         0 2
//         1 3
// 0 1
// 2 3
// dimsize = 2

void mm_kernel_accum(
        double *__restrict result, // row major
        const double *const __restrict a, // row major
        const double *const __restrict b, // column major
        const int dim_size) {
    int x, y, k;
    for (y = 0; y < dim_size; ++y) {
        for (x = 0; x < dim_size; ++x) {
            double r = result[y * dim_size + x];
            for (k = 0; k < dim_size; ++k) {
                r += a[y * dim_size + k] * b[x * dim_size + k];
            }
            result[y * dim_size + x] = r;
        }
    }
}

double matrix_sum(const double *const __restrict a, // column major
                  const int dim_size) {
    double accum;
    int x, y;
    for (y = 0; y < dim_size; ++y) {
        for (x = 0; x < dim_size; ++x) {
            accum += a[y * dim_size + x];
        }
    }
    return accum;
}

int calc_dims(const int dims[2]) {
    int total = dims[0] * dims[1];
    for (int i = dims[0]; i >= 0; --i) {
        if (i * i < total) {
            return i;
        }
    }
    return 1;
}

void print_matrix(double *result, int dim_size) {
    int x, y;
    for (y = 0; y < dim_size; ++y) {
        for (x = 0; x < dim_size; ++x) {
            printf("%f ", result[y * dim_size + x]);
        }
        printf("\n");
    }
    printf("\n");
}

#define SWP(u, i) do { tmp = (u); (u) = (i); (i) = tmp; } while (0);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rankme, rankleft, rankup, rankright, rankdown;
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
    int dims[2] = {0};
    int coords[2] = {0};
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Dims_create(p, 2, dims);
    // debug
    MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
    if (rankme == 0) {
        printf("dims %i %i\n", dims[0], dims[1]);
    }
    // end debug
    dims[0] = dims[1] = calc_dims(dims); // this makes it so it is guaranteed to be square
    int dim = dims[0];
    int periods[2] = {1, 1};
    MPI_Comm topocomm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &topocomm);
    if (topocomm == MPI_COMM_NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
        printf("Node at rank %i is unused\n", rankme);
        exit(0);
    }
    MPI_Comm_rank(topocomm, &rankme);
    MPI_Cart_shift(topocomm, 0, 1, &rankme, &rankdown);
    MPI_Cart_shift(topocomm, 0, -1, &rankme, &rankup);
    MPI_Cart_shift(topocomm, 1, 1, &rankme, &rankright);
    MPI_Cart_shift(topocomm, 1, -1, &rankme, &rankleft);

    MPI_Cart_coords(topocomm, rankme, 2, coords);

    MPI_Datatype dbl = MPI_DOUBLE;
    int xtag = 0xff;
    int ytag = 0xf00f;
    int endtag = 0xbeef;
    MPI_Request *rightsend;
    MPI_Request *downsend;
    MPI_Request *leftrecv;
    MPI_Request *uprecv;

    int each_matrixsize = matrix_dimension_size / dim;

    int mystartx = each_matrixsize * coords[0];
    int mystarty = each_matrixsize * coords[1];

    double *matrices[5];
    // at start
    // 0/1 are used for a(horizontal)
    // 2/3 are used for b(vertical)
    // 4 is used for output
    double *au;
    double *ai;
    double *bu;
    double *bi;
    double *o;
    double *tmp;

    // allocate arrays
    for (int i = 0; i < 5; ++i) {
        matrices[i] = (double *) my_calloc(each_matrixsize * each_matrixsize, sizeof(double));
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

    // get first two
    if (gen_sub_matrix(rankme, test_set, 0, au,
                       mystartx, mystartx + each_matrixsize - 1, 1,
                       mystarty, mystarty + each_matrixsize - 1, 1, 1) ==
        NULL) {
        printf("inconsistency in gen_sub_matrix for first matrix\n");
        exit(1);
    }
    for (int matrixnum = 1; matrixnum < num_arg_matrices; ++matrixnum) {
        if (gen_sub_matrix(rankme, test_set, matrixnum, bu,
                           mystartx, mystartx + each_matrixsize - 1, 1,
                           mystarty, mystarty + each_matrixsize - 1, 1, 0) ==
            NULL) {
            printf("inconsistency in gen_sub_matrix for %i matrix\n", matrixnum);
            exit(1);
        }
        for (int iteration = 0; iteration < dim; ++iteration) {
            printf("rank %i iteration %i\n", rankme, iteration);
            fflush(stdout);
            // start to send/receive
            MPI_Ibsend(au, each_matrixsize * each_matrixsize, MPI_DOUBLE,
                       rankright, xtag, topocomm, rightsend);
            MPI_Ibsend(bu, each_matrixsize * each_matrixsize, MPI_DOUBLE,
                       rankdown, ytag, topocomm, downsend);
            MPI_Irecv(ai, each_matrixsize * each_matrixsize, MPI_DOUBLE,
                      rankleft, xtag, topocomm, leftrecv);
            MPI_Irecv(bi, each_matrixsize * each_matrixsize, MPI_DOUBLE,
                      rankup, ytag, topocomm, uprecv);
            // do matrix multiply
            mm_kernel_accum(o, au, bu, each_matrixsize);

            // finish send/receive
            MPI_Wait(rightsend, MPI_STATUS_IGNORE);
            MPI_Wait(downsend, MPI_STATUS_IGNORE);
            MPI_Wait(rightsend, MPI_STATUS_IGNORE);
            MPI_Wait(rightsend, MPI_STATUS_IGNORE);
            // shuffle matrices
            // swap au and ai, bu and bi
            SWP(au, ai)
            SWP(bu, bi)
        }
        // assign output to a (b will be overwritten at top of loop)
        SWP(au, o)
    }
    o = au;
    if (coords[0] == 0 && coords[1] == 0) {
        // i am monarch
        if (debug_perf == 0) {
            // TODO
        } else {
            double accum = matrix_sum(o, each_matrixsize);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    if (i == 0 && j == 0) {
                        continue;
                    }
                    // TODO receive value from this node, add to accum
                }
            }
            printf("%f\n", accum);
        }
    } else {
        // i am non-monarch, need to send to monarch
        if (debug_perf == 0) {
            // TODO send things to monarch
        } else {
            double accum = matrix_sum(o, each_matrixsize);
            printf("%i %f\n", rankme, accum);
            // TODO send to monarch
        }
    }

//    if (debug_perf == 0) {
//        // print each of the sub matrices
//        for (i = 0; i < num_arg_matrices; ++i) {
//            printf("argument matrix %d\n", i);
//            print_matrix(r[i], matrix_dimension_size);
//        }
//        printf("result matrix\n");
//        print_matrix(result[n], matrix_dimension_size);
//    } else {
//        double sum = 0.0;
//
//        for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
//            sum += result[n][i];
//        }
//        printf("%f\n", sum);
//    }
    return 0;
}