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
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    int periods[2] = {1, 1};
#ifdef DEBUG
    printf("precarted %i\n", dim);
    fflush(stdout);
#endif
    MPI_Comm topocomm;
    int dims[2] = {n, 1};
    int periods[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &topocomm);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &topocomm);
    if (topocomm == MPI_COMM_NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rankme);
        printf("Node at rank %i is unused\n", rankme);
    }
#ifdef DEBUG
    printf("carted one\n");
    fflush(stdout);
#endif
    // TODO: fix this
    int myrank, nextrank, prevrank;
    MPI_Comm_rank(topocomm, &rankme);
    MPI_Cart_shift(topocomm, 1, 1, &rankme, &nextrank);

    MPI_Cart_coords(topocomm, rankme, 2, coords);
#ifdef DEBUG
    printf("mpi carted\n");
    fflush(stdout);
#endif

    MPI_Datatype dbl = MPI_DOUBLE;
    int xtag = 0xff;
    int ytag = 0xf00f;
    int endtag = 0xbeef;
    MPI_Request *rightsend;
    MPI_Request *downsend;
    MPI_Request *leftrecv;
    MPI_Request *uprecv;

    int each_matrixsize = matrix_dimension_size / dim;

    // horizontal, vertical, temp, output
    double *matrices[4];

    // allocate arrays
    for (int i = 0; i < 4; ++i) {
        matrices[i] = (double *) my_calloc(each_matrixsize * matrix_dimension_size), sizeof(double);
    }
#ifdef DEBUG
    printf("alloced\n");
    fflush(stdout);
#endif
    double *horizontal = matrices[0];
    double *vertical = matrix[1];
    double *tmp = matrix[2];
    double *output = matrix[3];

    // get first input matrices
    int horizx = each_matrixsize * coords[0];
    int horizy = each_matrixsize * coords[1];
    gen_sub_matrix(rankme, test_set, 0, horizontal,
                    horizx, horizx + matrix_dimension_size - 1, 1,
                    horizy, horizy + each_matrixsize - 1, 1, 1);
    for (i = 1; i < num_arg_matrices; ++i) {
      gen_sub_matrix(rankme, test_set, i, vertical,
                      horizx, horizx + each_matrixsize - 1, 1,
                      horizy, horizy + matrix_dimension-size - 1, 1, 0);
      for (int iteration = 0; iteration < dim; ++iteration) {
        MPI_Request *send, *rcv;
        MPI_Ibsend(vertical, each_matrixsize*matrix_dimension_size, MPI_DOUBLE,
                    nextrank, 0xdead, topcomm, send);
        MPI_Ircv(tmp, each_matrixsize*matrix_dimension_size, MPI_DOUBLE,
                  prvrank, 0xbeef, topcomm, rcv);
        // TODO: shift output region
        custom_mm(output, horizontal, vertical, each_matrixsize, matrix_dimension_size);

        MPI_WAIT(send, MPI_STATUS_IGNORE);
        MPI_WAIT(rcv, MPI_STATUS_IGNORE);
        // TODO: copy shouldn't be necessary, just alternate
        SWP(vertical, tmp);
      }
      // assign output to horizontal
      SWP(horizontal, output);
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
