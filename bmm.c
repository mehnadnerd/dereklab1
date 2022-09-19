#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "gen_matrix.h"
#include "my_malloc.h"

// row major
void mm_kernel_accum(
        double *__restrict result,
        const double *const __restrict a,
        const double *const __restrict b,
        const int dim_size) {
    int x, y, k;
    for (y = 0; y < dim_size; ++y) {
        for (x = 0; x < dim_size; ++x) {
            double r = result[y * dim_size + x];
            for (k = 0; k < dim_size; ++k) {
                r += a[y * dim_size + k] * b[k * dim_size + x];
            }
            result[y * dim_size + x] = r;
        }
    }
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

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    double* a[2];
    double* b[2];
    double* result[2];
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
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Dims_create(p, 2, dims);
    int periods[2] = {1,1};
    MPI_Comm topocomm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &topocomm);

    int each_matrixdimsize = matrix_dimension_size /

    MPI_Cart_create(comm_old, 2, const int dims[],
    const int periods[], int reorder, MPI_Comm *comm_cart);
    // allocate arrays
    for(int i = 0; i < 2; ++i) {
        a[i] = (double *)
    }
    result[0] = (double *) my_calloc(sizeof(double), matrix_dimension_size * matrix_dimension_size);
    result[1] = (double *) my_calloc(sizeof(double), matrix_dimension_size * matrix_dimension_size);

    // get sub matrices
    for (int i = 0; i < num_arg_matrices; ++i) {
        r[i] = (double *) my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
        if (gen_sub_matrix(0, test_set, i, r[i], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) ==
            NULL) {
            printf("inconsistency in gen_sub_matrix\n");
            exit(1);
        }
    }

    // perform matrix multiplies
    int n = 0;

    mm_kernel_accum(result[0], r[0], r[1], matrix_dimension_size);
    for (i = 2; i < num_arg_matrices; ++i) {
        mm_kernel_accum(result[n ^ 0x1], result[n], r[i], matrix_dimension_size);
        n = n ^ 0x1;
    }

    if (debug_perf == 0) {
        // print each of the sub matrices
        for (i = 0; i < num_arg_matrices; ++i) {
            printf("argument matrix %d\n", i);
            print_matrix(r[i], matrix_dimension_size);
        }
        printf("result matrix\n");
        print_matrix(result[n], matrix_dimension_size);
    } else {
        double sum = 0.0;

        for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
            sum += result[n][i];
        }
        printf("%f\n", sum);
    }
    return 0;
}