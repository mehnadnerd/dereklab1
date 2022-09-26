#include <stdio.h>
#include <stdlib.h>
#include <cilk.h>
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

void mm_kernel_accum_cilq(
        double *__restrict result, // row major
        const double *const __restrict a, // row major
        const double *const __restrict b, // column major
        const int fullsize,
        const int xdim_size,
        const int ydim_size,
        const int xnumber,
        const int ynumber
) {
    int x, y, k;
    for (y = ynumber * ydim_size; y < (ynumber + 1) * ydim_size; ++y) {
        for (x = xnumber * xdim_size; x < (xnumber + 1) * xdim_size; ++x) {
            double r = result[y * fullsize + x];
            for (k = 0; k < fullsize; ++k) {
                r += a[y * fullsize + k] * b[x * fullsize + k];
            }
            result[y * fullsize + x] = r;
//            printf("%i %i %f\n", x, y, r);
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

    int xdim, ydim;
    if (matrix_dimension_size < 32) {
        xdim = 2;
        ydim = 2;
    }
    else {
        xdim = 32; // number from ass
        ydim = 32;
    }

    //stolen from https://web.cels.anl.gov/~thakur/sc16-mpi-tutorial/slides.pdf

    int xdim_size = matrix_dimension_size / xdim;
    int ydim_size = matrix_dimension_size / ydim;

    int matsize = xdim_size * ydim_size;

    double *matrices[3];
    double *a;
    double *b;
    double *o;
    double *tmp;

    // allocate arrays
    for (int i = 0; i < 3; ++i) {
        matrices[i] = (double *) my_calloc(matrix_dimension_size * matrix_dimension_size, sizeof(double));
    }

    a = matrices[0];
    b = matrices[1];
    o = matrices[2];
    if (a == NULL || b == NULL || o == NULL) {
        printf("mem failed\n");
        exit(1);
    }

    // get first a matrix
    if (gen_sub_matrix(0, test_set, 0, a,
                       0, matrix_dimension_size - 1, 1,
                       0, matrix_dimension_size - 1, 1, 1) ==
        NULL) {
        printf("inconsistency in gen_sub_matrix for first matrix\n");
        exit(1);
    }
    if (debug_perf == 0) {
        printf("argument matrix 0\n");
        for (int i = 0; i < matrix_dimension_size; i++) {
            for (int j = 0; j < matrix_dimension_size; j++) {
                printf("%lf ", a[i*matrix_dimension_size + j]);
            }
            printf("\n");
        }
    }
//    printf("Begin first print for %i\n", rankme);
//    debug_print_matrix(au, xdim_size, ydim_size, rankme);
    for (int matrixnum = 1; matrixnum < num_arg_matrices; ++matrixnum) {
        if (gen_sub_matrix(0, test_set, matrixnum, b,
                           0, matrix_dimension_size - 1, 1,
                           0, matrix_dimension_size - 1, 1, 0) ==
            NULL) {
            printf("inconsistency in gen_sub_matrix for %i matrix\n", matrixnum);
            exit(1);
        }
        if (debug_perf == 0) {
            printf("argument matrix %d\n", matrixnum);
            for (int i = 0; i < matrix_dimension_size; i++) {
                for (int j = 0; j < matrix_dimension_size; j++) {
                    printf("%lf ", b[j*matrix_dimension_size + i]);
                }
                printf("\n");
            }
        }
        for (int i = 0; i < xdim; ++i) {
            for (int j = 0; j < ydim; ++j) {
                cilk_spawn mm_kernel_accum_cilq(o, a, b, matrix_dimension_size, xdim_size, ydim_size, i, j);
            }
        }
        cilk_sync;
        SWP(o, a)
        memset(o, 0, matrix_dimension_size * matrix_dimension_size * sizeof(double ));
//        print_matrix(a, matrix_dimension_size, matrix_dimension_size);
//        print_matrix(b, matrix_dimension_size, matrix_dimension_size);
    }
    o = a;
    if (debug_perf == 0) {
        printf("result matrix\n");
        for (int i = 0; i < matrix_dimension_size; i++) {
            for (int j = 0; j < matrix_dimension_size; j++) {
                printf("%lf ", o[i*matrix_dimension_size + j]);
            }
            printf("\n");
        }
    }
    else {
        double accum = matrix_sum(o, matrix_dimension_size, matrix_dimension_size);
        printf("result sum %f\n", accum);
    }
    return 0;
}
