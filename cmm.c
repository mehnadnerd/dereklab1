#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <assert.h>
#include "gen_matrix.h"
#include "my_malloc.h"

//         0 2
//         1 3
// 0 1
// 2 3
// dimsize = 2

// a is a row major matrix, b is a column major matrix
// emit a row major matrix as output subregion
// a = size0 x size1
static void custom_mm(double *res, double *a, double *b, int size0, int size1)
{
    int offs = 0;
    for (int i = 0; i < size0; i++) {
        int idx = i*offs;
        for (int j = 0; j < size1; j++) {
            res[offs + j] += a[i*size1 + j]*b[i*size1 + j];
        }
        offs += size1;
    }
}

double *buff;
#define SWP(u, i) do { buff = (u); (u) = (i); (i) = buff; } while (0);

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
    // construct grid
    MPI_Comm topocomm;
    int dims[2] = {size, 1};
    int periods[2] = {1, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &topocomm);
    assert(topocomm != MPI_COMM_NULL);

    int myrank, nextrank, prevrank;
    MPI_Comm_rank(topocomm, &myrank);
    MPI_Cart_shift(topocomm, 0, 1, &myrank, &nextrank);
    MPI_Cart_shift(topocomm, 0, -1, &myrank, &prevrank);

    int coords[2]; MPI_Cart_coords(topocomm, myrank, 2, coords);
    int split_size = matrix_dimension_size / size;

    // horizontal submatrix, vertical submatrix, temp submatrix, output submatrix
    double *matrices[4];

    // allocate arrays
    for (int i = 0; i < 4; ++i) {
        matrices[i] = (double *) my_calloc(split_size * matrix_dimension_size, sizeof(double));
    }

    // final result will be in horizontal submatrix
    double *horizontal = matrices[0];
    double *vertical = matrices[1];
    double *tmp = matrices[2];
    double *output = matrices[3];

    // flattened matrix
    int base = split_size * coords[0];
    // get first input matrix
    // horizontal
    //  each core responsible for one row
    gen_sub_matrix(myrank, test_set, 0, horizontal,
                    0, matrix_dimension_size - 1, 1,
                    base, split_size - 1, 1, 1);

    // buffer all prints
    if (debug_perf == 0) {
        bool dummy;
        if (coords[0] == 0) {
            // monarch
            printf("argument matrix 0\n");
        }
        else {
            MPI_Recv(&dummy, 1, MPI_C_BOOL, prevrank, 0xbeef, topocomm, MPI_STATUS_IGNORE);
        }

        for (int i = 0; i < split_size; i++) {
            for (int j = 0; j < matrix_dimension_size; j++) {
                printf("%lf ", horizontal[i*matrix_dimension_size + j]);
            }
            printf("\n");
        }
        MPI_Send(&dummy, 1, MPI_C_BOOL, nextrank, 0xbeef, topocomm);
    }

    for (int i = 1; i < num_arg_matrices; ++i) {
        // handles
        MPI_Request *send, *rcv;
	    // get next matrix vertically
      	gen_sub_matrix(myrank, test_set, i, vertical,
                      	base, split_size - 1, 1,
                      	0, matrix_dimension_size - 1, 1, 0);

        if (debug_perf == 0) {
            bool first = true; bool dummy;
            if (coords[0] == 0) {
                // monarch
                printf("argument matrix %d\n", i);
            }
            for (int i = 0; i < matrix_dimension_size; i++) {
                if (!first || (coords[0] != 0)) {
                    MPI_Recv(&dummy, 1, MPI_C_BOOL, prevrank, 0xbeef, topocomm, MPI_STATUS_IGNORE);
                }
                else {
                    first = false;
                }
                for (int j = 0; j < split_size; j++) {
                    printf("%lf ", vertical[j*matrix_dimension_size + i]);
                }
                if (coords[0] == size - 1) {
                    printf("\n");
                }
                MPI_Send(&dummy, 1, MPI_C_BOOL, nextrank, 0xbeef, topocomm);
            }
        }

        double *begin = output;
      	for (int iteration = 0; iteration < size; ++iteration) {
            // nonblocking send
            if (iteration != (size - 1)) {
                MPI_Ibsend(vertical, split_size*matrix_dimension_size, MPI_DOUBLE,
                            nextrank, 0xbeef, topocomm, send);
                // nonblocking receive
                MPI_Irecv(tmp, split_size*matrix_dimension_size, MPI_DOUBLE,
                          prevrank, 0xbeef, topocomm, rcv);
            }
            // matrix multiply into output subregion
            custom_mm(output, horizontal, vertical, split_size, matrix_dimension_size);
            output += split_size*matrix_dimension_size;

            if (iteration != (size - 1)) {
                MPI_Wait(send, MPI_STATUS_IGNORE);
                MPI_Wait(rcv, MPI_STATUS_IGNORE);
            }
            // alternate buffer
            SWP(vertical, tmp);
        }
        // reset output ptr
        output = begin;
        // row done - swap horizontal and output buffers
        SWP(horizontal, output);
        // zero out output matrix
        for (int j = 0; j < split_size*matrix_dimension_size; j++) {
            output[j] = 0;
        }
    }

    if (debug_perf == 0) {
        bool dummy;
        if (coords[0] == 0) {
            // monarch
            printf("result matrix\n");
        }
        else {
            MPI_Recv(&dummy, 1, MPI_C_BOOL, prevrank, 0xbeef, topocomm, MPI_STATUS_IGNORE);
        }

        for (int i = 0; i < split_size; i++) {
            for (int j = 0; j < matrix_dimension_size; j++) {
                printf("%lf ", horizontal[i*matrix_dimension_size + j]);
            }
            printf("\n");
        }
        MPI_Send(&dummy, 1, MPI_C_BOOL, nextrank, 0xbeef, topocomm);
    }
    return 0;
}
