#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define TOP 0
#define BOTTOM 1
/*
 * handle_args
 * root process must do a few extra things before starting
 */
void handle_args(int argc, char **argv, int *s, int *p, int *arrlen)
{
        if (argc != 3) {
                fprintf(stderr, "Error: Wrong number of arguments\n");
                exit(1);
        } else {
                *s = atoi(argv[1]);
                *p = atoi(argv[2]);
                *arrlen = *s * *s;
        }
}

/*
 * generate_array
 * generate an array and fill with random numbers
 */
void generate_array(int length, int randwanted, double *array)
{
        int i;
        if (randwanted) srand(time(NULL));
        for (i = 0; i < length; i++) {
                array[i] = (double) (rand() % 2);
        }
}

void generate_Larray(int size, double *array)
{
        int i, j;
        int count = 0;
        for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                        if (i == 0 || j == 0) {
                                array[count] = (double) 1;
                        } else {
                                array[count] = (double) 0;
                        }
                        count++;
                }
        }
}

/*
 * print_array
 * print to stdout in square form
 */
void print_array(double *array, int nrows, int ncols)
{
        int i, j;
        int index = 0;
        for (i = 0; i < nrows; i++) {
                for (j = 0; j < ncols; j++) {
                        if (i == 0 || i == nrows - 1) {
                                printf("\x1b[32m%f\t\x1b[0m", array[index]);
                        } else {
                                printf("%f\t", array[index]);
                        }
                        index++;
                }
                printf("\n");
        }
        printf("\n");
}

/*
 * get_indexes
 * given the size of the array and the number of processors
 * work out the send counts and displacemnets for the array
 */
void get_indexes(int ncols, int numprocs, int *send_count, int *send_displ, int *recv_count, int *recv_displ)
{
        int inner_ncols = ncols - 2;
        int chunk = inner_ncols / numprocs;
        int remainder = inner_ncols % numprocs;
        int i;

        /* allocate complete rows to processors */
        /* and add 2 overlapping rows for each processor */
        for (i = 0; i < numprocs; i++) {
                send_count[i] = (chunk * ncols) + (2 * ncols);
        }
        i = 0;
        while (remainder > 0) {
                send_count[i] += ncols;
                remainder--;
                i++;
                if (i == numprocs - 1) i = 0;
        }

        /* get the send_displ within the array */
        /* and account for the overlapped rows */
        for (i = 0; i < numprocs; i++) {
                if (i == 0) {
                        send_displ[i] = 0;
                } else {
                        send_displ[i] = send_displ[i - 1] +
                                send_count[i - 1] - (2 * ncols);
                }
        }

        /* work out the recv counts and displacements */
        for (i = 0; i < numprocs; i++) {
                if (i == 0 || i == numprocs - 1) {
                        recv_count[i] = send_count[i] - ncols;
                } else {
                        recv_count[i] = send_count[i] - (2 * ncols);
                }
                if (i == 0) {
                        recv_displ[i] = 0;
                } else {
                        recv_displ[i] = recv_displ[i - 1] + recv_count[i - 1];
                }
        }
}

/*
 * relax
 */
void relax(double *srcarr, double *resarr, int arrlen, int ncols)
{
        int i, j, position;
        int numrows = arrlen / ncols;
        double sum, avg;

        for (i = 1; i < numrows - 1; i++) {
                for (j = 1; j < ncols - 1; j++) {
                        position = (i * ncols) + j;
                        sum = srcarr[position - ncols] +
                                srcarr[position + 1] +
                                srcarr[position + ncols] +
                                srcarr[position - 1];
                        avg = sum / 4;
                        resarr[position] = avg;
                }
        }
}

/*
 * swap
 */
void swap(double **srcarr, double **resarr)
{
        double *tmp;
        tmp = *srcarr;
        *srcarr = *resarr;
        *resarr = tmp;
}

/*
 * send_receive
 * this can communicate all rows in two steps
 * using alternating pairs of processors
 */
void send_receive(double *arr, int *arrlen, int rank, int ncols, int numprocs)
{
        int i;
        MPI_Status STATUS;

        for (i = 0; i < 2; i++) {
                if ((rank + i) % 2 == 0) {
                        /* send/recv with next proc */
                        /*if (rank == numprocs - 1) break;*/
                        MPI_Sendrecv(
                                        &arr[arrlen[rank] - (2 * ncols) + 1],
                                        ncols - 2,
                                        MPI_DOUBLE,
                                        rank + 1,
                                        BOTTOM,
                                        &arr[arrlen[rank] - ncols + 1],
                                        ncols - 2,
                                        MPI_DOUBLE,
                                        rank + 1,
                                        TOP,
                                        MPI_COMM_WORLD,
                                        &STATUS);

                } else {
                        /* send/recv with previous proc */
                        /*if (rank == 0) break;*/
                        MPI_Sendrecv(
                                        &arr[ncols + 1],
                                        ncols - 2,
                                        MPI_DOUBLE,
                                        rank - 1,
                                        TOP,
                                        &arr[1],
                                        ncols - 2,
                                        MPI_DOUBLE,
                                        rank - 1,
                                        BOTTOM,
                                        MPI_COMM_WORLD,
                                        &STATUS);
                }
        }

}

int main(int argc, char **argv)
{
        int rc, numprocs, myrank, namelen;
        char name[MPI_MAX_PROCESSOR_NAME];
        /*MPI_Status STATUS;*/

        rc = MPI_Init(&argc, &argv);
        if (rc != MPI_SUCCESS) {
                fprintf(stderr, "Error starting MPI program\n");
                MPI_Abort(MPI_COMM_WORLD, rc);
                exit(1);
        }
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Get_processor_name(name, &namelen);

        int size, precision, arraylen;
        double *full_array;
        double *source_array;
        double *result_array;
        /* TODO: refactor these into a struct */
        int *send_count = malloc(numprocs * sizeof(int));
        int *recv_count = malloc(numprocs * sizeof(int));
        int *send_displ = malloc(numprocs * sizeof(int));
        int *recv_displ = malloc(numprocs * sizeof(int));

        handle_args(argc, argv, &size, &precision, &arraylen);
        get_indexes(size, numprocs, send_count, send_displ,
                        recv_count, recv_displ);
        if (myrank == 0) {
                full_array = malloc(arraylen * sizeof(double));
                /*generate_array(arraylen, 1, full_array);*/
                generate_Larray(size, full_array);
        } else {
                full_array = NULL;
        }

        if (myrank == 0) print_array(full_array, size, size);
        source_array = malloc(send_count[myrank] * sizeof(double));
        result_array = malloc(send_count[myrank] * sizeof(double));

        MPI_Scatterv(full_array,
                        send_count,
                        send_displ,
                        MPI_DOUBLE,
                        &source_array[0],
                        arraylen,
                        MPI_DOUBLE,
                        0,
                        MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        /* copy contents of source_array into result_array */
        memcpy(result_array, source_array,
                        (send_count[myrank] * sizeof(double)));

        int i = 0;
        do {
                relax(source_array, result_array, send_count[myrank], size);
                swap(&source_array, &result_array);
                if (numprocs > 1) {
                        send_receive(source_array, send_count,
                                        myrank, size, numprocs);
                }
                i++;
        } while (i < 10);

        /* reduce to check when finished? */

        MPI_Barrier(MPI_COMM_WORLD);

        /* Gather pieces of array back together again */
        int offset;
        if (myrank == 0) {
                offset = 0;
        } else {
                offset = size;
        }
        MPI_Gatherv(
                        &source_array[offset],
                        recv_count[myrank],
                        MPI_DOUBLE,
                        full_array,
                        recv_count,
                        recv_displ,
                        MPI_DOUBLE,
                        0,
                        MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        if (myrank == 0) print_array(full_array, size, size);

        free(source_array);
        free(result_array);
        free(send_count);
        free(recv_count);
        free(send_displ);
        free(recv_displ);
        free(full_array);

        MPI_Finalize();

        return 0;
}
