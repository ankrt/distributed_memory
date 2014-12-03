#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

/*
 * handle_args
 * root process must do a few extra things before starting
 */
void handle_args(int argc, char **argv, int* s, int* p, int* alen)
{
        if (argc != 3) {
                fprintf(stderr, "Error: Wrong number of arguments\n");
                exit(1);
        } else {
                *s = atoi(argv[1]);
                *p = atoi(argv[2]);
                *alen = *s * *s;
        }
}

/*
 * generate_array
 * generate an array and fill with random numbers
 */
void generate_array(int length, int randwanted, double* array)
{
        int i;
        if (randwanted) srand(time(NULL));
        for (i = 0; i < length; i++) {
                array[i] = (double) (rand() % 2);
        }
}

/*
 * print_array
 * print to stdout in square form
 */
void print_array(double* array, int nrow, int ncol)
{
        int i, j;
        int index = 0;
        for (i = 0; i < nrow; i++) {
                for (j = 0; j < ncol; j++) {
                        printf("%f\t", array[index]);
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
void get_indexes(int rowlen, int numprocs, int* send_count, int* displacement)
{
        int inner_rowlen = rowlen - 2;
        int chunk = inner_rowlen / numprocs;
        int remainder = inner_rowlen % numprocs;
        int i;

        /* allocate complete rows to processors */
        for (i = 0; i < numprocs; i++) {
                send_count[i] = chunk * rowlen;
        }
        i = 0;
        while (remainder > 0) {
                send_count[i] += rowlen;
                remainder--;
                i++;
                if (i == numprocs - 1) i = 0;
        }
        /* first and last ranked will have direct access to their */
        /* top and bottom rows respectively - no messages needed */
        send_count[0] += rowlen;
        send_count[numprocs - 1] += rowlen;

        /* get the displacement within the array */
        for (i = 0; i < numprocs; i++) {
                if (i == 0) {
                        displacement[i] = 0;
                } else {
                        displacement[i] = displacement[i - 1] +
                                send_count[i - 1];
                }
        }
}

/*
 * relax
 */
        /* THIS WILL NOT WORK UNTIL MESSAGE PASSING IS IMPLEMENTED */
void relax(double* srcarr, double* resarr, int arrlen, int numcols)
{
        int i, start, end;
        int numrows = arrlen / numcols;
        double sum, avg;

        /* THIS WILL NOT WORK UNTIL MESSAGE PASSING IS IMPLEMENTED */

        for (i = start; i < end; i++) {
                sum = srcarr[(i / numcols - 1) + (i % numcols)] +
                        srcarr[(i / numcols) + (i % numcols + 1)] +
                        srcarr[(i / numcols + 1) + (i % numcols)] +
                        srcarr[(i / numcols) + (i % numcols - 1)];
                avg = sum / 4;
        }
}

int main(int argc, char **argv)
{
        int rc, numprocs, myrank, namelen;
        char name[MPI_MAX_PROCESSOR_NAME];

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
        double *full_array, *source_array;
        int *send_count = malloc(numprocs * sizeof(int));
        int *displacement = malloc(numprocs * sizeof(int));

        handle_args(argc, argv, &size, &precision, &arraylen);
        get_indexes(size, numprocs, send_count, displacement);
        source_array = malloc(send_count[myrank] * sizeof(double));
        if (myrank == 0) {
                full_array = malloc(arraylen * sizeof(double));
                generate_array(arraylen, 1, full_array);
                print_array(full_array, size, size);
        } else {
                full_array = NULL;
        }

        MPI_Scatterv(full_array,
                        send_count,
                        displacement,
                        MPI_DOUBLE,
                        &source_array[0],
                        arraylen,
                        MPI_DOUBLE,
                        0,
                        MPI_COMM_WORLD);
        free(full_array);

        /* wait for everyone to recieve their pieces */
        MPI_Barrier(MPI_COMM_WORLD);

        /* make a copy of initial array, this will be the results array */
        double *result_array = malloc(send_count[myrank] * sizeof(double));

        /*relax(source_array,*/
                        /*result_array,*/
                        /*send_count[myrank],*/
                        /*size);*/
        free(source_array);
        free(result_array);
        free(displacement);
        free(send_count);

        // enter loop, calculating averages
        // gather up pieces of array
        // finish

        MPI_Finalize();

        return 0;
}
