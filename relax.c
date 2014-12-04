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
void print_array(double *array, int nrow, int ncol)
{
        int i, j;
        int index = 0;
        for (i = 0; i < nrow; i++) {
                for (j = 0; j < ncol; j++) {
                        if (i == 0 || i == nrow - 1) {
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
void get_indexes(int rowlen, int numprocs, int *arrlen, int *dsp)
{
        int inner_rowlen = rowlen - 2;
        int chunk = inner_rowlen / numprocs;
        int remainder = inner_rowlen % numprocs;
        int i;

        /* allocate complete rows to processors */
        for (i = 0; i < numprocs; i++) {
                arrlen[i] = (chunk * rowlen) + (2 * rowlen);
        }
        i = 0;
        while (remainder > 0) {
                arrlen[i] += rowlen;
                remainder--;
                i++;
                if (i == numprocs - 1) i = 0;
        }

        /* get the dsp within the array */
        /* and account for the overlapped rows */
        for (i = 0; i < numprocs; i++) {
                if (i == 0) {
                        dsp[i] = 0;
                } else {
                        dsp[i] = dsp[i - 1] +
                                arrlen[i - 1] - (2 * rowlen);
                }
        }
}

/*
 * relax
 */
void relax(double *srcarr, double *resarr, int arrlen, int numcols)
{
        int i, j, position;
        int numrows = arrlen / numcols;
        double sum, avg;

        for (i = 1; i < numrows - 1; i++) {
                for (j = 1; j < numcols - 1; j++) {
                        position = (i * numcols) + j;
                        sum = srcarr[position - numcols] +
                                srcarr[position + 1] +
                                srcarr[position + numcols] +
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
 * get_overlap
 * return the number of 'extra' elements (rows) that a process has
 */
int get_overlap(int rank, int numprocs, int ncols)
{
        if (rank == 0 || rank == numprocs - 1) {
                return ncols;
        } else {
                return 2 * ncols;
        }
}

/*
 * send_receive
 */
void send_receive(double *arr, int *arrlen, int rank, int ncols, int numprocs)
{
        int i;
        MPI_Status STATUS;

        for (i = 0; i < 2; i++) {
                if ((rank + i) % 2 == 0) {
                        // send/rec with next proc
                        if (rank == numprocs - 1) break;
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
                        // send/rec with previous proc
                        if (rank == 0) break;
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
        int *myarraylen = malloc(numprocs * sizeof(int));
        int *displacement = malloc(numprocs * sizeof(int));

        handle_args(argc, argv, &size, &precision, &arraylen);
        get_indexes(size, numprocs, myarraylen, displacement);
        if (myrank == 0) {
                full_array = malloc(arraylen * sizeof(double));
                /*generate_array(arraylen, 1, full_array);*/
                generate_Larray(size, full_array);
        } else {
                full_array = NULL;
        }

        source_array = malloc(myarraylen[myrank] * sizeof(double));
        result_array = malloc(myarraylen[myrank] * sizeof(double));

        MPI_Scatterv(full_array,
                        myarraylen,
                        displacement,
                        MPI_DOUBLE,
                        &source_array[0],
                        arraylen,
                        MPI_DOUBLE,
                        0,
                        MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        /* copy contents of source_array into result_array */
        memcpy(result_array, source_array,
                        (myarraylen[myrank] * sizeof(double)));

        /*if (myrank == 0) print_array(source_array, (myarraylen[myrank] / size), size);*/
        int i = 0;
        do {
                relax(source_array, result_array, myarraylen[myrank], size);
                swap(&source_array, &result_array);
                print_array(source_array, (myarraylen[myrank] / size), size);
                if (numprocs > 1) {
                        send_receive(source_array, myarraylen,
                                        myrank, size, numprocs);
                }
                i++;
        } while (i < 1);

        /*[>print_array(source_array, (myarraylen[myrank] / size), size);<]*/
        /*// reduce to*/

        /*MPI_Barrier(MPI_COMM_WORLD);*/

        /*[> Gather pieces of array back together again <]*/
        /*[> update own indexing <]*/
        /*int overlap = get_overlap(myrank, numprocs, size);*/


        free(source_array);
        free(result_array);
        free(displacement);
        free(myarraylen);
        free(full_array);

        MPI_Finalize();

        return 0;
}
