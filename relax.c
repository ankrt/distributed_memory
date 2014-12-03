#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

// displacements and chunk size
struct scatter_spec {
        int *displacement, *send_count;
};

/*
 * handle_args
 * root process must do a few extra things before starting
 */
void handle_args(int argc, char **argv, int* size, int* precision, int* arraylen)
{
        if (argc != 3) {
                fprintf(stderr, "Error: Wrong number of arguments\n");
                exit(1);
        } else {
                *size = atoi(argv[1]);
                *precision = atoi(argv[2]);
                *arraylen = *size * *size;
        }
}

/*
 * generate_array
 * generate an array and fill with random numbers
 */
double* generate_array(int length, int randwanted)
{
        int i;
        double *array = malloc(length * sizeof(double));
        if (randwanted) srand(time(NULL));
        for (i = 0; i < length; i++) {
                array[i] = (double) (rand() % 2);
        }
        return array;
}

/*
 * print_array
 * print to stdout in square form
 */
void print_array(double *array, int size)
{
        int i, j;
        int index = 0;
        for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                        printf("%f\t", array[index]);
                        index++;
                }
                printf("\n");
        }
        printf("\n");
}

/*
 * get_ranges
 * given the size of the array and the number of processors
 * work out the row ranges that each processor will be sent
 */
struct scatter_spec* get_scatter_spec(int rowlen, int numprocs)
{
        int inner_rowlen = rowlen - 2;
        int chunk = inner_rowlen / numprocs;
        int remainder = inner_rowlen % numprocs;
        int i;
        struct scatter_spec *scsp = malloc(sizeof(struct scatter_spec));

        int *allocation = malloc(numprocs * sizeof(int));
        int *dsp = malloc(numprocs * sizeof(int));

        // allocate complete rows to processors
        for (i = 0; i < numprocs; i++) {
                allocation[i] = chunk * rowlen;
        }
        i = 0;
        while (remainder > 0) {
                allocation[i] += rowlen;
                remainder--;
                i++;
                if (i == numprocs - 1) i = 0;
        }
        allocation[0] += rowlen;
        allocation[numprocs - 1] += rowlen;

        // get the displacement within the array
        for (i = 0; i < numprocs; i++) {
                if (i == 0) {
                        dsp[i] = 0;
                } else {
                        dsp[i] = dsp[i - 1] + allocation[i - 1];
                }
        }
        scsp->displacement = dsp;
        scsp->send_count = allocation;
        return scsp;
}

/*void lol()*/
/*{*/
        /*array = generate_array(arraylen, 1);*/
        /*scsp = get_scatter_spec(size, numprocs);*/

        /*MPI_Scatterv(*/
                        /*array,*/
                        /*scsp->send_count,*/
                        /*scsp->displacement,*/
                        /*MPI_DOUBLE,*/
                        /*0,*/
/*}*/

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

        /*printf("numprocs: %d, my rank: %d, running on: %s\n",*/
                        /*numprocs,*/
                        /*myrank,*/
                        /*name);*/

        int size, precision, arraylen;
        double *sen_array; // send buffer
        double *rec_array; // receive buffer
        struct scatter_spec *scsp; // contains send counts and displacements
        handle_args(argc, argv, &size, &precision, &arraylen);
        scsp = get_scatter_spec(size, numprocs);
        rec_array = malloc(scsp->send_count[myrank] * sizeof(double));
        if (myrank == 0) {
                sen_array = generate_array(size, 1);
        } else {
                sen_array = NULL;
        }

        MPI_Scatterv(sen_array, scsp->send_count, scsp->displacement, MPI_DOUBLE,
                        &rec_array[0], arraylen, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        printf("I will get %d elements\n", scsp->send_count[myrank]);

        free(sen_array);
        free(scsp->displacement);
        free(scsp->send_count);
        free(scsp);
        // receive scattered pieces

        // enter loop, calculating averages
        // gather up pieces of array
        // finish

        MPI_Finalize();

        return 0;
}
