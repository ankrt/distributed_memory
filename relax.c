#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

// displacements and chunk size
struct scatter_spec {
        int *displs, *allocation;
};

/*
 * generate_array
 * generate an array and fill with random numbers
 * TODO: choose whether seed is on or off
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
        struct scatter_spec *sspc = malloc(sizeof(struct scatter_spec));

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

        sspc->displs = dsp;
        sspc->allocation = allocation;
        return sspc;
}

/*
 * root_pretask
 * root process must do a few extra things before starting
 */
void root_pretask(int argc, char **argv, int numprocs)
{
        double *array;
        int size, arraylen, precision;
        struct scatter_spec *sspc;

        // handle arguments
        if (argc != 3) {
                fprintf(stderr, "Error: Wrong number of arguments\n");
                exit(1);
        } else {
                size = atoi(argv[1]);
                precision = atoi(argv[2]);
                arraylen = size * size;
        }

        array = generate_array(arraylen, 1);
        // get the displacements to pass to the scatter function
        sspc = get_scatter_spec(size, numprocs);
        print_array(array, size);

        int i;
        for (i = 0; i < numprocs; i++) {
                printf("proc %d starts at %d, using %d elements\n",
                                i,
                                sspc->displs[i],
                                sspc->allocation[i]);
        }

        free(array);
        free(sspc->displs);
        free(sspc->allocation);
        free(sspc);
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

        /*printf("numprocs: %d, my rank: %d, running on: %s\n",*/
                        /*numprocs,*/
                        /*myrank,*/
                        /*name);*/

        if (myrank == 0) {
                root_pretask(argc, argv, numprocs);
        }
        // enter loop, calculating averages
        // gather up pieces of array
        // finish

        MPI_Finalize();

        return 0;
}
