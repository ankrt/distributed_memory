#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

/*
 * Generate array and fill with random numbers
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
 * Print array to stdout in square form
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
 * Initialisation that the root process does
 * at the start of the program
 */
void root_pretask(int argc, char **argv)
{
        double *array;
        int size,
            arraylen,
            precision;

        // handle arguments
        if (argc != 3) {
                fprintf(stderr, "Error: Wrong number of arguments\n");
                exit(1);
        } else {
                size = atoi(argv[1]);
                precision = atoi(argv[2]);
                arraylen = size * size;
        }

        // generate array and fill with numbers
        array = generate_array(arraylen, 1);
        print_array(array, size);
        free(array);
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

        printf("numprocs: %d, my rank: %d, running on: %s\n",
                        numprocs,
                        myrank,
                        name);

        // root process needs to do extra initialisation
        //      - create array
        //      - partition array
        //      - scatter array
        if (myrank == 0) {
                root_pretask(argc, argv);
        }
        // enter loop, calculating averages
        // gather up pieces of array
        // finish

        MPI_Finalize();

        return 0;
}
