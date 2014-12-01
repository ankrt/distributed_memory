#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define LEN (128)

int main(int argc, char **argv)
{
        int rc, myrank, nproc, namelen;
        char name[LEN];

        rc = MPI_Init(&argc, &argv);
        if (rc != MPI_SUCCESS) {
                printf("Error starting MPI program\n");
                MPI_Abort(MPI_COMM_WORLD, rc);
        }

        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);

        if (myrank == 0) {
                printf("main reports %d procs\n", nproc);
        }

        namelen = LEN;
        MPI_Get_processor_name(name, &namelen);
        printf("hello world %d\n", myrank);

        /* implicit barrier in Finalize */
        /* MPI_Barrier(MPI_COMM_WORLD); */

        MPI_Finalize();
        return 0;
}
