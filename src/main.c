

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#include "global.h"

/*! \file main.c
 *  \brief contains the main simulation loop
 */

/*!
 * This function intializes MPI, calls Timothy initialization and allocation functions.
 * It also contains the main simulation loop where all important simulation steps are called.
 */

int main(int argc, char **argv)
{
        system_t system;
        settings_t settings;
        celltype_t *celltype;
        environment_t *environment;
        cellsinfo_t cellsinfo;
        grid_t grid;
        commdata_t commdata;
        statistics_t statistics;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &system.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &system.rank);
        system.nthreads = omp_get_max_threads();

        getsysteminfo(&system);
        printinfo(system);
        initialisation(argc,argv,&system,&settings,&celltype,&environment);
        allocatecells(system,settings,celltype,&cellsinfo);
        allocategrid(system,settings,&grid);
        allocatefields(system,settings,&environment);
        lbinit(argc,argv,MPI_COMM_WORLD,system,&cellsinfo);

        for (settings.step = 0; settings.step < settings.numberofsteps; settings.step++) {
                updateglobalcounts(&cellsinfo);
                lbexchange();
                octbuild(&cellsinfo,celltype);
                createexportlist(system,settings,cellsinfo,celltype,&commdata);
                singlestep(system,&cellsinfo,celltype,&commdata);
                exchangecleanup(system,cellsinfo,&commdata);
                printstatistics(system,settings,cellsinfo,&statistics);
                cellsupdate(settings,&cellsinfo);
                writevtk(system,settings,cellsinfo);
                octfree(&cellsinfo);
        }

        MPI_Abort(MPI_COMM_WORLD,-1);


        for (step = 0; step < nsteps; step++) {

//    ioWriteStepVTK(step);

                if (!(step % statOutStep))
                        printStepNum();

                //decompositionExecute();
                //octBuild();
                //createExportList();
                //computeStep();

                //if (!(step % statOutStep))
                //statisticsPrint();

                if (simStart)
                        simTime += secondsPerStep / 3600.0; /* biological process time in hours */

                if (!(step % vtkOutStep)) {
                        if (vtkout)
                                ioWriteStepVTK(step);
                        if (povout)
                                ioWriteStepPovRay(step, 0);
//      if (vnfout)
//        ioWriteFields(step);
                }

                //updateCellPositions(statistics);
                //updateCellStates();
                //commCleanup();
                octfree();

                //if (!(step % rstOutStep))
                //        saveRstFile();
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //decompositionFinalize();
        randomstreamfree(&settings);
        cellsCleanup();

        if (MPIrank == 0)
                printf("\nEnd of simulation run.\n");

        MPI_Finalize();

        return 0;
}
