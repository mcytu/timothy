

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#include "global.h"

#include "lb.h"
#include "io.h"
#include "cells.h"
#include "grid.h"
#include "step.h"
#include "fields.h"
#include "utils.h"
#include "octree.h"
#include "exchange.h"
#include "initialisation.h"

/*! \file main.c
 *  \brief contains the main simulation loop
 */

/*!
 * This function intializes MPI, calls Timothy initialization and allocation functions.
 * It also contains the main simulation loop where all important simulation steps are called.
 */

int main(int argc, char **argv)
{
        systeminfo_t systeminfo;
        settings_t settings;
        celltype_t *celltype;
        environment_t *environment;
        cellsinfo_t cellsinfo;
        grid_t grid;
        commdata_t commdata;
        statistics_t statistics;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &systeminfo.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &systeminfo.rank);
        systeminfo.nthreads = omp_get_max_threads();

        getsysteminfo(&systeminfo);
        printinfo(systeminfo);
        initialisation(argc,argv,&systeminfo,&settings,&celltype,&environment);
        allocatecells(systeminfo,settings,celltype,&cellsinfo);
        allocategrid(systeminfo,settings,&grid);
        allocatefields(systeminfo,settings,grid,&environment);
        lbinit(argc,argv,MPI_COMM_WORLD,systeminfo,&cellsinfo);

        for (settings.step = 0; settings.step < settings.numberofsteps; settings.step++) {
                updateglobalcounts(&cellsinfo);
                lbexchange(systeminfo);
                octbuild(systeminfo,&cellsinfo,celltype);
                createexportlist(systeminfo,settings,cellsinfo,celltype,&commdata);
                singlestep(systeminfo,&cellsinfo,celltype,&commdata);
                exchangecleanup(systeminfo,cellsinfo,&commdata);
                printstatistics(systeminfo,settings,cellsinfo,&statistics);
                cellsupdate(settings,&cellsinfo);
                writevtk(systeminfo,settings,cellsinfo);
                octfree(&cellsinfo);
        }

        MPI_Abort(MPI_COMM_WORLD,-1);


        for (step = 0; step < nsteps; step++) {

//    ioWriteStepVTK(step);

                //  if (!(step % statOutStep))
                //          printStepNum();

                //decompositionExecute();
                //octBuild();
                //createExportList();
                //computeStep();

                //if (!(step % statOutStep))
                //statisticsPrint();

                if (simStart)
                        simTime += secondsPerStep / 3600.0; /* biological process time in hours */

                /*      if (!(step % vtkOutStep)) {
                              if (vtkout)
                                      ioWriteStepVTK(step);
                              if (povout)
                                      ioWriteStepPovRay(step, 0);
                   //      if (vnfout)
                   //        ioWriteFields(step);
                      }*/

                //updateCellPositions(statistics);
                //updateCellStates();
                //commCleanup();
                //octfree();

                //if (!(step % rstOutStep))
                //        saveRstFile();
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //decompositionFinalize();
        cellsdestroy();
        lbdestroy();

        if (MPIrank == 0)
                printf("\nEnd of simulation run.\n");

        MPI_Finalize();

        return 0;
}
