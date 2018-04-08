

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
        cellcommdata_t cellcommdata;
        fieldcommdata_t fieldcommdata;
        interpdata_t interpdata;
        statistics_t statistics;
        solverdata_t solverdata;
        cellenvdata_t **cellenvdata;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &systeminfo.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &systeminfo.rank);
        systeminfo.nthreads = omp_get_max_threads();

        getsysteminfo(&systeminfo);
        printinfo(systeminfo);
        initialisation(argc,argv,&systeminfo,&settings,&celltype,&environment);
        allocatecells(systeminfo,settings,celltype,&cellsinfo);
        allocategrid(systeminfo,settings,&grid);
        //envalloc(systeminfo,settings,grid,&environment,&solverdata);
        //envinit(systeminfo,settings,grid,&environment);
        allocatefields(systeminfo,settings,grid,&environment,&solverdata);
        initfields(systeminfo,settings,grid,&environment);
        lbinit(argc,argv,MPI_COMM_WORLD,systeminfo,&cellsinfo);

        for (settings.step = 0; settings.step < settings.numberofsteps; settings.step++) {
                updateglobalcounts(&cellsinfo);
                lbexchange(systeminfo);
                octbuild(systeminfo,&cellsinfo,celltype);
                createexportlist(systeminfo,settings,cellsinfo,grid,celltype,&cellcommdata,&fieldcommdata);
                singlestep(systeminfo,settings,&cellsinfo,celltype,&grid,&environment,&cellcommdata,&interpdata,&cellenvdata);
                exchangecleanup(systeminfo,cellsinfo,&cellcommdata,&fieldcommdata);
                printstatistics(systeminfo,settings,cellsinfo,&statistics);
                cellsupdate(settings,&cellsinfo);
                writevtk(systeminfo,settings,cellsinfo);
                octfree(&cellsinfo);
                cleanstep(settings,&cellenvdata);
        }

        MPI_Abort(MPI_COMM_WORLD,-1);

        MPI_Barrier(MPI_COMM_WORLD);

        //decompositionFinalize();
        cellsdestroy();
        lbdestroy();

        if (MPIrank == 0)
                printf("\nEnd of simulation run.\n");

        MPI_Finalize();

        return 0;
}
