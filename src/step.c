/* **************************************************************************
 * This file is part of Timothy
 *
 * Copyright (c) 2014/15 Maciej Cytowski
 * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "potential.h"
#include "fields.h"
//#include "inline.h"
#include "exchange.h"

/*! \file compute.c
 *  \brief contains main computational function called in each time step of the simulation
 */

/*!
 * This function calls all important simulation steps (cellular dynamics and global fields computations).
 */
int singlestep(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, celltype_t* celltype, grid_t *grid, environment_t **environment,cellcommdata_t *cellcommdata,interpdata_t *interpdata,cellenvdata_t ***cellenvdata, solverdata_t *solverdata,solversettings_t *solversettings)
{
        int p;
        double sf;
        patches_t patches;
        /* 0. Initialization */
        patchesalloc(systeminfo,settings,&patches,cellsinfo,grid);
        cells2envinit(systeminfo,settings,&patches,cellsinfo,grid);

        /* initiate asynchronous data transfers between processors */
        cellssendrecv(systeminfo,*cellsinfo,cellcommdata);

        /* 1. Compute potential for local cells */

        /* compute potential for local cells */
        computepotential(cellsinfo,celltype,*cellcommdata);

        /* 2. Solve global fields */

        if(settings.step>0) {
                cells2envwait(systeminfo,settings,&patches,grid,environment);
                envcompute(systeminfo,settings,grid,environment,solverdata,solversettings);
        }

        /* wait for data transfers to finish */
        cellswait(systeminfo,*cellsinfo,cellcommdata);

        /* 3. Compute potential for remote cells */
        computeremotepotential(cellsinfo,celltype,*cellcommdata);

        /* 4. Add chemotactic term to potential */

        /* add chemotaxis term to potential */

        /* 5. Compute gradient of the potential for local cells */
        /* initiate transfer of the density and potential data from remote cells */
        datasendrecv(systeminfo,*cellsinfo,cellcommdata);
        /* compute gradient of the potential for local cells */
        computegradient(cellsinfo,celltype,*cellcommdata);

        /* 6. Interpolate global fields and compute gradient */

        /* interpolate data */
        env2cellsinit(systeminfo,settings,&patches,grid,environment);

        env2cellswait(systeminfo,settings,&patches,cellsinfo,grid,cellenvdata);

        patchesfree(&patches);
        /* compute gradient of global fields */
        //fieldGradient();

        /* 7. Compute gradient of the potential for remote cells */
        /* wait for density and potential data from remote cells */
        datawait(systeminfo,*cellsinfo,cellcommdata);

        /* compute gradient of the potential for remote cells */
        computeremotegradient(cellsinfo,celltype,*cellcommdata);

        return 0;
}

void cleanstep(settings_t settings,cellenvdata_t ***cellenvdata) {
        int f;
        for(f=0; f<settings.numberoffields; f++) {
                free((*cellenvdata)[f]);
        }
        free((*cellenvdata));
}
