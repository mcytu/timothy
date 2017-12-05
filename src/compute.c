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
#include <float.h>
#include <math.h>

#include "global.h"
#include "potential.h"
#include "fields.h"
#include "inline.h"

/*! \file compute.c
 *  \brief contains main computational function called in each time step of the simulation
 */

/*!
 * This function calls all important simulation steps (cellular dynamics and global fields computations).
 */
int computestep(system_t system, cellsinfo_t *cellsinfo, commdata_t *commdata)
{
        int p;
        double dvel,sf;

        /* 0. Initialization */

        /* initialize statistics data */
        if(nc>1)
                statistics.mindist=DBL_MAX;
        else
                statistics.mindist=0.0;
        statistics.minvel=DBL_MAX;
        statistics.maxvel=0.0;
        statistics.minsize=DBL_MAX;
        statistics.maxsize=0.0;

        //initCellsToGridExchange();

        /* initiate asynchronous data transfers between processors */
        cellscomminit(system,*cellsinfo,commdata);

        /* 1. Compute potential for local cells */

        /* compute potential for local cells */
        computepotential(cellsinfo,*commdata);

        /* 2. Solve global fields */

        //if(step>0) {
        //        waitCellsToGridExchange();
        //        fieldsSolve();
        //}


        /* wait for data transfers to finish */
        cellscommwait(system,*cellsinfo,commdata);

        /* 3. Compute potential for remote cells */
        computeremotepotential(cellsinfo,*commdata);

        /* 4. Add chemotactic term to potential */

        /* add chemotaxis term to potential */

        /* 5. Compute gradient of the potential for local cells */
        /* initiate transfer of the density and potential data from remote cells */
        dpcomminit(system,*cellsinfo,commdata);
        /* compute gradient of the potential for local cells */
        computegradient(cellsinfo,*commdata);

        /* 6. Interpolate global fields and compute gradient */

        /* interpolate data */
        //interpolateFieldsToCells();
        /* compute gradient of global fields */
        //fieldGradient();

        /* 7. Compute gradient of the potential for remote cells */
        /* wait for density and potential data from remote cells */
        dpcommwait(system,*cellsinfo,commdata);

        /* compute gradient of the potential for remote cells */
        computeremotegradient(cellsinfo,*commdata);

        /* 8. Correct forces for various cell types */
/*        for(p=0; p<lnc; p++)
                if(cells[p].ctype!=1) {
                        cellsinfo->forces[p].x += 0.01*cellFields[NFIELDS][p];
                        cellsinfo->forces[p].y += 0.01*cellFields[NFIELDS+1][p];
                        cellsinfo->forces[p].z += 0.01*cellFields[NFIELDS+2][p];
                }
 */
        /* 9. Compute and collect statistical data */
        for (p = 0; p < cellsinfo->localcount.n; p++) {
                dvel =
                        sqrt(cellsinfo->forces[p].x * cellsinfo->forces[p].x +
                             cellsinfo->forces[p].y * cellsinfo->forces[p].y +
                             cellsinfo->forces[p].z * cellsinfo->forces[p].z);
                printf("dvel=%f\n",dvel);
                if (dvel < statistics.minvel)
                        statistics.minvel = dvel;
                if (dvel > statistics.maxvel)
                        statistics.maxvel = dvel;
                if (cellsinfo->cells[p].size < statistics.minsize)
                        statistics.minsize = cellsinfo->cells[p].size;
                if (cellsinfo->cells[p].size > statistics.maxsize)
                        statistics.maxsize = cellsinfo->cells[p].size;
        }

        /* this should be removed soon (do we really need to reduceall here?) */
        MPI_Allreduce(&statistics.minvel, &globalMinVel, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&statistics.maxvel, &globalMaxVel, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);

        if (globalMaxVel == 0.0)
                sf = 0.0;
        else
                sf =  maxSpeedInUnits * secondsPerStep / globalMaxVel;

        //printf("sf=%f\n",sf);
        //sf=0.001;

        printf("globalMaxVel=%f sf=%.16f\n",globalMaxVel,sf);

        //printf("maxspeed=%f\n", maxSpeedInUnits * secondsPerStep);

//i        statistics.minvel = DBL_MAX; /* minimal velocity is set to DBL_MAX */
//i        statistics.maxvel = 0;  /* maximal velocity is set to zero */

/*i        for (p = 0; p < lnc; p++) {
   //    velocity[p].x *= sf;
   //    velocity[p].y *= sf;
   //    velocity[p].z *= sf;
                dvel =
                        sqrt(velocity[p].x * velocity[p].x +
                             velocity[p].y * velocity[p].y +
                             velocity[p].z * velocity[p].z);
                if(dvel > maxSpeedInUnits*secondsPerStep) {
                        velocity[p].x *= maxSpeedInUnits*secondsPerStep/dvel;
                        velocity[p].y *= maxSpeedInUnits*secondsPerStep/dvel;
                        velocity[p].z *= maxSpeedInUnits*secondsPerStep/dvel;
                }
                dvel =
                        sqrt(velocity[p].x * velocity[p].x +
                             velocity[p].y * velocity[p].y +
                             velocity[p].z * velocity[p].z);

                if (dvel < statistics.minvel)
                        statistics.minvel = dvel;
                if (dvel > statistics.maxvel)
                        statistics.maxvel = dvel;
        }
 */
}
