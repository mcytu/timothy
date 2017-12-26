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
#include <mpi.h>
#include <math.h>
#include <inttypes.h>
#include <float.h>


#include "global.h"

/*! \file stats.c
 *  \brief contains functions computing and printing simulation statisctical information
 */

void printstatistics(system_t system,settings_t settings,cellsinfo_t cellsinfo,statistics_t* statistics)
{
        int p;
        statistics->mindist=DBL_MAX;
        statistics->minspeed=DBL_MAX;
        statistics->maxspeed=0.0;
        statistics->minsize=DBL_MAX;
        statistics->maxsize=0.0;
        statistics->mindens=DBL_MAX;
        statistics->maxdens=0.0;
        for (p = 0; p < cellsinfo.localcount.n; p++) {
                double speed;
                if(cellsinfo.cells[p].size==0) printf("ZERO?? %d\n",p);
                speed =  sqrt(cellsinfo.forces[p].x * cellsinfo.forces[p].x +
                              cellsinfo.forces[p].y * cellsinfo.forces[p].y +
                              cellsinfo.forces[p].z * cellsinfo.forces[p].z);
                statistics->minspeed = (speed < statistics->minspeed ? speed : statistics->minspeed);
                statistics->maxspeed = (speed > statistics->maxspeed ? speed : statistics->maxspeed);
                statistics->minsize = (cellsinfo.cells[p].size < statistics->minsize ? cellsinfo.cells[p].size : statistics->minsize);
                statistics->maxsize = (cellsinfo.cells[p].size > statistics->maxsize ? cellsinfo.cells[p].size : statistics->maxsize);
                statistics->mindist = (cellsinfo.cells[p].mindist < statistics->mindist ? cellsinfo.cells[p].mindist : statistics->mindist);
                statistics->mindens = (cellsinfo.cells[p].density < statistics->mindens ? cellsinfo.cells[p].density : statistics->mindens);
                statistics->maxdens = (cellsinfo.cells[p].density > statistics->maxdens ? cellsinfo.cells[p].density : statistics->maxdens);
        }

        MPI_Allreduce(MPI_IN_PLACE,&(statistics->minspeed), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->maxspeed), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->minsize), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->maxsize), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->mindist), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        if(system.rank==0) {
                printf("\n+++ simulation step %12d\n",settings.step);
                printf("%12s%10s%10s\n", "", "min", "max");
                printf("%12s%10.4lf%10.4lf\n", "size       ",statistics->minsize,statistics->maxsize);
                printf("%12s%10.4lf%10.4lf\n", "density    ",statistics->mindens,statistics->maxdens);
                printf("%12s%10.4lf%10.4lf\n", "speed      ",statistics->minspeed,statistics->maxspeed);
                printf("%12s%10.4lf%10s\n", "distance   ", statistics->mindist, "N/A");
        }

        return;
}
