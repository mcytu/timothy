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

void computestatistics(cellsinfo_t cellsinfo,statistics_t* statistics)
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
        return;
}

void printstatistics(system_t system,statistics_t statistics)
{
        if(system.rank==0) {
                printf("%12s%10s%10s\n", "", "Min", "Max");
                printf("%12s%10.4lf%10.4lf\n", "Size       ",statistics.minsize,statistics.maxsize);
                printf("%12s%10.4lf%10.4lf\n", "Density    ",statistics.mindens,statistics.maxdens);
                printf("%12s%10.4lf%10.4lf\n", "Speed      ",statistics.minspeed,statistics.maxspeed);
                printf("%12s%10.4lf%10s\n", "Distance   ", statistics.mindist, "N/A");
        }
        return;
}

/*!
 * This function computes and prints out phases statistics.
 */
/*void statisticsPhases()
   {
        if (MPIrank == 0) {
                printf("%12s%12s%12s%12s%12s%12s\n", "Cell phase ", "G0", "G1", "S",
                       "G2", "M");
                printf("%12s%12" PRId64 "%12" PRId64 "%12" PRId64 "%12" PRId64 "%12"
                       PRId64 "\n", "N. of cells", g0nc, g1nc, snc, g2nc, mnc);
                printf("\n%16s%12" PRId64 "\n", "Healthy cells  ", nc - cnc - nnc);
                printf("%16s%12" PRId64 "\n", "Cancer cells   ", cnc);
                printf("%16s%12" PRId64 "\n", "Necrotic cells ", nnc);
                if(scsim) printf("%16s%12" PRId64 "\n", "Blood cells ", globalbc);
        }
   }*/
