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
int computeStep()
{
  int p;
  double dvel,sf;
#ifdef __MIC__
  double *phiwork1,*phiwork2;
#endif

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

  initCellsToGridExchange();

  /* initiate asynchronous data transfers between processors */
  cellsExchangeInit();

  /* 1. Compute potential for local cells */

  /* offload trasfer: host -> accelerator  */
#ifdef __MIC__
#pragma offload_transfer target(mic:MPIrank%2) in(cells:length(lnc) alloc_if(1) free_if(0)) in(octree:length(octSize) alloc_if(1) free_if(0)) nocopy(velocity:length(lnc) alloc_if(1) free_if(0)) in(statistics)
#endif
  /* compute potential for local cells */
#ifdef __MIC__
#pragma offload target(mic:MPIrank%2) signal(phiwork1) nocopy(cells:length(lnc) alloc_if(0) free_if(0)) nocopy(octree:length(octSize) alloc_if(0) free_if(0)) in(affShift,affScale,h,tnc,lnc,nx,sdim,h2,h3,csize)
#endif
  compPot();

  /* 2. Solve global fields */

  if(step>0) {
    waitCellsToGridExchange();
    fieldsSolve();
  }
#ifdef __MIC__
#pragma offload_wait target(mic:MPIrank%2) wait(phiwork1)
#endif
  /* wait for data transfers to finish */
  cellsExchangeWait();

  /* 3. Compute potential for remote cells */

#ifdef __MIC__
#pragma offload target(mic:MPIrank%2) nocopy(cells:length(lnc) alloc_if(0) free_if(0)) nocopy(octree:length(octSize) alloc_if(0) free_if(0)) in(recvData:length(numImp) alloc_if(1) free_if(0))  in(affShift,affScale,h,tnc,lnc,nx,sdim,h2,h3,csize,numImp)
#endif
  compRPot();

  /* 4. Add chemotactic term to potential */

  /* add chemotaxis term to potential */
/*  if(bvsim) {
    int p;
    for(p=0; p<lnc; p++)
      cells[p].v+=100*sqrt(pow(cellFields[NFIELDS][p],2)+pow(cellFields[NFIELDS+1][p],2)+pow(cellFields[NFIELDS+2][p],2));
  }*/

  /* 5. Compute gradient of the potential for local cells */

#ifdef __MIC__
#pragma offload_transfer target(mic:MPIrank%2) out(cells:length(lnc) alloc_if(0) free_if(0)) out(statistics)
#endif
  /* initiate transfer of the density and potential data from remote cells */
  densPotExchangeInit();
  /* compute gradient of the potential for local cells */
#ifdef __MIC__
#pragma offload target(mic:MPIrank%2) signal(phiwork2) nocopy(cells:length(lnc) alloc_if(0) free_if(0)) nocopy(octree:length(octSize) alloc_if(0) free_if(0)) in(affShift,affScale,h,tnc,lnc,nx,sdim,h2,h3,h4,csize) nocopy(velocity:length(lnc) alloc_if(0) free_if(0))
#endif
  compPotGrad();

  /* 6. Interpolate global fields and compute gradient */

  /* interpolate data */
  interpolateFieldsToCells();
  /* compute gradient of global fields */
  fieldGradient();
#ifdef __MIC__
#pragma offload_wait target(mic:MPIrank%2) wait(phiwork2)
#endif

  /* 7. Compute gradient of the potential for remote cells */
  /* wait for density and potential data from remote cells */
  densPotExchangeWait();
#ifdef __MIC__
#pragma offload_transfer target(mic:MPIrank%2) out(cells:length(lnc) alloc_if(0) free_if(1)) nocopy(octree:length(octSize) alloc_if(0) free_if(1)) out(velocity:length(lnc) alloc_if(0) free_if(1))
#endif
#ifdef __MIC__
#pragma offload_transfer target(mic:MPIrank%2) nocopy(recvData:length(numImp) alloc_if(0) free_if(1))
#endif
  /* compute gradient of the potential for remote cells */
  compRPotGrad();

  /* 8. Correct velocity for various cell types */
  for(p=0; p<lnc; p++)
    if(cells[p].ctype!=1) {
      velocity[p].x += 0.01*cellFields[NFIELDS][p];
      velocity[p].y += 0.01*cellFields[NFIELDS+1][p];
      velocity[p].z += 0.01*cellFields[NFIELDS+2][p];
    }

  /* 9. Compute and collect statistical data */
  for (p = 0; p < lnc; p++) {
    dvel =
      sqrt(velocity[p].x * velocity[p].x +
           velocity[p].y * velocity[p].y +
           velocity[p].z * velocity[p].z);
    if (dvel < statistics.minvel)
      statistics.minvel = dvel;
    if (dvel > statistics.maxvel)
      statistics.maxvel = dvel;
    if (cells[p].size < statistics.minsize)
      statistics.minsize = cells[p].size;
    if (cells[p].size > statistics.maxsize)
      statistics.maxsize = cells[p].size;
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

  //printf("globalMaxVel=%f sf=%.16f\n",globalMaxVel,sf);

  printf("maxspeed=%f\n", maxSpeedInUnits * secondsPerStep);

  statistics.minvel = DBL_MAX;  /* minimal velocity is set to DBL_MAX */
  statistics.maxvel = 0;        /* maximal velocity is set to zero */

  for (p = 0; p < lnc; p++) {
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

}
