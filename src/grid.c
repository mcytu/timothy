/* **************************************************************************
 * Timothy - Tissue Modelling Framework
 * Copyright (C) 2014-2018 Maciej Cytowski
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
#include <inttypes.h>

#include "global.h"
//#include "fields.h"

/*! \file grid.c
 *  \brief contains functions that build the computational grid
 */

 /*!
  * This function computes the sizes of the grid.
  */
 void computegridsize(system_t system,settings_t settings,grid_t *grid)
 {
   float3dv_t globalgridsize;
   grid->resolution=settings.gfh;

   grid->lowercorner.x = -BOXSIZEX/2.0;
   grid->lowercorner.y = -BOXSIZEY/2.0;
   if(settings.dimension==3) grid->lowercorner.z = -BOXSIZEZ/2.0;
   else grid->lowercorner.z = 0.0;

   grid->uppercorner.x = (BOXSIZEX/2.0);//-1.0;
   grid->uppercorner.y = (BOXSIZEY/2.0);//-1.0;
   if(settings.dimension==3) grid->uppercorner.z = (BOXSIZEZ/2.0);//-1.0;
   else grid->uppercorner.z = 0.0;

   globalgridsize.x = grid->uppercorner.x - grid->lowercorner.x + 1;
   globalgridsize.y = grid->uppercorner.y - grid->lowercorner.y + 1;
   if(settings.dimension==3) globalgridsize.z = grid->uppercorner.z - grid->lowercorner.z + 1;
   else globalgridsize.z = 0.0;

   grid->globalsize.x = (int64_t)((globalgridsize.x+1)/grid->resolution);
   grid->globalsize.y = (int64_t)((globalgridsize.y+1)/grid->resolution);
   if(settings.dimension==3) grid->globalsize.z = (int64_t)((globalgridsize.z+1)/grid->resolution);
   else grid->globalsize.z = 0.0;

   grid->globalsize.x = grid->globalsize.x + (system.dim[0] - grid->globalsize.x % system.dim[0]);
   grid->globalsize.y = grid->globalsize.y + (system.dim[1] - grid->globalsize.y % system.dim[1]);
   if(settings.dimension==3) grid->globalsize.z = grid->globalsize.z + (system.dim[2] - grid->globalsize.z % system.dim[2]);

   grid->localsize.x = grid->globalsize.x / system.dim[0];
   grid->localsize.y = grid->globalsize.y / system.dim[1];
   if(settings.dimension==3) grid->localsize.z = grid->globalsize.z / system.dim[2];
   else grid->localsize.z = 1;

   if(system.rank==0) {
     printf("environment grid size = %" PRId64 "x%" PRId64 "x%" PRId64 "\n", grid->globalsize.x, grid->globalsize.y,grid->globalsize.z);
     fflush(stdout);
   }

   return;
 }


 void allocategrid(system_t system,settings_t settings,grid_t *grid) {
   int i,j,k;
   #define gridnode(i,j,k) (grid->data[grid->localsize.y*grid->localsize.z*i+grid->localsize.z*j+k])

   computegridsize(system,settings,grid);

   if(!(grid->loweridx=(int643dv_t*)malloc(system.size*sizeof(int643dv_t))))
     terminate(system,"cannot allocate grid->loweridx", __FILE__, __LINE__);
   if(!(grid->upperidx=(int643dv_t*)malloc(system.size*sizeof(int643dv_t))))
     terminate(system,"cannot allocate grid->upperidx", __FILE__, __LINE__);

   for(i=0;i<system.size;i++) {
     grid->loweridx[i].x = grid->localsize.x * system.coords[i][0];
     grid->loweridx[i].y = grid->localsize.y * system.coords[i][1];
     if(settings.dimension==3) grid->loweridx[i].z = grid->localsize.z * system.coords[i][2];
     else grid->loweridx[i].z = 0;
     grid->upperidx[i].x = grid->loweridx[i].x + grid->localsize.x - 1;
     grid->upperidx[i].y = grid->loweridx[i].y + grid->localsize.y - 1;
     if(settings.dimension==3) grid->upperidx[i].z = grid->loweridx[i].z + grid->localsize.z - 1;
     else grid->upperidx[i].z = 0;
   }
   if(!(grid->data=(double3dv_t*)calloc(grid->localsize.x*grid->localsize.y*grid->localsize.z,sizeof(double3dv_t))))
     terminate(system,"cannot allocate grid->data", __FILE__, __LINE__);

   for(i=0;i<grid->localsize.x;i++)
     for(j=0;j<grid->localsize.y;j++)
       for(k=0;k<grid->localsize.z;k++) {
         gridnode(i,j,k).x = grid->lowercorner.x + grid->resolution * (grid->loweridx[system.rank].x + i);
         gridnode(i,j,k).y = grid->lowercorner.y + grid->resolution * (grid->loweridx[system.rank].y + j);
         gridnode(i,j,k).z = grid->lowercorner.z + grid->resolution * (grid->loweridx[system.rank].z + k);
       }

   #undef grid_node
   return;
 }
