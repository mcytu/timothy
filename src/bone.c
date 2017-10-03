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
#include <sprng.h>
#include <float.h>

#include "global.h"

/*! \file bone.c
 *  \brief contains functions defining virtual bone structure
 */

/*!
 * This functions builds a bone structure from available data.
 */
int initBone()
{
  int i,j,k;
  int m;
  int ax,ay,az;
  double ox,oy,oz;
  double v0x,v0y,v0z;
  double v1x,v1y,v1z;
  double v2x,v2y,v2z;
  char text[256];
  int ***data;
  int d1,d2;
  int ret;

  double minx,miny,minz;
  double maxx,maxy,maxz;

  bnc=0;
  lbnc=0;

  if(MPIrank==0) {
    FILE *fh1,*fh2;
    if(!(fh1=fopen("bones.txt","r"))) {
      printf("File bones.txt not found\n");
      exit(1);
    }
    if(!(fh2=fopen("bones.vnf","r"))) {
      printf("File bones.vnf not found\n");
      exit(1);
    }
    minx=DBL_MAX;
    miny=DBL_MAX;
    minz=DBL_MAX;
    maxx=-DBL_MAX;
    maxy=-DBL_MAX;
    maxz=-DBL_MAX;

    /* start with the .vnf file */
    /* skip the header */
    ret=fscanf(fh2,"%*[^\n]\n");
    /* read dimensions */
    ret=fscanf(fh2,"%s",text);
    if(strcmp(text,"dims")==0) {
      ret=fscanf(fh2,"%d %d %d\n",&ax,&ay,&az);
    }
    /* read 3d box */
    ret=fscanf(fh2,"%s",text);
    if(strcmp(text,"origin")==0) {
      ret=fscanf(fh2,"%lf %lf %lf\n",&ox,&oy,&oz);
    }
    ret=fscanf(fh2,"%s",text);
    if(strcmp(text,"v0")==0) {
      ret=fscanf(fh2,"%lf %lf %lf\n",&v0x,&v0y,&v0z);
    }
    ret=fscanf(fh2,"%s",text);
    if(strcmp(text,"v1")==0) {
      ret=fscanf(fh2,"%lf %lf %lf\n",&v1x,&v1y,&v1z);
    }
    ret=fscanf(fh2,"%s",text);
    if(strcmp(text,"v2")==0) {
      ret=fscanf(fh2,"%lf %lf %lf\n",&v2x,&v2y,&v2z);
    }
    /* print info */
    printf("Dimensions: %d %d %d\nOrigin: %lf %lf %lf\n",ax,ay,az,ox,oy,oz);
    printf("v0: %lf %lf %lf\nv1: %lf %lf %lf\nv2: %lf %lf %lf\n",v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z);
    /* allocate data table */
    data=(int***)malloc(ax*sizeof(int**));
    for(i=0; i<ax; i++) {
      data[i]=(int**)malloc(ay*sizeof(int*));
      for(j=0; j<ay; j++) {
        data[i][j]=(int*)malloc(az*sizeof(int));
      }
    }

    /* now .txt file */
    /* skip the header */
    ret=fscanf(fh1,"%*[^\n]\n");
    /* read data */
    for(k=0; k<az; k++)
      for(j=0; j<ay; j++)
        for(i=0; i<ax; i++) {
          ret=fscanf(fh1,"%d %d\n",&d1,&d2);
          data[i][j][k]=d1;
        }

    /* compute coordinates - use origin,v0,v1,v2 */
    for(i=0; i<ax; i++)
      for(j=0; j<ay; j++)
        for(k=0; k<az; k++) {
          if (data[i][j][k]==2) {
            cells[lnc].size=pow(2.0, -(1.0 / 3.0)) * csize;//csize;
            cells[lnc].x=ox + i*v0x + j*v1x + k*v2x;
            cells[lnc].y=oy + i*v0y + j*v1y + k*v2y;
            cells[lnc].z=oz + i*v0z + j*v1z + k*v2z;
            maxx=(cells[lnc].x>maxx?cells[lnc].x:maxx);
            maxy=(cells[lnc].y>maxy?cells[lnc].y:maxy);
            maxz=(cells[lnc].z>maxz?cells[lnc].z:maxz);
            minx=(cells[lnc].x<minx?cells[lnc].x:minx);
            miny=(cells[lnc].y<miny?cells[lnc].y:miny);
            minz=(cells[lnc].z<minz?cells[lnc].z:minz);

            cells[lnc].gid =
              (unsigned long long int) MPIrank *(unsigned long long int)
              maxCellsPerProc + (unsigned long long int) lnc;
            cells[lnc].v = 0.0;
            cells[lnc].density = 0.0;
            cells[lnc].h = 3*csize;//h;
            cells[lnc].young = 2100.0 + sprng(stream) * 100.0;
            cells[lnc].halo = 0;
            cells[lnc].phase = 0;
            cells[lnc].g1 = g1 * (1 + (sprng(stream) * 2 - 1) * v);
            cells[lnc].g2 = g2 * (1 + (sprng(stream) * 2 - 1) * v);
            cells[lnc].s = s * (1 + (sprng(stream) * 2 - 1) * v);
            cells[lnc].m = m * (1 + (sprng(stream) * 2 - 1) * v);
            cells[lnc].phasetime = 0.0;
            cells[lnc].age = 0;
            cells[lnc].death = 0;
            cells[lnc].tumor = 0;
            /* tag bone cell */
            cells[lnc].ctype = 3;

            lnc++;
            lbnc++;
            localID++;

          }
        }
    fclose(fh1);
    fclose(fh2);
  }

  for(i=0; i<lnc; i++) {
    cells[i].x-=minx;
    cells[i].y-=miny;
    cells[i].z-=minz;
  }

  MPI_Allreduce(localCellCount, totalCellCount, numberOfCounts,
                MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  /* free data */
  for(i=0; i<ax; i++) {
    for(j=0; j<ay; j++)
      free(data[i][j]);
    free(data[i]);
  }
  free(data);

  return 0;
}

