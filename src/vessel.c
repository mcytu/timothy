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

/*! \file vessel.c
 *  \brief contains functions defining virtual vessels
 */

/*!
 * This functions builds a tube vessel from available data.
 */
void initVessel()
{
  int i;

  double minx,miny,minz;
  double maxx,maxy,maxz;

  double dist;
  double mindist;
  double middle[3];

  vc=0;
  lvc=0;

  lnc=0;

  //printf("csize=%f\n",csize);

  //csize=0.5;
  h=2*csize;
  if (MPIrank == 0) {
    FILE *fh;
    double ax,ay,az;
    if(!(fh=fopen("vessel.txt","r"))) {
      printf("File vessel.txt not found\n");
      exit(1);
    }
    minx=DBL_MAX;
    miny=DBL_MAX;
    minz=DBL_MAX;
    maxx=-DBL_MAX;
    maxy=-DBL_MAX;
    maxz=-DBL_MAX;
    while(fscanf(fh,"%lf %lf %lf\n",&ax,&ay,&az)!=EOF) {
      int skip=0;
      for(i=0; i<lnc; i++) {
        dist=sqrt(pow(cells[i].x-ax,2.0)+pow(cells[i].y-ay,2.0)+pow(cells[i].z-az,2.0));
        if(dist<0.01) {
          skip=1;
          break;
        }
      }
      if(skip) continue; 
      cells[lnc].size=csize;//csize;
      cells[lnc].x=ax;
      cells[lnc].y=ay;
      cells[lnc].z=az;
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
      cells[lnc].h = h;//h;
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
      /* tag vessel cell */
      cells[lnc].ctype = 1;

      lnc++;
      lvc++;
      localID++;
    }
    /* scale */
    for(i=0; i<lnc; i++) {
      cells[i].x-=minx;
      cells[i].y-=miny;
      cells[i].z-=minz;
    }
    /* find a single cell in the middle */
    middle[0]=0.0;
    middle[1]=0.0;
    middle[2]=0.0;
    for(i=0; i<lnc; i++) {
      middle[0]+=cells[i].x;
      middle[1]+=cells[i].y;
      middle[2]+=cells[i].z;
    }
    middle[0]/=lnc;
    middle[1]/=lnc;
    middle[2]/=lnc;
    mindist=DBL_MAX;
    for(i=0; i<lnc; i++) {
      dist=sqrt(pow(cells[i].x-middle[0],2.0)+pow(cells[i].y-middle[1],2.0)+pow(cells[i].z-middle[2],2.0));
      if(dist<mindist) {
        middleCellIdx=i;
        mindist=dist;
      }
    }
    /* compare with box size */

    lnc=lnc+1;

    fclose(fh);
  }
  
  //csize=0.5;
  //h=3*csize;
  //printf("h=%f\n",h);
  //if(lnc==0) lnc=1;

  

  MPI_Allreduce(localCellCount, totalCellCount, numberOfCounts,
                MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  printf("lnc=%ld nc=%ld lvc=%ld vc=%ld\n",lnc,nc,lvc,vc);
}
