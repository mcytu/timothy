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
#include "fields.h"

/*! \file
 *  \brief contains driving functions for the global fields
 */

/* Global fields' IDs:
 * 0 - density field
 * 1 - temperature
 * 2 - ??
 * 3 - ??
 */

/*!
 * This function intializes the fields.
 * Various parameters are set and the Hypre initialization of the system is executed.
 */
void fieldsInit()
{
  int nf, i;
  int chf;

  if (!gfields)
    return;

  nf = 0;

  strcpy(fieldName[nf], "tissue");
  fieldType[nf] = SCALAR_FIELD;
  fieldAddr[nf] =
    (double *) calloc(gridSize.x * gridSize.y * gridSize.z,
                      sizeof(double));
  tissueField = (double *) fieldAddr[nf];
  for (i = 0; i < gridSize.x * gridSize.y * gridSize.z; i++)
    tissueField[i] = 0.0;
  nf++;

  if(bvsim) {
    strcpy(fieldName[nf], "vessel");
    fieldType[nf] = SCALAR_FIELD;
    fieldAddr[nf] =
      (double *) calloc(gridSize.x * gridSize.y * gridSize.z,
                        sizeof(double));
    vesselField = (double *) fieldAddr[nf];
    for (i = 0; i < gridSize.x * gridSize.y * gridSize.z; i++)
      vesselField[i] = 0.0;
  }
  nf++;

  if(temperature) {
    strcpy(fieldName[nf], "temp");
    fieldType[nf] = SCALAR_FIELD;
    fieldAddr[nf] =
      (double *) calloc(gridSize.x * gridSize.y * gridSize.z,
                        sizeof(double));
    tempField = (double *) fieldAddr[nf];
    for (i = 0; i < gridSize.x * gridSize.y * gridSize.z; i++)
      tempField[i] = 0.0;
    fieldDiffCoef[nf] = 1 * 1e-5;
    fieldBC[nf] = 42.0;
    fieldICMean[nf] = 36.0;
    fieldICVar[nf] = 0;
  }
  nf++;

  /* set default values for chemical parameters */
  for (chf = 0; chf < NCHEM; chf++) {
    if(chf==0 && !oxygen) {
      nf++;
      continue;
    }
    if(chf==1 && !glucose) {
      nf++;
      continue;
    }
    if(chf==2 && !hydrogenIon) {
      nf++;
      continue;
    }
    if (chf == 0) {
      strcpy(fieldName[nf], "oxygen");
      fieldDt[nf] = gfDt;
      fieldCriticalLevel1[nf] *= gfDt;	//*(boxVolume/cellVolume);
      fieldCriticalLevel2[nf] *= gfDt;	//*(boxVolume/cellVolume);
    }
    if (chf == 1) {
      strcpy(fieldName[nf], "glucose");
      fieldDt[nf] = secondsPerStep;
    }
    if (chf == 2) {
      strcpy(fieldName[nf], "hydrogenIon");
      fieldDt[nf] = secondsPerStep;
    }
    fieldType[nf] = SCALAR_FIELD;
    fieldAddr[nf] =
      (double *) calloc(gridSize.x * gridSize.y * gridSize.z,
                        sizeof(double));
    chemField[chf] = (double *) fieldAddr[nf];

    /* gradient allocation */
    gradAddr[chf] = (double*) calloc(gridSize.x * gridSize.y * gridSize.z * 3, sizeof(double));

    for (i = 0; i < gridSize.x * gridSize.y * gridSize.z; i++)
      chemField[chf][i] = 0.0;

    for (i=0; i<gridSize.x*gridSize.y*gridSize.z*3; i++)
      gradAddr[chf][i] = 0.0;

    nf++;
  }

  /* initialize temperature field */
  if (temperature) {
    tempEnvInitSystem();
    tempEnvInitBC();
    tempEnvInitSolver();
  }

  /* initialize chemical fields */
  for (i = 0; i < NCHEM; i++) {
    if (i == 0 && !oxygen)
      continue;
    if (i == 1 && !glucose)
      continue;
    if (i == 2 && !hydrogenIon)
      continue;
    chemEnvInitSystem(i);
    chemEnvInitBC(i);
    chemEnvInitSolver(i);
  }
}

/*!
 * This is a driving function for global field solver.
 * It is called in each step of the simulation.
 */
void fieldsSolve()
{
  int i;
  if (!gfields)
    return;
  for (gfIter = 0; gfIter < gfIterPerStep; gfIter++) {
    /* update cell state (if more than one iteration) */
    if (gfIter > 0)
      updateCellStates();
    /* solve temperature field */
    if (temperature)
      tempEnvSolve();
    /* solve chemical fields */
    for (i = 0; i < NCHEM; i++) {
      if (i == 0 && !oxygen)
        continue;
      if (i == 1 && !glucose)
        continue;
      if (i == 2 && !hydrogenIon)
        continue;
      chemEnvSolve(i);
    }
  }
}

void allocateFieldGradient()
{

  if (!gfields)
    return;

  MPI_Cart_shift(MPI_CART_COMM,0,1,&rX0,&rX1);
  MPI_Cart_shift(MPI_CART_COMM,1,1,&rY0,&rY1);
  MPI_Cart_shift(MPI_CART_COMM,2,1,&rZ0,&rZ1);

  /* allocate send buffers */
  if(rX0!=MPI_PROC_NULL) haloSX0=(double*)calloc(gridSize.y*gridSize.z,sizeof(double));
  if(rX1!=MPI_PROC_NULL) haloSX1=(double*)calloc(gridSize.y*gridSize.z,sizeof(double));
  if(rY0!=MPI_PROC_NULL) haloSY0=(double*)calloc(gridSize.x*gridSize.z,sizeof(double));
  if(rY1!=MPI_PROC_NULL) haloSY1=(double*)calloc(gridSize.x*gridSize.z,sizeof(double));
  if(rZ0!=MPI_PROC_NULL) haloSZ0=(double*)calloc(gridSize.y*gridSize.x,sizeof(double));
  if(rZ1!=MPI_PROC_NULL) haloSZ1=(double*)calloc(gridSize.y*gridSize.x,sizeof(double));

  /* allocate receive buffers */
  if(rX0!=MPI_PROC_NULL) haloRX0=(double*)calloc(gridSize.y*gridSize.z,sizeof(double));
  if(rX1!=MPI_PROC_NULL) haloRX1=(double*)calloc(gridSize.y*gridSize.z,sizeof(double));
  if(rY0!=MPI_PROC_NULL) haloRY0=(double*)calloc(gridSize.x*gridSize.z,sizeof(double));
  if(rY1!=MPI_PROC_NULL) haloRY1=(double*)calloc(gridSize.x*gridSize.z,sizeof(double));
  if(rZ0!=MPI_PROC_NULL) haloRZ0=(double*)calloc(gridSize.y*gridSize.x,sizeof(double));
  if(rZ1!=MPI_PROC_NULL) haloRZ1=(double*)calloc(gridSize.y*gridSize.x,sizeof(double));

  return;
}

void initFieldHaloExchange(int chf)
{
  int i,j,k;
  if (!gfields)
    return;

  for(i=0; i<gridSize.x; i++)
    for(j=0; j<gridSize.y; j++)
      for(k=0; k<gridSize.z; k++) {
        double val;
        val=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+k];
        if(i==0 && rX0!=MPI_PROC_NULL) haloSX0[gridSize.z*j+k]=val;
        if(i==gridSize.x-1 && rX1!=MPI_PROC_NULL) haloSX1[gridSize.z*j+k]=val;
        if(j==0 && rY0!=MPI_PROC_NULL) haloSY0[gridSize.z*i+k]=val;
        if(j==gridSize.y-1 && rY1!=MPI_PROC_NULL) haloSY1[gridSize.z*i+k]=val;
        if(k==0 && rZ0!=MPI_PROC_NULL) haloSZ0[gridSize.y*i+j]=val;
        if(k==gridSize.z-1 && rZ1!=MPI_PROC_NULL) haloSZ1[gridSize.y*i+j]=val;
      }

  if(rX0!=MPI_PROC_NULL) {
    MPI_Isend(haloSX0,gridSize.y*gridSize.z,MPI_DOUBLE,rX0,MPIrank,MPI_CART_COMM,&reqFGSend[0]);
    MPI_Irecv(haloRX0,gridSize.y*gridSize.z,MPI_DOUBLE,rX0,MPIsize+rX0,MPI_CART_COMM,&reqFGRecv[0]);
  }
  if(rX1!=MPI_PROC_NULL) {
    MPI_Isend(haloSX1,gridSize.y*gridSize.z,MPI_DOUBLE,rX1,MPIsize+MPIrank,MPI_CART_COMM,&reqFGSend[1]);
    MPI_Irecv(haloRX1,gridSize.y*gridSize.z,MPI_DOUBLE,rX1,rX1,MPI_CART_COMM,&reqFGRecv[1]);
  }
  if(rY0!=MPI_PROC_NULL) {
    MPI_Isend(haloSY0,gridSize.x*gridSize.z,MPI_DOUBLE,rY0,2*MPIsize+MPIrank,MPI_CART_COMM,&reqFGSend[2]);
    MPI_Irecv(haloRY0,gridSize.x*gridSize.z,MPI_DOUBLE,rY0,3*MPIsize+rY0,MPI_CART_COMM,&reqFGRecv[2]);
  }
  if(rY1!=MPI_PROC_NULL) {
    MPI_Isend(haloSY1,gridSize.x*gridSize.z,MPI_DOUBLE,rY1,3*MPIsize+MPIrank,MPI_CART_COMM,&reqFGSend[3]);
    MPI_Irecv(haloRY1,gridSize.x*gridSize.z,MPI_DOUBLE,rY1,2*MPIsize+rY1,MPI_CART_COMM,&reqFGRecv[3]);
  }
  if(rZ0!=MPI_PROC_NULL) {
    MPI_Isend(haloSZ0,gridSize.y*gridSize.x,MPI_DOUBLE,rZ0,4*MPIsize+MPIrank,MPI_CART_COMM,&reqFGSend[4]);
    MPI_Irecv(haloRZ0,gridSize.y*gridSize.x,MPI_DOUBLE,rZ0,5*MPIsize+rZ0,MPI_CART_COMM,&reqFGRecv[4]);
  }
  if(rZ1!=MPI_PROC_NULL) {
    MPI_Isend(haloSZ1,gridSize.y*gridSize.x,MPI_DOUBLE,rZ1,5*MPIsize+MPIrank,MPI_CART_COMM,&reqFGSend[5]);
    MPI_Irecv(haloRZ1,gridSize.y*gridSize.x,MPI_DOUBLE,rZ1,4*MPIsize+rZ1,MPI_CART_COMM,&reqFGRecv[5]);
  }

  return;
}

void waitFieldHaloExchange()
{
  if (!gfields)
    return;
  MPI_Status status;
  if(rX0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(rX1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(rY0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(rY1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(rZ0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(rZ1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
    if (MPI_Wait(&reqFGRecv[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }

  return;
}


void computeFieldGradient(int chf)
{
  int i,j,k;
  if (!gfields)
    return;
  /* compute internal data */
  for(i=0; i<gridSize.x; i++)
    for(j=0; j<gridSize.y; j++)
      for(k=0; k<gridSize.z; k++) {
        /* x coord gradient */
        if(i!=0 && i!=gridSize.x-1) {
          gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*j+3*k+0]=
            (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*(i+1)+gridSize.z*j+k]
             -fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*(i-1)+gridSize.z*j+k])/gridResolution;

        }
        /* y coord gradient */
        if(j!=0 && j!=gridSize.y-1) {
          gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*j+3*k+1]=
            (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*(j+1)+k]
             -fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*(j-1)+k])/gridResolution;
        }
        /* z coord gradient */
        if(k!=0 && k!=gridSize.z-1) {
          gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*j+3*k+2]=
            (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+k+1]
             -fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+k-1])/gridResolution;

        }
      }
  /* wait for boundary data to arrive */
  waitFieldHaloExchange();
  /* update with boundary data */
  for(j=0; j<gridSize.y; j++)
    for(k=0; k<gridSize.z; k++) {
      double x0,x1;
      if(rX0!=MPI_PROC_NULL) x0=haloRX0[gridSize.z*j+k];
      else x0=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*0+gridSize.z*j+k];
      if(rX1!=MPI_PROC_NULL) x1=haloRX1[gridSize.z*j+k];
      else x1=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*(gridSize.x-1)+gridSize.z*j+k];
      gradAddr[chf][3*gridSize.z*gridSize.y*0+3*gridSize.z*j+3*k+0]=
        (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*1+gridSize.z*j+k]
         - x0)/gridResolution;
      gradAddr[chf][3*gridSize.z*gridSize.y*(gridSize.x-1)+3*gridSize.z*j+3*k+0]=
        (x1
         - fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*(gridSize.x-2)+gridSize.z*j+k])/gridResolution;
    }
  for(i=0; i<gridSize.x; i++)
    for(k=0; k<gridSize.z; k++) {
      double y0,y1;
      if(rY0!=MPI_PROC_NULL) y0=haloRY0[gridSize.z*i+k];
      else y0=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*0+k];
      if(rY1!=MPI_PROC_NULL) y1=haloRY1[gridSize.z*i+k];
      else y1=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*(gridSize.y-1)+k];
      gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*0+3*k+1]=
        (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*1+k]
         - y0)/gridResolution;
      gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*(gridSize.y-1)+3*k+1]=
        (y1
         - fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*(gridSize.y-2)+k])/gridResolution;
    }
  for(i=0; i<gridSize.x; i++)
    for(j=0; j<gridSize.y; j++) {
      double z0,z1;
      if(rZ0!=MPI_PROC_NULL) z0=haloRZ0[gridSize.y*i+j];
      else z0=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+0];
      if(rZ1!=MPI_PROC_NULL) z1=haloRZ1[gridSize.y*i+j];
      else z1=fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+(gridSize.z-1)];
      gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*j+3*0+2]=
        (fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+1]
         - z0)/gridResolution;
      gradAddr[chf][3*gridSize.z*gridSize.y*i+3*gridSize.z*j+3*(gridSize.z-1)+2]=
        (z1
         - fieldAddr[NGLOB+chf][gridSize.z*gridSize.y*i+gridSize.z*j+(gridSize.z-2)])/gridResolution;
    }

  if(rX0!=MPI_PROC_NULL) {
    free(haloSX0);
    free(haloRX0);
  }
  if(rX1!=MPI_PROC_NULL) {
    free(haloSX1);
    free(haloRX1);
  }
  if(rY0!=MPI_PROC_NULL) {
    free(haloSY0);
    free(haloRY0);
  }
  if(rY1!=MPI_PROC_NULL) {
    free(haloSY1);
    free(haloRY1);
  }
  if(rZ0!=MPI_PROC_NULL) {
    free(haloSZ0);
    free(haloRZ0);
  }
  if(rZ1!=MPI_PROC_NULL) {
    free(haloSZ1);
    free(haloRZ1);
  }

  return;
}

void fieldGradient()
{
  int chf;
  if (!gfields)
    return;
  allocateFieldGradient();
  for(chf=0; chf<NCHEM; chf++) {
    if(chf==OXYG-NGLOB && !oxygen) continue;
    if(chf==GLUC-NGLOB && !glucose) continue;
    if(chf==HYDR-NGLOB && !hydrogenIon) continue;
    initFieldHaloExchange(chf);
    //waitFieldHaloExchange();
    computeFieldGradient(chf);

  }
  return;
}


