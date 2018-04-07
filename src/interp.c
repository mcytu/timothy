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
#include <limits.h>

#include "global.h"
#include "fields.h"
#include "utils.h"

/*! \file interp.c
 *  \brief contains grid to cellular data interpolation functions
 */

double **cicPatch;
int *cicIntersect;
int *cicReceiver;
int *cicSender;
double **cicRecvPatch;

MPI_Request *cicReqSend, *cicReqRecv;

int *recvP;

int643dv_t *lowerPatchCorner, *upperPatchCorner;
int643dv_t *lowerPatchCornerR, *upperPatchCornerR;
int643dv_t *patchSize;
//int643dv_t *patchSizeR;
#define patch(p,f,i,j,k) (cicPatch[p][patchSize[p].x*patchSize[p].y*patchSize[p].z*f+patchSize[p].y*patchSize[p].z*i+patchSize[p].z*j+k])

/*
   void createinterpdata(systeminfo_t systeminfo,settings_t settings,grid_t grid,cellsinfo_t cellsinfo) {
        int i,j;
        int nexport[systeminfo.size];
        for(i=0; i<systeminfo.size; i++)
                nexport[i]=0;

        for(i=0; i<cellsinfo.localcount.n; i++) {
                double x,y,z;
                x=cellsinfo.cells[i].x;
                y=cellsinfo.cells[i].y;
                z=cellsinfo.cells[i].z;
                for(j=0; j<systeminfo.size; j++) {
                        if(x>=grid.lowleftnear[j].x && x<grid.uprightfar[j].x &&
                           y>=grid.lowleftnear[j].y && y<grid.uprightfar[j].y &&
                           z>=grid.lowleftnear[j].z && z<grid.uprightfar[j].z ) {
                                nexport[j]++;
                                break;
                        }
                }
        }
        printf("%d n=%d \n",systeminfo.rank,cellsinfo.localcount.n);
        for(i=0; i<systeminfo.size; i++)
                printf("%d %d:%d \n",systeminfo.rank,i,nexport[i]);


   }/*



   /*!
 * For each local cell we check its position in the grid.
 * If cell is located in the grid partition of other process
 * then the information about this cell should be send to
 * remote process by MPI communication. This is done by preparing
 * special patches and sending the information on all bunch of cells
 * rather than sendind it one by one.
 * This function identifies patches and allocate memory buffers for patches.
 */
void findpatches(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, grid_t *grid)
{

        int c,p;
        int643dv_t cellIdx;

        cicPatch = (double **) calloc(systeminfo.size, sizeof(double *));
        cicIntersect = (int *) calloc(systeminfo.size, sizeof(int));
        lowerPatchCorner = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        upperPatchCorner = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        lowerPatchCornerR = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        upperPatchCornerR = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        patchSize = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));

        for (p = 0; p < systeminfo.size; p++) {
                cicIntersect[p] = 0;
                patchSize[p].x = 0;
                patchSize[p].y = 0;
                patchSize[p].z = 0;
                lowerPatchCorner[p].x = INT_MAX;
                lowerPatchCorner[p].y = INT_MAX;
                if (cellsinfo->dimension == 3)
                        lowerPatchCorner[p].z = INT_MAX;
                else
                        lowerPatchCorner[p].z = 0;
                upperPatchCorner[p].x = INT_MIN;
                upperPatchCorner[p].y = INT_MIN;
                if (cellsinfo->dimension == 3)
                        upperPatchCorner[p].z = INT_MIN;
                else
                        upperPatchCorner[p].z = 0;
        }

        //#pragma omp parallel for default(none) private(p,c,cellIdx) shared(cells,gridResolution,lnc,MPIsize,gridStartIdx,gridEndIdx,lowerPatchCorner,upperPatchCorner,cicIntersect,sdim,lowerGridCorner)
        for (p = 0; p < systeminfo.size; p++) {

                for (c = 0; c < cellsinfo->localcount.n; c++) {

                        int ax, ay, az;

                        cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                        cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                        cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);

                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {

                                                if (cellIdx.x + ax >= grid->loweridx[p].x
                                                    && cellIdx.y + ay >= grid->loweridx[p].y
                                                    && cellIdx.z + az >= grid->loweridx[p].z
                                                    && cellIdx.x + ax <= grid->upperidx[p].x
                                                    && cellIdx.y + ay <= grid->upperidx[p].y
                                                    && cellIdx.z + az <= grid->upperidx[p].z) {
                                                        lowerPatchCorner[p].x =
                                                                (lowerPatchCorner[p].x >
                                                                 cellIdx.x + ax ? cellIdx.x +
                                                                 ax : lowerPatchCorner[p].x);
                                                        lowerPatchCorner[p].y =
                                                                (lowerPatchCorner[p].y >
                                                                 cellIdx.y + ay ? cellIdx.y +
                                                                 ay : lowerPatchCorner[p].y);
                                                        if (cellsinfo->dimension == 3)
                                                                lowerPatchCorner[p].z =
                                                                        (lowerPatchCorner[p].z >
                                                                         cellIdx.z + az ? cellIdx.z +
                                                                         az : lowerPatchCorner[p].z);
                                                        upperPatchCorner[p].x =
                                                                (upperPatchCorner[p].x <
                                                                 cellIdx.x + ax ? cellIdx.x +
                                                                 ax : upperPatchCorner[p].x);
                                                        upperPatchCorner[p].y =
                                                                (upperPatchCorner[p].y <
                                                                 cellIdx.y + ay ? cellIdx.y +
                                                                 ay : upperPatchCorner[p].y);
                                                        if (cellsinfo->dimension == 3)
                                                                upperPatchCorner[p].z =
                                                                        (upperPatchCorner[p].z <
                                                                         cellIdx.z + az ? cellIdx.z +
                                                                         az : upperPatchCorner[p].z);

                                                        cicIntersect[p] = 1;
                                                }
                                        }
                }
        }

        for (p = 0; p < systeminfo.size; p++)
                if (cicIntersect[p]) {
                        patchSize[p].x = upperPatchCorner[p].x - lowerPatchCorner[p].x + 1;
                        patchSize[p].y = upperPatchCorner[p].y - lowerPatchCorner[p].y + 1;
                        if (cellsinfo->dimension == 3)
                                patchSize[p].z = upperPatchCorner[p].z - lowerPatchCorner[p].z + 1;
                        else
                                patchSize[p].z = 1;
                        cicPatch[p] =
                                (double *) calloc(patchSize[p].x * patchSize[p].y *
                                                  patchSize[p].z * settings.numberoffields, sizeof(double));
                }

        return;
}

/*!
 * For each local cell its density value is interpolated accross
 * neighbouring grid vertices with the use of Cloud-In-Cell method.
 * Computed values are stored in patches instead in field buffers.
 * No additional memory allocations are made here.
 */
void cells2envinfo(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, grid_t *grid)
{

        int c, i, p, j, k, f;
        int643dv_t idx;
        double3dv_t d, t;
        int643dv_t cellIdx;
        double3dv_t cicCoord;

        for (p = 0; p < systeminfo.size; p++)
                for(f = 0; f < settings.numberoffields; f++)
                        for (i = 0; i < patchSize[p].x; i++)
                                for (j = 0; j < patchSize[p].y; j++)
                                        for (k = 0; k < patchSize[p].z; k++)
                                                patch(p, f, i, j, k) = 0.0;

        for (c = 0; c < cellsinfo->localcount.n; c++) {

                cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);

                for (p = 0; p < systeminfo.size; p++) {
                        int ax, ay, az;
                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {
                                                if (cellIdx.x + ax >= grid->loweridx[p].x
                                                    && cellIdx.y + ay >= grid->loweridx[p].y
                                                    && cellIdx.z + az >= grid->loweridx[p].z
                                                    && cellIdx.x + ax <= grid->upperidx[p].x
                                                    && cellIdx.y + ay <= grid->upperidx[p].y
                                                    && cellIdx.z + az <= grid->upperidx[p].z) {

                                                        idx.x = (cellIdx.x + ax) - lowerPatchCorner[p].x;
                                                        idx.y = (cellIdx.y + ay) - lowerPatchCorner[p].y;
                                                        idx.z = (cellIdx.z + az) - lowerPatchCorner[p].z;

                                                        cicCoord.x = grid->lowercorner.x + cellIdx.x * grid->resolution;
                                                        cicCoord.y = grid->lowercorner.y + cellIdx.y * grid->resolution;
                                                        cicCoord.z = grid->lowercorner.z + cellIdx.z * grid->resolution;

                                                        d.x = (cellsinfo->cells[c].x - cicCoord.x) / grid->resolution;
                                                        d.y = (cellsinfo->cells[c].y - cicCoord.y) / grid->resolution;
                                                        d.z = (cellsinfo->cells[c].z - cicCoord.z) / grid->resolution;

                                                        t.x = 1.0 - d.x;
                                                        t.y = 1.0 - d.y;
                                                        t.z = 1.0 - d.z;

                                                        /*if(cells[c].ctype==1) { endothelial cell - production
                                                                patch(p,1,idx.x,idx.y,idx.z) +=
                                                                        1.0 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                             (1 -
                                                                                                              ay) * t.y) *
                                                                        (az * d.z + (1 - az) * t.z);
                                                           } else */
                                                        if (cellsinfo->cells[c].phase != 5) { /* if not in necrotic phase */
                                                                if (cellsinfo->cells[c].phase == 0) { /* if in G0 phase - lower consumption */
                                                                        patch(p,0, idx.x, idx.y, idx.z) +=
                                                                                0.75 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                                      (1 -
                                                                                                                       ay) * t.y) *
                                                                                (az * d.z + (1 - az) * t.z);
                                                                } else { /* if not in G0 phase - normal consumption */
                                                                        patch(p,0, idx.x, idx.y, idx.z) +=
                                                                                1.0 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                                     (1 -
                                                                                                                      ay) * t.y) *
                                                                                (az * d.z + (1 - az) * t.z);
                                                                }
                                                        }

                                                }
                                        }
                }
        }
        return;
}

/*!
 *  Patches are being sent to receiving processes with non-blocking
 * communication scheme (Isend, Irecv).
 * Receiving patches are allocated here.
 * MPI_Request tables are allocated here.
 */
void commcells2envinit(systeminfo_t systeminfo, settings_t settings)
{
        int p;

        cicReceiver = (int *) calloc(systeminfo.size, sizeof(int));
        cicSender = (int *) calloc(systeminfo.size, sizeof(int));

        for (p = 0; p < systeminfo.size; p++)
                cicReceiver[p] = cicIntersect[p];

        MPI_Alltoall(cicReceiver, 1, MPI_INT, cicSender, 1, MPI_INT,
                     MPI_COMM_WORLD);

        MPI_Alltoall(lowerPatchCorner, sizeof(int643dv_t), MPI_BYTE,
                     lowerPatchCornerR, sizeof(int643dv_t), MPI_BYTE,
                     MPI_COMM_WORLD);
        MPI_Alltoall(upperPatchCorner, sizeof(int643dv_t), MPI_BYTE,
                     upperPatchCornerR, sizeof(int643dv_t), MPI_BYTE,
                     MPI_COMM_WORLD);

        cicRecvPatch = (double **) calloc(systeminfo.size, sizeof(double *));
        cicReqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        cicReqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        for (p = 0; p < systeminfo.size; p++) {
                if (cicReceiver[p]) {
                        MPI_Isend(&(cicPatch[p][0]),
                                  patchSize[p].x * patchSize[p].y * patchSize[p].z * settings.numberoffields,
                                  MPI_DOUBLE, p, MPIrank, MPI_COMM_WORLD, &cicReqSend[p]);
                }
                if (cicSender[p]) {
                        int recvSize;
                        recvSize =
                                (upperPatchCornerR[p].x - lowerPatchCornerR[p].x +
                                 1) * (upperPatchCornerR[p].y - lowerPatchCornerR[p].y +
                                       1) * (upperPatchCornerR[p].z - lowerPatchCornerR[p].z +
                                             1);
                        if (!(cicRecvPatch[p] = (double *) calloc(recvSize*settings.numberoffields, sizeof(double))))
                                exit(0);
                        MPI_Irecv(&(cicRecvPatch[p][0]), recvSize * settings.numberoffields, MPI_DOUBLE, p, p,
                                  MPI_COMM_WORLD, &cicReqRecv[p]);
                }
        }
        return;
}

/*!
 * Wait for communication to finish.
 * MPI_Request tables are deallocated here.
 */
int commcells2envwait(systeminfo_t systeminfo)
{
        int p;
        MPI_Status status;

        for (p = 0; p < systeminfo.size; p++) {
                if (!cicReceiver[p])
                        continue;
                if (MPI_Wait(&cicReqSend[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commcells2envwait 1", __FILE__, __LINE__);
        }

        for (p = 0; p < systeminfo.size; p++) {
                if (!cicSender[p])
                        continue;
                if (MPI_Wait(&cicReqRecv[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commcells2envwait 2", __FILE__, __LINE__);
        }

        free(cicReqSend);
        free(cicReqRecv);

        return 0;
}

/*!
 * Update local tissueField part with information from remote processes
 * received in patches. Receiveing patches are deallocated here.
 */
int applycells2envpatches(systeminfo_t systeminfo, settings_t settings, grid_t *grid,environment_t **environment)
{

        int i,p,f;
        int index;

        for(i=0; i<settings.numberoffields; i++) {
                memset((*environment)[i].production,0,grid->localsize.x*grid->localsize.y*grid->localsize.z*sizeof(double));
        }

        for (p = 0; p < systeminfo.size; p++) {
                int i, j, k;
                if (!cicSender[p])
                        continue;
                for (i = lowerPatchCornerR[p].x; i <= upperPatchCornerR[p].x; i++)
                        for (j = lowerPatchCornerR[p].y; j <= upperPatchCornerR[p].y; j++)
                                for (k = lowerPatchCornerR[p].z; k <= upperPatchCornerR[p].z; k++) {
                                        int643dv_t c, g, size;
                                        size.x = upperPatchCornerR[p].x - lowerPatchCornerR[p].x + 1;
                                        size.y = upperPatchCornerR[p].y - lowerPatchCornerR[p].y + 1;
                                        size.z = upperPatchCornerR[p].z - lowerPatchCornerR[p].z + 1;
                                        c.x = i - lowerPatchCornerR[p].x;
                                        c.y = j - lowerPatchCornerR[p].y;
                                        c.z = k - lowerPatchCornerR[p].z;
                                        if (i >= grid->loweridx[systeminfo.rank].x && i <= grid->upperidx[systeminfo.rank].x
                                            && j >= grid->loweridx[systeminfo.rank].y && j <= grid->upperidx[systeminfo.rank].y
                                            && k >= grid->loweridx[systeminfo.rank].z && k <= grid->upperidx[systeminfo.rank].z) {
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index=grid->localsize.z * grid->localsize.y * g.x + grid->localsize.z * g.y + g.z;
                                                for(f=0; f<settings.numberoffields; f++)
                                                        (*environment)[f].production[index] +=
                                                                cicRecvPatch[p][f*(size.x*size.y*size.z)+size.z * size.y * c.x + size.z * c.y + c.z];
                                        }
                                }
                free(cicRecvPatch[p]);
        }

        free(cicRecvPatch);

        for (p = 0; p < systeminfo.size; p++)
                if (cicIntersect[p])
                        free(cicPatch[p]);

        free(cicPatch);
        return 0;
}

/*!
 * Now each local cell should receive information about values
 * of global fields in the grid.
 * Field patches are filled with appropriate values.
 * Non-blocking communication is initiated.
 * This is executed after global fields are computed.
 * Field patches buffers are allocated here.
 * Sizes of the patches are the same as those from previous CIC communication.
 * Receiving field patches are also allocated here.
 * MPI_Request tables are allocated here.
 */
void env2cellsinfo(systeminfo_t systeminfo, settings_t settings, grid_t *grid,environment_t **environment)
{

        int f; /* fields index */
        int p; /* process index */
        int643dv_t idx, g;
        int643dv_t size;

        fieldsPatchesCommBuff = (double **) calloc(systeminfo.size, sizeof(double *));
        for (p = 0; p < systeminfo.size; p++) {
                int fieldPatchSize;
                if (!cicSender[p])
                        continue; /* continue to next process if current process do not overlap domain */
                size.x = upperPatchCornerR[p].x - lowerPatchCornerR[p].x + 1;
                size.y = upperPatchCornerR[p].y - lowerPatchCornerR[p].y + 1;
                size.z = upperPatchCornerR[p].z - lowerPatchCornerR[p].z + 1;
                fieldPatchSize = size.x * size.y * size.z;
                fieldsPatchesCommBuff[p] =
                        (double *) calloc((settings.numberoffields*4) * fieldPatchSize, sizeof(double));
                /* fields */
                for (f = 0; f < settings.numberoffields; f++) {
                        int64_t i, j, k;
                        for (i = lowerPatchCornerR[p].x; i <= upperPatchCornerR[p].x; i++)
                                for (j = lowerPatchCornerR[p].y; j <= upperPatchCornerR[p].y; j++)
                                        for (k = lowerPatchCornerR[p].z; k <= upperPatchCornerR[p].z; k++) {
                                                int index1,index2;
                                                idx.x = i - lowerPatchCornerR[p].x;
                                                idx.y = j - lowerPatchCornerR[p].y;
                                                idx.z = k - lowerPatchCornerR[p].z;
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index1 = f * fieldPatchSize + size.z * size.y * idx.x + size.z * idx.y + idx.z;
                                                index2 = grid->localsize.z * grid->localsize.y * g.x + grid->localsize.z * g.y + g.z;
                                                fieldsPatchesCommBuff[p][index1] = (*environment)[f].data[index2];
                                        }
                }
                /* field gradients */
                for (f = 0; f < settings.numberoffields*3; f++) {
                        int64_t i, j, k;
                        for (i = lowerPatchCornerR[p].x; i <= upperPatchCornerR[p].x; i++)
                                for (j = lowerPatchCornerR[p].y; j <= upperPatchCornerR[p].y; j++)
                                        for (k = lowerPatchCornerR[p].z; k <= upperPatchCornerR[p].z; k++) {
                                                int index1,index2;
                                                idx.x = i - lowerPatchCornerR[p].x;
                                                idx.y = j - lowerPatchCornerR[p].y;
                                                idx.z = k - lowerPatchCornerR[p].z;
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index1 = fieldPatchSize*settings.numberoffields + 3* (f * fieldPatchSize + size.z * size.y * idx.x + size.z * idx.y + idx.z);
                                                index2 = 3*(grid->localsize.z*grid->localsize.y*g.x + grid->localsize.z*g.y + g.z);
                                                fieldsPatchesCommBuff[p][index1] = (*environment)[f].gradient[index2];
                                                fieldsPatchesCommBuff[p][index1+1] = (*environment)[f].gradient[index2+1];
                                                fieldsPatchesCommBuff[p][index1+2] = (*environment)[f].gradient[index2+2];
                                        }
                }

        }
        return;
}

void commenv2cellsinit(systeminfo_t systeminfo, settings_t settings){

        int p; /* process index */

        if (!(fieldsPatches = (double **) calloc(systeminfo.size, sizeof(double *))))
                terminate(systeminfo,"cannot allocate fieldsPatches", __FILE__, __LINE__);

        cicReqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        cicReqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        for (p = 0; p < systeminfo.size; p++) {
                if (cicSender[p]) {
                        int sendSize;
                        sendSize =  (upperPatchCornerR[p].x - lowerPatchCornerR[p].x +
                                     1) * (upperPatchCornerR[p].y - lowerPatchCornerR[p].y +
                                           1) * (upperPatchCornerR[p].z - lowerPatchCornerR[p].z +
                                                 1) * (settings.numberoffields*4);
                        MPI_Isend(&(fieldsPatchesCommBuff[p][0]), sendSize, MPI_DOUBLE, p,
                                  MPIrank, MPI_COMM_WORLD, &cicReqSend[p]);
                }
                if (cicReceiver[p]) {
                        int recvSize;
                        recvSize = patchSize[p].x * patchSize[p].y * patchSize[p].z * (settings.numberoffields*4);
                        if (!(fieldsPatches[p] = (double *) calloc(recvSize, sizeof(double))))
                                terminate(systeminfo,"cannot allocate fieldsPatches[p]", __FILE__, __LINE__);
                        MPI_Irecv(&(fieldsPatches[p][0]), recvSize, MPI_DOUBLE, p, p,
                                  MPI_COMM_WORLD, &cicReqRecv[p]);
                }
        }
        return;
}

/*!
 * Wait for communication to finish.
 * MPI_Request tables are deallocated here.
 * Field patches buffers are deallocated here.
 */
void commenv2cellswait(systeminfo_t systeminfo)
{
        int p;
        MPI_Status status;

        for (p = 0; p < systeminfo.size; p++) {
                if (!cicSender[p])
                        continue;
                if (MPI_Wait(&cicReqSend[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commenv2cellswait 1", __FILE__, __LINE__);
        }

        for (p = 0; p < systeminfo.size; p++) {
                if (!cicReceiver[p])
                        continue;
                if (MPI_Wait(&cicReqRecv[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commenv2cellswait 2", __FILE__, __LINE__);
        }

        free(cicReqSend);
        free(cicReqRecv);

        for (p = 0; p < systeminfo.size; p++)
                free(fieldsPatchesCommBuff[p]);
        free(fieldsPatchesCommBuff);

        return;
}

/*!
 * Update local cell fields with information from
 * remote processes received in patches.
 * Receiveing field patches are deallocated here.
 */
void applyenv2cellspatches(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo,grid_t *grid,cellenvdata_t ***cellenvdata)
{

        int p,c,f;
        int643dv_t idx;
        double3dv_t d, t;
        int643dv_t cellIdx;
        double3dv_t cicCoord;

        /* allocate cellenvdata */
        if(!((*cellenvdata)=(cellenvdata_t**)calloc(settings.numberoffields,sizeof(cellenvdata_t*))))
                terminate(systeminfo,"cannot allocate cellenvdata", __FILE__, __LINE__);
        for(f=0; f<settings.numberoffields; f++)
                if(!((*cellenvdata)[f]=(cellenvdata_t*)calloc(cellsinfo->localcount.n,sizeof(cellenvdata_t))))
                        terminate(systeminfo,"cannot allocate cellenvdata[f]", __FILE__, __LINE__);

        /* reset fields */
        for (f = 0; f < settings.numberoffields; f++)
                for (c = 0; c < cellsinfo->localcount.n; c++) {
                        (*cellenvdata)[f][c].value=0.0;
                        (*cellenvdata)[f][c].gx=0.0;
                        (*cellenvdata)[f][c].gy=0.0;
                        (*cellenvdata)[f][c].gz=0.0;
                }

        for (c = 0; c < cellsinfo->localcount.n; c++) { /* for every cell */
                cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);
                for (p = 0; p < systeminfo.size; p++) { /* for each process */
                        int ax, ay, az;
                        if (!cicReceiver[p])
                                continue; /* there is no patch from this process */
                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {
                                                if (cellIdx.x + ax >= gridStartIdx[p].x
                                                    && cellIdx.y + ay >= gridStartIdx[p].y
                                                    && cellIdx.z + az >= gridStartIdx[p].z
                                                    && cellIdx.x + ax <= gridEndIdx[p].x
                                                    && cellIdx.y + ay <= gridEndIdx[p].y
                                                    && cellIdx.z + az <= gridEndIdx[p].z) {

                                                        idx.x = (cellIdx.x + ax) - lowerPatchCorner[p].x;
                                                        idx.y = (cellIdx.y + ay) - lowerPatchCorner[p].y;
                                                        idx.z = (cellIdx.z + az) - lowerPatchCorner[p].z;

                                                        cicCoord.x = grid->lowercorner.x + cellIdx.x * grid->resolution;
                                                        cicCoord.y = grid->lowercorner.y + cellIdx.y * grid->resolution;
                                                        cicCoord.z = grid->lowercorner.z + cellIdx.z * grid->resolution;

                                                        d.x = (cellsinfo->cells[c].x - cicCoord.x) / grid->resolution;
                                                        d.y = (cellsinfo->cells[c].y - cicCoord.y) / grid->resolution;
                                                        d.z = (cellsinfo->cells[c].z - cicCoord.z) / grid->resolution;

                                                        t.x = 1.0 - d.x;
                                                        t.y = 1.0 - d.y;
                                                        t.z = 1.0 - d.z;

                                                        /* interpolating back to cells */
                                                        /* scaling from mol/cm^3 to mol/cell */
                                                        for (f = 0; f < settings.numberoffields; f++) {
                                                                int index1 = f * patchSize[p].x * patchSize[p].y * patchSize[p].z + patchSize[p].y * patchSize[p].z * idx.x + patchSize[p].z * idx.y + idx.z;
                                                                (*cellenvdata)[f][c].value += fieldsPatches[p][index1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z); //*cellVolume;
                                                        }
                                                        for (f = 0; f < settings.numberoffields; f++) {
                                                                int index1=settings.numberoffields * patchSize[p].x * patchSize[p].y * patchSize[p].z + 3* (f * patchSize[p].x * patchSize[p].y * patchSize[p].z + patchSize[p].y * patchSize[p].z * idx.x + patchSize[p].z * idx.y + idx.z);
                                                                (*cellenvdata)[f][c].gx += fieldsPatches[p][index1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                                (*cellenvdata)[f][c].gy += fieldsPatches[p][index1 + 1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                                (*cellenvdata)[f][c].gz += fieldsPatches[p][index1 + 2] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                        }
                                                } // if
                                        }// az
                } // p
        } // c

        for (p = 0; p < systeminfo.size; p++)
                free(fieldsPatches[p]);
        free(fieldsPatches);

        return;
}

/*!
 * This is a driving function interpolating
 * cellular data to grid data.
 * This function does not enable overlapping communication
 * and computations.
 */
/*void interpolateCellsToGrid()
   {
        findPatches();
        doInterpolation();
        initPatchExchange();
        waitPatchExchange();
        applyPatches();
   }*/

/*!
 * This function initializes data exchange between
 * processes required in cells-to-grid interpolation.
 * This function enables overlapping communication and
 * computations.
 */
void cells2envinit(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, grid_t *grid)
{
        if (settings.numberoffields==0)
                return;
        findpatches(systeminfo,settings,cellsinfo,grid);
        cells2envinfo(systeminfo,settings,cellsinfo,grid);
        commcells2envinit(systeminfo,settings);
        return;
}

/*!
 * This function wait for patches communication to finish
 * in cells-to-grid interpolation.
 * This function enables overlapping communication and
 * computations.
 */
void cells2envwait(systeminfo_t systeminfo, settings_t settings, grid_t *grid, environment_t **environment)
{
        if (settings.numberoffields==0)
                return;
        commcells2envwait(systeminfo);
        applycells2envpatches(systeminfo,settings,grid,environment);
        return;
}

/*!
 * This function is used to interpolate field data back to
 * cells. Overlapping of communication and computations is not
 * implemented here.
 * This function deallocates all important arrays used in interpolation.
 */
void env2cells(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, grid_t *grid, environment_t **environment,cellenvdata_t ***cellenvdata)
{
        if (settings.numberoffields==0)
                return;

        env2cellsinfo(systeminfo,settings,grid,environment);
        commenv2cellsinit(systeminfo,settings);
        commenv2cellswait(systeminfo);
        applyenv2cellspatches(systeminfo,settings,cellsinfo,grid,cellenvdata);

        free(cicReceiver);
        free(cicSender);
        free(cicIntersect);

        free(lowerPatchCorner);
        free(upperPatchCorner);
        free(lowerPatchCornerR);
        free(upperPatchCornerR);
        free(patchSize);
        return;
}
