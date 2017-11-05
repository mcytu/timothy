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

#include <stdlib.h>

#include "global.h"

/*! \file comm.c
 *  \brief contains communication functions
 */

/* data type for IDs of exported cells */
struct expData {
        int cell;
        int proc;
};

MPI_Request *reqSend;
MPI_Request *reqRecv;

int64_t *sendOffset;
int64_t *recvOffset;

struct expData *expList;

int *recvCount;
int *sendCount;

#define MAX_EXPORTED_PER_PROC 2*maxCellsPerProc

/*!
 * This function is a comparison function used
 * for sorting export list table.
 */
int explistcompare(const void *a, const void *b)
{
        return ((explist_t *) a)->proc - ((explist_t *) b)->proc;
}

/*!
 * This function uses Zoltan's library function Zoltan_LB_Box_Assign
 * to find possible intersections of cells' neighbourhoods
 * and other processes' geometries.
 */
void createexportlist(system_t system,settings_t settings,cellsinfo_t cellsinfo,commdata_t *commdata)
{

        int i, p;
        int procs[system.size];
        int numprocs;

        commdata->numexp = 0;
        commdata->numimp = 0;

        commdata->explistmaxsize=settings.maxlocalcells;

        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;

        if(!(commdata->explist = (explist_t*) malloc(sizeof(explist_t) * commdata->explistmaxsize)))
                terminate(system,"cannot allocate commdata->explist", __FILE__, __LINE__);
        if(!(commdata->recvcount = (int *) calloc(system.size, sizeof(int))))
                terminate(system,"cannot allocate commdata->recvcount", __FILE__, __LINE__);
        if(!(commdata->sendcount = (int *) calloc(system.size, sizeof(int))))
                terminate(system,"cannot allocate commdata->sendcount", __FILE__, __LINE__);
        if(!(commdata->sendoffset = (int64_t *) calloc(system.size, sizeof(int64_t))))
                terminate(system,"cannot allocate commdata->sendoffset", __FILE__, __LINE__);
        if(!(commdata->recvoffset = (int64_t *) calloc(system.size, sizeof(int64_t))))
                terminate(system,"cannot allocate commdata->recvoffset", __FILE__, __LINE__);

        /* loop over local cells */
        /*#pragma omp parallel for private(procs) */
        for (p = 0; p < cellsinfo.localcount.n; p++) {
                double xmin, xmax, ymin, ymax, zmin, zmax;
                double r;

                if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC)
                        continue;

                r = cellsinfo.cells[p].h * 1.5;

                /* compute neighbourhood box */
                xmin = cellsinfo.cells[p].x - r;
                xmax = cellsinfo.cells[p].x + r;
                ymin = cellsinfo.cells[p].y - r;
                ymax = cellsinfo.cells[p].y + r;
                if (cellsinfo.dimension == 3) {
                        zmin = cellsinfo.cells[p].z - r;
                        zmax = cellsinfo.cells[p].z + r;
                } else {
                        zmin = 0.0;
                        zmax = 0.0;
                }

                /* look for possible neighbours */
                Zoltan_LB_Box_Assign(ztn, xmin, ymin, zmin, xmax, ymax, zmax, procs, &numprocs);

                /* loop over receivers */
                for (i = 0; i < numprocs; i++) {
                        if (procs[i] == system.rank || cellsinfo.cellsperproc[procs[i]] == 0)
                                continue;
                        commdata->explist[numexp].cell = p;
                        commdata->explist[numexp].proc = procs[i];
                        commdata->sendcount[procs[i]]++;
                        commdata->numexp++;
                        /* too many refugees  - reallocate */
                        if (commdata->numexp >= commdata->explistmaxsize) {
                                commdata->explistmaxsize+=64;
                                if(!(commdata->explist=(explist_t*)realloc(commdata->explist,sizeof(explist_t)*commdata->explistmaxsize))) {
                                        terminate(system,"cannot reallocate commdata->explist", __FILE__, __LINE__);
                                }
                        }
                }
        }

        /* sort export list with respect to process number */
        qsort(comdata->explist, commdata->numexp, sizeof(explist_t), explistcompare);

        /* distribute the information on transfer sizes between each process */
        MPI_Alltoall(commdata->sendcount, 1, MPI_INT, commdata->recvcount, 1, MPI_INT, MPI_COMM_WORLD);

        /* compute send offsets */
        commdata->sendoffset[0] = 0;
        for (i = 1; i < system.size; i++)
                commdata->sendoffset[i] = commdata->sendoffset[i - 1] + commdata->sendcount[i - 1];

        /* compute receive offsets */
        commdata->recvoffset[0] = 0;
        for (i = 1; i < system.size; i++)
                commdata->recvoffset[i] = commdata->recvoffset[i - 1] + commdata->recvcount[i - 1];

        /* count cells to be imported */
        for (i = 0; i < system.size; i++)
                commdata->numimp += commdata->recvcount[i];

        return;
}

/*!
 * This function deallocates all communication buffers and
 * auxiliary tables.
 */
void commcleanup(system_t system,cellsinfo_t cellsinfo,commdata_t *commdata)
{
        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;
        free(commdata->recvdata);
        free(commdata->recvdpdata);
        free(commdata->explist);
        free(commdata->recvcount);
        free(commdata->sendcount);
        free(commdata->sendoffset);
        free(commdata->recvoffset);
        return;
}

/*!
 * This function initiate sending and receiving cells' data between processes.
 */
void cellsExchangeInit()
{
        int i;

        if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
                return;

        /* allocate communication buffers */
        sendData = (struct partData *) malloc(numExp * sizeof(struct partData));
#ifdef __MIC__
        recvData = (struct partData *) _mm_malloc(numImp * sizeof(struct partData),64);
#else
        recvData = (struct partData *) malloc(numImp * sizeof(struct partData));
#endif
        reqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);
        reqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);

        /* create reduced particle data buffer for exporting */
        for (i = 0; i < numExp; i++) {
                sendData[i].x = cells[expList[i].cell].x;
                sendData[i].y = cells[expList[i].cell].y;
                sendData[i].z = cells[expList[i].cell].z;
                sendData[i].size = cells[expList[i].cell].size;
                sendData[i].h = h;
                sendData[i].young = (double) cells[expList[i].cell].young;
                sendData[i].ctype = cells[expList[i].cell].ctype;
        }

        /* send cells - asynchronous MPI call */
        for (i = 0; i < MPIsize; i++) {
                if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                MPI_Isend(&sendData[sendOffset[i]],
                          sendCount[i] * sizeof(struct partData), MPI_BYTE, i, MPIrank,
                          MPI_COMM_WORLD, &reqSend[i]);
        }

        /* receive cells - asynchronous MPI call */
        for (i = 0; i < MPIsize; i++) {
                if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                MPI_Irecv(&recvData[recvOffset[i]],
                          recvCount[i] * sizeof(struct partData), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &reqRecv[i]);
        }

}

/*!
 * This function waits for cells' data exchange completion.
 */
void cellsExchangeWait()
{
        int i;
        MPI_Status status;

        if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
                return;

        /* wait for send completion */
        for (i = 0; i < MPIsize; i++) {
                if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                if (MPI_Wait(&reqSend[i], &status) != MPI_SUCCESS)
                        stopRun(103, "reqSend", __FILE__, __LINE__);
        }

        /* wait for receive completion */
        for (i = 0; i < MPIsize; i++) {
                if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                if (MPI_Wait(&reqRecv[i], &status) != MPI_SUCCESS)
                        stopRun(103, "reqRecv", __FILE__, __LINE__);
        }

        /* some of the buffers can be deallocated here */
        free(sendData);
        free(reqSend);
        free(reqRecv);
}

/*!
 * This function initiate sending and receiving density
 * and potential values between processes.
 */
void densPotExchangeInit()
{
        int i;

        if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
                return;

        /* allocate communication buffers */
        sendDensPotData =
                (struct densPotData *) malloc(numExp * sizeof(struct densPotData));
#ifdef __MIC__
        recvDensPotData =
                (struct densPotData *) _mm_malloc(numImp * sizeof(struct densPotData),64);
#else
        recvDensPotData =
                (struct densPotData *) malloc(numImp * sizeof(struct densPotData));
#endif
        reqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);
        reqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);

        /* create density and potential buffer for exporting */
        for (i = 0; i < numExp; i++) {
                sendDensPotData[i].v = cells[expList[i].cell].v;
                sendDensPotData[i].density = cells[expList[i].cell].density;
        }

        /* send data - asynchronous MPI call */
        for (i = 0; i < MPIsize; i++) {
                if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                MPI_Isend(&sendDensPotData[sendOffset[i]],
                          sendCount[i] * sizeof(struct densPotData), MPI_BYTE, i,
                          MPIrank, MPI_COMM_WORLD, &reqSend[i]);
        }

        /* receive data - asynchronous MPI call */
        for (i = 0; i < MPIsize; i++) {
                if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                MPI_Irecv(&recvDensPotData[recvOffset[i]],
                          recvCount[i] * sizeof(struct densPotData), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &reqRecv[i]);
        }

}

/*!
 * This function waits for density and potential data exchange completion.
 */
void densPotExchangeWait()
{
        int i;
        MPI_Status status;

        if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
                return;

        // Wait for send completion
        for (i = 0; i < MPIsize; i++) {
                if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                if (MPI_Wait(&reqSend[i], &status) != MPI_SUCCESS)
                        stopRun(103, "sending", __FILE__, __LINE__);
        }

        // Wait for receive completion
        for (i = 0; i < MPIsize; i++) {
                if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
                        continue;
                if (MPI_Wait(&reqRecv[i], &status) != MPI_SUCCESS)
                        stopRun(103, "receiving", __FILE__, __LINE__);
        }

        /* some of the buffers can be deallocated */
        free(sendDensPotData);
        free(reqSend);
        free(reqRecv);
}
