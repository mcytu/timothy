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
#include "utils.h"

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
void createexportlist(system_t system,settings_t settings,cellsinfo_t cellsinfo,celltype_t* celltype,commdata_t *commdata)
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

                r = celltype[cellsinfo.cells[p].ctype].h * 1.5;
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
                        commdata->explist[commdata->numexp].cell = p;
                        commdata->explist[commdata->numexp].proc = procs[i];
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
        qsort(commdata->explist, commdata->numexp, sizeof(explist_t), explistcompare);
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
void exchangecleanup(system_t system,cellsinfo_t cellsinfo,commdata_t *commdata)
{
        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;
        free(commdata->recvcelldata);
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
void cellssendrecv(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata)
{
        int i;

        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;

        expcelldata_t *sendcelldata;
        expcelldata_t *recvcelldata;
        expdpdata_t *senddpdata;
        expdpdata_t *recvdpdata;
        /* allocate communication buffers */
        commdata->sendcelldata = (expcelldata_t *) malloc(commdata->numexp * sizeof(expcelldata_t));
        commdata->recvcelldata = (expcelldata_t *) malloc(commdata->numimp * sizeof(expcelldata_t));
        commdata->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * system.size);
        commdata->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * system.size);

        /* create reduced particle data buffer for exporting */
        for (i = 0; i < commdata->numexp; i++) {
                int64_t cellidx=commdata->explist[i].cell;
                commdata->sendcelldata[i].x = cellsinfo.cells[cellidx].x;
                commdata->sendcelldata[i].y = cellsinfo.cells[cellidx].y;
                commdata->sendcelldata[i].z = cellsinfo.cells[cellidx].z;
                commdata->sendcelldata[i].size = cellsinfo.cells[cellidx].size;
                commdata->sendcelldata[i].young = (double) cellsinfo.cells[cellidx].young;
                commdata->sendcelldata[i].ctype = cellsinfo.cells[cellidx].ctype;
        }

        /* send cells - asynchronous MPI call */
        for (i = 0; i < system.size; i++) {
                if (commdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Isend(&(commdata->sendcelldata[commdata->sendoffset[i]]),
                          commdata->sendcount[i] * sizeof(expcelldata_t), MPI_BYTE, i, system.rank,
                          MPI_COMM_WORLD, &(commdata->reqsend[i]));
        }

        /* receive cells - asynchronous MPI call */
        for (i = 0; i < system.size; i++) {
                if (commdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Irecv(&(commdata->recvcelldata[commdata->recvoffset[i]]),
                          commdata->recvcount[i] * sizeof(expcelldata_t), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &(commdata->reqrecv[i]));
        }

        return;
}

/*!
 * This function waits for cells' data exchange completion.
 */
void cellswait(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata)
{
        int i;
        MPI_Status status;

        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;

        /* wait for send completion */
        for (i = 0; i < system.size; i++) {
                if (commdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(commdata->reqsend[i]), &status) != MPI_SUCCESS)
                        terminate(system,"communication error", __FILE__, __LINE__);
        }

        /* wait for receive completion */
        for (i = 0; i < system.size; i++) {
                if (commdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(commdata->reqrecv[i]), &status) != MPI_SUCCESS)
                        terminate(system,"communication error", __FILE__, __LINE__);
        }

        /* some of the buffers can be deallocated here */
        free(commdata->sendcelldata);
        free(commdata->reqsend);
        free(commdata->reqrecv);
}

/*!
 * This function initiate sending and receiving density
 * and potential values between processes.
 */
void datasendrecv(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata)
{
        int i;

        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;

        /* allocate communication buffers */
        commdata->senddpdata =
                (expdpdata_t *) malloc(commdata->numexp * sizeof(expdpdata_t));
        commdata->recvdpdata =
                (expdpdata_t *) malloc(commdata->numimp * sizeof(expdpdata_t));
        commdata->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * system.size);
        commdata->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * system.size);

        /* create density and potential buffer for exporting */
        for (i = 0; i < commdata->numexp; i++) {
                int64_t cellidx=commdata->explist[i].cell;
                commdata->senddpdata[i].v = cellsinfo.cells[cellidx].v;
                commdata->senddpdata[i].density = cellsinfo.cells[cellidx].density;
        }

        /* send data - asynchronous MPI call */
        for (i = 0; i < system.size; i++) {
                if (commdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Isend(&(commdata->senddpdata[commdata->sendoffset[i]]),
                          commdata->sendcount[i] * sizeof(expdpdata_t), MPI_BYTE, i,
                          system.rank, MPI_COMM_WORLD, &(commdata->reqsend[i]));
        }

        /* receive data - asynchronous MPI call */
        for (i = 0; i < system.size; i++) {
                if (commdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Irecv(&(commdata->recvdpdata[commdata->recvoffset[i]]),
                          commdata->recvcount[i] * sizeof(expdpdata_t), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &(commdata->reqrecv[i]));
        }

        return;
}

/*!
 * This function waits for density and potential data exchange completion.
 */
void datawait(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata)
{
        int i;
        MPI_Status status;

        if (cellsinfo.globalcount.n < system.size*MIN_CELLS_PER_PROC || system.size == 1)
                return;

        // Wait for send completion
        for (i = 0; i < system.size; i++) {
                if (commdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(commdata->reqsend[i]), &status) != MPI_SUCCESS)
                        terminate(system,"communication error", __FILE__, __LINE__);
        }

        // Wait for receive completion
        for (i = 0; i < system.size; i++) {
                if (commdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(commdata->reqrecv[i]), &status) != MPI_SUCCESS)
                        terminate(system,"communication error", __FILE__, __LINE__);
        }
        /* some of the buffers can be deallocated */
        free(commdata->senddpdata);
        free(commdata->reqsend);
        free(commdata->reqrecv);
        return;
}
