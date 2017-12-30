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

#include "utils.h"
#include "inline.h"
#include "fields.h"

/*! \file cells.c
 *  \brief contains functions which control current states and evolution of cells
 */

/*!
 * This function checks whether the cell p is outside the computational box.
 */
static inline int outsidethebox(cellsinfo_t *cellsinfo,int c)
{
        double x, y, z, r;

        x = cellsinfo->cells[c].x;
        y = cellsinfo->cells[c].y;
        z = cellsinfo->cells[c].z;
        r = cellsinfo->cells[c].size;

        if (x - r < -BOXSIZEX/2.0 || x + r > BOXSIZEX/2.0 )
                return 1;
        if (y - r < -BOXSIZEY/2.0 || y + r > BOXSIZEY/2.0 )
                return 1;
        if ( cellsinfo->dimension == 3 && (z - r < -BOXSIZEZ/2.0 || z + r > BOXSIZEZ/2.0 ) )
                return 1;

        return 0;
}

/* normal distribution (Box-Muller transform) */
void cellsrandominit(int nrandom,int ctype,systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo) {
        int i;
        double x1, x2, x3;
        double z1, z2, z3;
        double r1, r2;
        double l;
        double D=1.0;
        uint64_t idx;

        idx=cellsinfo->localcount.n-1;

        if (settings.dimension == 2)
                D = celltype[ctype].size*pow(8.0 * nrandom, 1.0 / 2.0);
        if (settings.dimension == 3)
                D = celltype[ctype].size*pow(8.0 * nrandom, 1.0 / 3.0);

        simTime = 0;
        for (i = 0; i < nrandom; i++) {

                r2 = 1.1;
                while (r2 >= 1.0) {
                        r1 = 1.1;
                        while (r1 == 0 || r1 >= 1.0) {
                                x1 = ((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0;
                                x2 = ((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0;
                                x3 = ((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0;
                                r1 = x1 * x1 + x2 * x2 + x3 * x3;
                        }
                        l = sqrt(-2 * log(r1) / r1);
                        z1 = x1 * l;
                        z2 = x2 * l;
                        z3 = x3 * l;

                        r2 = z1 * z1 + z2 * z2 + z3 * z3;
                }

                cellsinfo->cells[cellsinfo->localcount.n].x=z1 * D + cellsinfo->cells[idx].x;
                cellsinfo->cells[cellsinfo->localcount.n].y=z2 * D + cellsinfo->cells[idx].y;
                if(settings.dimension>2)
                        cellsinfo->cells[cellsinfo->localcount.n].z=z3 * D + cellsinfo->cells[idx].z;

                cellsinfo->cells[cellsinfo->localcount.n].ctype=ctype;
                cellsinfo->cells[cellsinfo->localcount.n].density=0.0;
                cellsinfo->cells[cellsinfo->localcount.n].size=celltype[ctype].size;
                cellsinfo->localcount.n+=1;
                cellsinfo->localcount.g0phase+=1;
                cellsinfo->localtypecount[ctype].n+=1;
                cellsinfo->localtypecount[ctype].g0phase+=1;
                cellsinfo->cells[cellsinfo->localcount.n].young=2100.0 + (float)rand_r(&(settings.rseed)) * 100.0;


                if(cellsinfo->localcount.n==settings.maxlocalcells) {
                        printf("warning: too many local cells, skipping rest of file %s\n",celltype[i].inputfile);
                        break;
                }

        }
        return;
}

/*!
 * This function implements mitosis of cells.
 */
void mitosis(int c, unsigned char *removecell,int64_t *removecount)
{

        double sc;
        double shift[3];

        if (lnc + 1 > maxCellsPerProc)
                stopRun(109, NULL, __FILE__, __LINE__);

        if(cells[c].scstage<nscstages-1) {
                if(sprng(stream)>sctprob[cells[c].scstage]) {
                        cells[lnc].scstage=cells[c].scstage+1;
                        nscinst[cells[lnc].scstage]+=1;
                } else {
                        cells[lnc].scstage=cells[c].scstage;
                        nscinst[cells[lnc].scstage]+=1;
                }
                if(sprng(stream)>sctprob[cells[c].scstage]) {
                        nscinst[cells[c].scstage]-=1;
                        cells[c].scstage=cells[c].scstage+1;
                        nscinst[cells[c].scstage]+=1;
                }
                /*if(sprng(stream)>sctprob[cells[c].scstage]) {
                   cells[lnc].scstage=cells[c].scstage+1;
                   } else {
                   cells[lnc].scstage=cells[c].scstage;
                   }*/
                //nscinst[cells[lnc].scstage]+=1;
        }

        sc = sqrt(velocity[c].x * velocity[c].x + velocity[c].y * velocity[c].y +
                  velocity[c].z * velocity[c].z);

        /* daughter cells are shifted away from the center of parent cell */
        if (sc > 0 && mitrand == 0) { /* direction of shift related to velocity vector */
                sc = cells[c].size / (2 * sc);
                shift[0] = sc * velocity[c].x;
                shift[1] = sc * velocity[c].y;
                if (sdim == 3)
                        shift[2] = sc * velocity[c].z;
                else
                        shift[2] = 0.0;
        } else { /* direction of shift chosen randomly */
                int accept = 0;
                while (accept == 0) {
                        shift[0] = sprng(stream) * 2.0 - 1.0;
                        shift[1] = sprng(stream) * 2.0 - 1.0;
                        if (sdim == 3)
                                shift[2] = sprng(stream) * 2.0 - 1.0;
                        else
                                shift[2] = 0.0;
                        sc = sqrt(pow(shift[0], 2) + pow(shift[1], 2) + pow(shift[2], 2));
                        if (sc == 0)
                                continue;
                        sc = cells[c].size / (2 * sc);
                        shift[0] = sc * shift[0];
                        shift[1] = sc * shift[1];
                        shift[2] = sc * shift[2];
                        accept = 1;
                }
        }
        /* 1st daughter cell position, size, type and age */
        cells[lnc].x = cells[c].x + shift[0];
        cells[lnc].y = cells[c].y + shift[1];
        cells[lnc].z = cells[c].z + shift[2];
        cells[lnc].size = pow(2.0, -(1.0 / 3.0)) * cells[c].size;;
        cells[lnc].tumor = cells[c].tumor;
        cells[lnc].age = cells[c].age + 1;

        /* 2nd daughter cell position, size, type and age */
        cells[c].x -= shift[0];
        cells[c].y -= shift[1];
        cells[c].z -= shift[2];
        cells[c].size = cells[lnc].size;;
        cells[c].age += 1;

        /* 2nd daughter cell cycle phases lenghts */
        if (cells[c].tumor == 1) {
                cells[c].g1 = cg1 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].g2 = cg2 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].s = cs * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].m = cm * (1 + (sprng(stream) * 2 - 1) * v);
        } else {
                cells[c].g1 = g1 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].g2 = g2 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].s = s * (1 + (sprng(stream) * 2 - 1) * v);
                cells[c].m = m * (1 + (sprng(stream) * 2 - 1) * v);
        }
        /* 1st daughter cell global ID */
        cells[lnc].gid =
                (unsigned long long int) MPIrank *(unsigned long long int)
                maxCellsPerProc + (unsigned long long int) lnc;

        /* 1st daughter cell parameters */
        cells[lnc].v = 0.0;
        cells[lnc].density = cells[c].density;
        cells[lnc].ctype = cells[c].ctype;
        cells[lnc].young = 2100.0 + sprng(stream) * 100.0;
        cells[lnc].halo = 0;
        cells[lnc].phase = 1;
        cells[lnc].death = 0;
        cells[lnc].phasetime = 0.0;
        /* 1st daughter cell cycle phases lenghts */
        if (cells[lnc].tumor == 1) {
                cells[lnc].g1 = cg1 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].g2 = cg2 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].s = cs * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].m = cm * (1 + (sprng(stream) * 2 - 1) * v);
        } else {
                cells[lnc].g1 = g1 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].g2 = g2 * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].s = s * (1 + (sprng(stream) * 2 - 1) * v);
                cells[lnc].m = m * (1 + (sprng(stream) * 2 - 1) * v);
        }

        /* update local cell counters */
        if (cells[lnc].tumor == 1)
                lcnc += 1;
        lnc = lnc + 1;
        lg1nc += 1;
        /* increment local ID */
        localID++;

}

/*!
 * This function finds locates cell closest to the center of mass of the systeminfo
 * and marks this cell as a cancer cell.
 */
void markMiddleCancerCell()
{
        int c;
        int middle = 0;
        double dist;
        struct {
                double val;
                int rank;
        } lmdist, gmdist;
        double center[3];
        double gcenter[3];

        /* each process computes its local center of mass */
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
        for (c = 0; c < lnc; c++) {
                center[0] += cells[c].x / nc;
                center[1] += cells[c].y / nc;
                center[2] += cells[c].z / nc;
        }

        /* MPI Reduce operation computes global center of mass */
        MPI_Allreduce(center, gcenter, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* intialization */
        lmdist.rank = MPIrank;
        lmdist.val = INT_MAX;

        /* each process finds local cell closest to the global center of mass */
        for (c = 0; c < lnc; c++) {
                dist =
                        sqrt((cells[c].x - gcenter[0]) * (cells[c].x - gcenter[0]) +
                             (cells[c].y - gcenter[1]) * (cells[c].y - gcenter[1]) +
                             (cells[c].z - gcenter[2]) * (cells[c].z - gcenter[2]));
                if (dist < lmdist.val) {
                        lmdist.val = dist;
                        middle = c;
                }
        }

        /* MPI_Allreduce locates the cell closest to the global center of mass */
        MPI_Allreduce(&lmdist, &gmdist, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                      MPI_COMM_WORLD);
        /* mark the found cell as cancer one */
        if (MPIrank == gmdist.rank) {
                cells[middle].tumor = 1;
                cells[middle].phase = 1;
                lg0nc--;
                lg1nc++;
                lcnc++;
        }

        /* indicate that there is a cancer cell in the systeminfo */
        cancer = 1;
}

/*!
 * This function dealocates all tables allocated during initialization of cell data
 */
void cellsdestroy()
{
        int f;
        free(tlnc);
#ifdef __MIC__
        _mm_free(cells);
#else
        free(cells);
#endif
        for (f = 0; f < NFIELDS; f++)
                free(cellFields[f]);
        free(cellFields);
#ifdef __MIC__
        _mm_free(velocity);
#else
        free(velocity);
#endif
}

/*!
 * This function removes a dead cell from the simulation.
 */
void cellsDeath(int lnc_old,unsigned char *removecell,int64_t removecount)
{
        int c, pos;

        pos = 0;
        for (c = 0; c < lnc; c++) {
                /* shift cells after dead cell removal */
                if (c >= lnc_old) {
                        cells[pos] = cells[c];
                        pos++;
                        continue;
                }
                if (c != pos && removecell[c] == 0)
                        cells[pos] = cells[c];
                if (removecell[c] == 0)
                        pos++;
                /* update cell counters */
                if (removecell[c] == 1) {
                        switch (cells[c].phase) {
                        case 0:
                                lg0nc--;
                                break;
                        case 1:
                                lg1nc--;
                                break;
                        case 2:
                                lsnc--;
                                break;
                        case 3:
                                lg2nc--;
                                break;
                        case 4:
                                lmnc--;
                                break;
                        }
                        if (cells[c].tumor == 1)
                                lcnc--;
                }
        }
        lnc -= removecount;
}

/*!
 * This function updates cell counters.
 */
void updateCellCounters()
{
        MPI_Allgather(&lnc, 1, MPI_INT64_T, tlnc, 1, MPI_INT64_T,
                      MPI_COMM_WORLD);
        MPI_Allreduce(localCellCount, totalCellCount, numberOfCounts,
                      MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(nscinst, gnscinst, nscstages,
                      MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localbc, &globalbc, 1,
                      MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
}

void updateChemotaxis()
{
        int c;
        for(c=0; c<lnc; c++) {


        }
}


/*!
 * This function updates cells' positions.
 */
void updatepositions(settings_t settings,cellsinfo_t *cellsinfo,unsigned char *removecell, int64_t *removecount)
{
        int c;
        double alpha=0.00001;

        /* move cells */
        for (c = 0; c < cellsinfo->localcount.n; c++) {
                cellsinfo->cells[c].x += cellsinfo->forces[c].x;
                cellsinfo->cells[c].y += cellsinfo->forces[c].y;
                cellsinfo->cells[c].z += cellsinfo->forces[c].z;
                /* random movement */
                cellsinfo->cells[c].x += alpha*(((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0);
                cellsinfo->cells[c].y += alpha*(((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0);
                cellsinfo->cells[c].z += alpha*(((float)rand_r(&(settings.rseed))/RAND_MAX)*2.0 - 1.0);

                /* mark cells outside the box for removal */
                if( outsidethebox(cellsinfo,c) ) {
                        removecell[c]=1;
                        (*removecount)=(*removecount)+1;
                }
        }
        return;
}

/*!
 * This function updates cells' cycle phases.
 */
int updateCellCycles(cellsinfo_t *cellsinfo,unsigned char *removecell,int64_t *removecount)
{

        int c;
        double eps, epsCancer;
        int lncAtThisStep;

        eps = densityCriticalLevel1;
        epsCancer = densityCriticalLevel2;

        lncAtThisStep = lnc;

        for (c = 0; c < lncAtThisStep; c++) {

                if(cells[c].ctype==1) continue;

                if ( outsidethebox(cellsinfo,c) ) {
                        removecell[c] = 1;
                        (*removecount)=(*removecell)+1;
                        continue;
                }

                if (removecell[c])
                        continue;

                if (simStart) {

                        if (cells[c].phase != 0
                            && ((cells[c].tumor == 0 && cells[c].density <= eps)
                                || (cells[c].tumor == 1 && cells[c].density <= epsCancer)))
                                cells[c].phasetime += gfDt / 3600.0;

                        switch (cells[c].phase) {

                        case 0: /* G0 phase */
                                if (gfields && oxygen
                                    && cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
                                        cells[c].phase = 5;
                                        cells[c].phasetime = 0;
                                        lg0nc--;
                                        lnnc++;
                                        break;
                                }
                                /* transition to G1 phase */
                                if ((cells[c].tumor == 0 && cells[c].density <= eps) || /* enough space for healthy cell */
                                    (cells[c].tumor == 1 && cells[c].density <= epsCancer) || /* enough space for tumor cell */
                                    nc == 1 || /* only single cell in the simulation */
                                    (gfields && oxygen && cellFields[OXYG][c] >= fieldCriticalLevel1[OXYG])) { /* sufficient level of oxygen */
                                        cells[c].phase = 1;
                                        lg0nc--;
                                        lg1nc++;
                                        break;
                                }
                                break;
                        case 1: /* G1 phase */
                                /* transition to G0 or Necrotic phase */
                                if ((cells[c].tumor == 0 && cells[c].density > eps) || /* too crowdy for healthy cell */
                                    (cells[c].tumor == 1 && cells[c].density > epsCancer) || /* too crowdy for tumor cell */
                                    (gfields && oxygen && cellFields[OXYG][c] < fieldCriticalLevel1[OXYG])) { /* too low oxygen level */
                                        if (gfields && oxygen && cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) { /* transition to Necrotic phase */
                                                cells[c].phase = 5;
                                                cells[c].phasetime = 0;
                                                lg1nc--;
                                                lnnc++;
                                        } else { /* transition to G0 phase */
                                                cells[c].phase = 0;
                                                lg1nc--;
                                                lg0nc++;
                                        }
                                        break;
                                }
                                /* cells grow in phase G1 */
                                if (cells[c].size < csize) {
                                        cells[c].size +=
                                                (csize -
                                                 pow(2.0,
                                                     -(1.0 / 3.0)) * csize) * (gfDt) / (3600.0 *
                                                                                        cells[c].g1);
                                }
                                if (cells[c].size > csize)
                                        cells[c].size = csize;
                                if (cells[c].phasetime >= cells[c].g1) {
                                        int death;
                                        cells[c].phase = 2;
                                        cells[c].phasetime = 0;
                                        lg1nc--;
                                        lsnc++;
                                        if (cells[c].tumor == 0) {
                                                death = (sprng(stream) < rd ? 1 : 0);
                                                if (death) {
                                                        removecell[c] = 1;
                                                        (*removecount)=(*removecount)+1;
                                                }
                                        }
                                }
                                break;
                        case 2: /* S phase */
                                if (gfields && oxygen
                                    && cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
                                        cells[c].phase = 5;
                                        cells[c].phasetime = 0;
                                        lsnc--;
                                        lnnc++;
                                        break;
                                }
                                if (cells[c].phasetime >= cells[c].s) {
                                        cells[c].phase = 3;
                                        cells[c].phasetime = 0;
                                        lsnc--;
                                        lg2nc++;
                                        break;
                                }
                                break;
                        case 3: /* G2 phase */
                                if (gfields && oxygen
                                    && cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
                                        cells[c].phase = 5;
                                        cells[c].phasetime = 0;
                                        lg2nc--;
                                        lnnc++;
                                        break;
                                }
                                if (cells[c].phasetime >= cells[c].g2) {
                                        int death;
                                        cells[c].phase = 4;
                                        cells[c].phasetime = 0;
                                        lg2nc--;
                                        lmnc++;
                                        if (cells[c].tumor == 0) {
                                                death = (sprng(stream) < rd ? 1 : 0);
                                                if (death) {
                                                        removecell[c] = 1;
                                                        (*removecount)=(*removecount)+1;
                                                }
                                        }
                                        break;
                                }
                                break;
                        case 4: /* M phase */
                                if (gfields && oxygen
                                    && cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
                                        cells[c].phase = 5;
                                        cells[c].phasetime = 0;
                                        lmnc--;
                                        lnnc++;

                                } else if (cells[c].phasetime >= cells[c].m) {
                                        mitosis(c,removecell,removecount);
                                        cells[c].phase = 1;
                                        cells[c].phasetime = 0;
                                        lmnc--;
                                        lg1nc++;
                                }
                                break;
                        } // switch
                } // if
        } // for loop

        /* update global number of cells */
        MPI_Allreduce(&lnc, &nc, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

        return 0;
}

/*!
 * This function fills the scalarOutput field of each cell.
 * It can be modified to output any float scalars that
 * user would like to analyze or visualize after simulation.
 * This field is printed to the output VTK files.
 */
void additionalScalarField()
{
        int c;
        for (c = 0; c < lnc; c++) {
                if (cells[c].tumor == 1)
                        cells[c].scalarField = 8.0;
                else
                        cells[c].scalarField = cells[c].density;
        }
}

/*!
 * This function drives the whole cell cycle update.
 */
void cellsupdate(settings_t settings,cellsinfo_t *cellsinfo)
{
        int lnc_old;
        unsigned char *removecell;
        int64_t removecount;
        /* number of local cells might change during the update */
        lnc_old = lnc;
        removecell = (unsigned char *) calloc(lnc_old, sizeof(unsigned char));
        removecount = 0;

        updatepositions(settings,cellsinfo,removecell,&removecount);
/*
        updateCellCycles();
        if (nhs > 0 && nc > nhs && tgs == 1 && cancer == 0)
                markMiddleCancerCell();
        if (nhs > 0 && nc > nhs)
                cellsDeath(lnc_old);
        updateCellCounters();
        additionalScalarField();*/
        free(removecell);
}
