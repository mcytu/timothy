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
#include "octree.h"
#include "fields.h"
#include "inline.h"

/*! \file potential.c
 *  \brief contains functions that compute the potential
 */

/*!
 * This function computes potential for two neighbour cells.
 * mode=0 - computations for two local cells
 * mode=1 - computations for one local and one remote cell
 */
MIC_ATTR double ccpot(int p1, int p2, int mode,double *mindist,cellsinfo_t *cellsinfo,commdata_t commdata)
{
        double pot;
        double dist=0.0;
        double size;
        double x, y, z, h;
        double xc;
        double D;
        double poisson = 0.33;
        double young;
        int ctype;

        if (mode == 0 && p1 == p2)
                return 0.0;

        if (mode == 0) {
                x = cellsinfo->cells[p1].x;
                y = cellsinfo->cells[p1].y;
                z = cellsinfo->cells[p1].z;
                size = cellsinfo->cells[p1].size;
                young = cellsinfo->cells[p1].young;
                ctype = cellsinfo->cells[p1].ctype;
                h = cellsinfo->cells[p1].h;
        }
        if (mode == 1) {
                x = commdata.recvcelldata[p1].x;
                y = commdata.recvcelldata[p1].y;
                z = commdata.recvcelldata[p1].z;
                size = commdata.recvcelldata[p1].size;
                young = commdata.recvcelldata[p1].young;
                ctype = commdata.recvcelldata[p1].ctype;
                h = commdata.recvcelldata[p1].h;
        }

        // this should be replaced by check (static cells or moving cells)
        if(ctype==1 && cellsinfo->cells[p2].ctype==1) return 0.0;

        /* compute the distance between two cells */
        if (cellsinfo->dimension == 2)
                dist =
                        sqrt((x - cellsinfo->cells[p2].x) * (x - cellsinfo->cells[p2].x) +
                             (y - cellsinfo->cells[p2].y) * (y - cellsinfo->cells[p2].y));
        if (cellsinfo->dimension == 3)
                dist =
                        sqrt((x - cellsinfo->cells[p2].x) * (x - cellsinfo->cells[p2].x) +
                             (y - cellsinfo->cells[p2].y) * (y - cellsinfo->cells[p2].y) +
                             (z - cellsinfo->cells[p2].z) * (z - cellsinfo->cells[p2].z));

        if (mindist[0] > dist ) {
                mindist[0] = dist;
        }

        if (dist <= h) {
//  if (dist <= size + cells[p2].size) {

                double r01, r02;
                double area;
                double sc=1.0;

                if (cellsinfo->dimension == 2)
                        sc = h*h;
                if (cellsinfo->dimension == 3)
                        sc = h*h*h;

                if (mode == 0) {
                        cellsinfo->cells[p1].density +=
                                sc * (cellsinfo->cells[p2].size / csize) * sph_kernel(dist);
                }
                if (mode == 1) {
      #pragma omp atomic
                        cellsinfo->cells[p2].density += sc * (size / csize) * sph_kernel(dist);
                }

                xc = size + cellsinfo->cells[p2].size - dist;

                if (xc <= 0.0)
                        return 0.0;

                D = 0.75 * ((1.0 - poisson * poisson) / young +
                            (1.0 - poisson * poisson / cellsinfo->cells[p2].young));

                /* adhesion */
                r01 =
                        (size * size - cellsinfo->cells[p2].size * cellsinfo->cells[p2].size +
                         dist * dist) / (2 * dist);
                r02 = dist - r01;

                area =
                        M_PI *
                        ((size * size * (size - r01) -
                          (size * size * size - r01 * r01 * r01) / 3) +
                         (cellsinfo->cells[p2].size * cellsinfo->cells[p2].size * (cellsinfo->cells[p2].size - r02) -
                          (cellsinfo->cells[p2].size * cellsinfo->cells[p2].size * cellsinfo->cells[p2].size -
                           r02 * r02 * r02) / 3));

                /* compute potential */
                pot =
                        (2.0 * pow(xc, 5 / 2) / (5.0 * D)) * sqrt((size * cellsinfo->cells[p2].size) /
                                                                  (size +
                                                                   cellsinfo->cells[p2].size)) +
                        area * 0.1;

                return pot;

        } else
                return 0.0;

}

/*!
 * This function implements tree walk algorithm for each local cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
void computepotential(cellsinfo_t *cellsinfo,commdata_t commdata)
{
        if(cellsinfo->localcount.n<=1) return;
  #pragma omp parallel
        {
                int p;
                double mindist;
                mindist=nx;

    #pragma omp for schedule(dynamic,64)
                for (p = 0; p < cellsinfo->localcount.n; p++) {
                        int64_t cellIdx;
                        int newIdx;
                        int s;
                        int octIndex;
                        uint3dv_t minLocCode,maxLocCode;
                        octheap_t octh;

                        octheapinit(&octh);

                        cellsinfo->cells[p].density = 0.0;
                        cellsinfo->cells[p].v = 0.0;

                        octcomputebox(p,&minLocCode,&maxLocCode,*cellsinfo);
                        octIndex=octlocateregion(minLocCode,maxLocCode,*cellsinfo);

                        for(s=0; s<8; s++) {
                                int idx;
                                idx=cellsinfo->octree[octIndex].child[s];
                                if(idx!=-1 && octnodeintersection(idx,minLocCode,maxLocCode,*cellsinfo))
                                        octheappush(&octh,idx);
                        }

                        while(octh.count>0) {
                                int idx;
                                idx=octheappop(&octh);
                                cellIdx=cellsinfo->octree[idx].data;
                                if(cellIdx>=0) {
                                        cellsinfo->cells[p].v += ccpot(p,cellIdx,0,&mindist,cellsinfo,commdata);
                                } else {
                                        for(s=0; s<8; s++) {
                                                newIdx=cellsinfo->octree[idx].child[s];
                                                if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                        octheappush(&octh,newIdx);
                                        }
                                }
                        }
                        octheapfree(&octh);
                }
    #pragma omp critical
                statistics.mindist=(statistics.mindist>mindist ? mindist : statistics.mindist);
        }
        return;
}

/*!
 * This function implements tree walk algorithm for each remote cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
void computeremotepotential(cellsinfo_t *cellsinfo,commdata_t commdata)
{
        int rp;
        if(commdata.numimp<=0) return;
  #pragma omp parallel
        {
                double mindist;
                mindist=nx;
    #pragma omp for schedule(dynamic,64)
                for (rp = 0; rp < commdata.numimp; rp++) {
                        int64_t cellIdx;
                        int newIdx;
                        int s;
                        double v;
                        uint3dv_t minLocCode,maxLocCode;
                        octheap_t octh;
                        octheapinit(&octh);
                        octcomputeboxr(rp,&minLocCode,&maxLocCode,commdata);
                        octheappush(&octh,0);
                        while(octh.count>0) {
                                int idx;
                                idx=octheappop(&octh);
                                cellIdx=cellsinfo->octree[idx].data;
                                if(cellIdx>=0) {
                                        v=ccpot(rp,cellIdx,1,&mindist,cellsinfo,commdata);
          #pragma omp atomic
                                        cellsinfo->cells[cellIdx].v += v;
                                } else {
                                        for(s=0; s<8; s++) {
                                                newIdx=cellsinfo->octree[idx].child[s];
                                                if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                        octheappush(&octh,newIdx);
                                        }
                                }
                        }
                        octheapfree(&octh);
                }
    #pragma omp critical
                statistics.mindist=(statistics.mindist>mindist ? mindist : statistics.mindist);
        }
}

/*!
 *  This function computes potential gradient for two neighbour cells.
 * mode=0 - computations for local cells
 * mode=1 - computations for remote cells
 */
void ccpotgrad(int p1, int p2, int mode,cellsinfo_t *cellsinfo,commdata_t commdata)
{
        double grad[3];
        /* we assume that cells' mass is always 1.0 */
        double m = 1.0;
        double sc;
        double x1, x2, y1, y2, z1, z2;
        double v, density, size;
        double dist;
        int ctype;

        if (p1 == p2 && mode == 0)
                return;

        if (mode == 0) {
                x1 = cellsinfo->cells[p1].x;
                x2 = cellsinfo->cells[p2].x;
                y1 = cellsinfo->cells[p1].y;
                y2 = cellsinfo->cells[p2].y;
                z1 = cellsinfo->cells[p1].z;
                z2 = cellsinfo->cells[p2].z;
                v = cellsinfo->cells[p2].v;
                density = cellsinfo->cells[p2].density;
                size = cellsinfo->cells[p2].size;
                ctype = cellsinfo->cells[p1].ctype;
        }
        if (mode == 1) {
                x1 = commdata.recvcelldata[p1].x;
                x2 = cellsinfo->cells[p2].x;
                y1 = commdata.recvcelldata[p1].y;
                y2 = cellsinfo->cells[p2].y;
                z1 = commdata.recvcelldata[p1].z;
                z2 = cellsinfo->cells[p2].z;
                v = commdata.recvdpdata[p1].v;
                density = commdata.recvdpdata[p1].density;
                size = commdata.recvcelldata[p1].size;
                ctype = commdata.recvcelldata[p1].ctype;
        }

        if (sdim == 2)
                dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        if (sdim == 3)
                dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) +
                            (z1 - z2) * (z1 - z2));

        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;

        /* compute the gradient of SPH kernel function */
        sph_kernel_gradient(p1, p2, grad, mode, dist, *cellsinfo, commdata);

        if (density == 0.0)
                sc = 0.0;
        else {
                m = size / csize;
                sc = (m / density) * v;
        }

        /* update forces */
        if (mode == 0) {
                velocity[p1].x += sc * grad[0];
                velocity[p1].y += sc * grad[1];
                velocity[p1].z += sc * grad[2];
        } else {
 #pragma omp atomic
                velocity[p2].x -= sc * grad[0];
 #pragma omp atomic
                velocity[p2].y -= sc * grad[1];
 #pragma omp atomic
                velocity[p2].z -= sc * grad[2];
        }
}

/*!
 * This function implements tree walk algorithm for each local cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
void computegradient(cellsinfo_t *cellsinfo,commdata_t commdata)
{
        int p;
        if(lnc<=1) return;
  #pragma omp parallel for schedule(dynamic,64)
        for(p=0; p<lnc; p++) {
                int64_t cellIdx;
                int newIdx;
                int s;
                int octIndex;
                uint3dv_t minLocCode,maxLocCode;
                octheap_t octh;

                octheapinit(&octh);
                velocity[p].x=0.0;
                velocity[p].y=0.0;
                velocity[p].z=0.0;

                //if(cells[p].ctype==1) continue;

                octcomputebox(p,&minLocCode,&maxLocCode,*cellsinfo);
                octIndex=octlocateregion(minLocCode,maxLocCode,*cellsinfo);

                for(s=0; s<8; s++) {
                        int idx;
                        idx=cellsinfo->octree[octIndex].child[s];
                        if(idx!=-1 && octnodeintersection(idx,minLocCode,maxLocCode,*cellsinfo))
                                octheappush(&octh,idx);
                }

                while(octh.count>0) {
                        int idx;
                        idx=octheappop(&octh);
                        cellIdx=cellsinfo->octree[idx].data;
                        if(cellIdx>=0) {
                                ccpotgrad(p,cellIdx,0,cellsinfo,commdata);
                        } else {
                                for(s=0; s<8; s++) {
                                        newIdx=cellsinfo->octree[idx].child[s];
                                        if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                octheappush(&octh,newIdx);
                                }
                        }
                }
                octheapfree(&octh);
        }
        return;
}

/*!
 * This function implements tree walk algorithm for each remote cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
void computeremotegradient(cellsinfo_t *cellsinfo,commdata_t commdata)
{
        int rp;
        if(commdata.numimp<=0) return;
  #pragma omp parallel for schedule(dynamic,64)
        for (rp = 0; rp < commdata.numimp; rp++) {
                int64_t cellIdx;
                int newIdx;
                int s;
                uint3dv_t minLocCode,maxLocCode;
                octheap_t octh;
                octheapinit(&octh);
                octcomputeboxr(rp,&minLocCode,&maxLocCode,commdata);
                octheappush(&octh,0);
                while(octh.count>0) {
                        int idx;
                        idx=octheappop(&octh);
                        cellIdx=cellsinfo->octree[idx].data;
                        if(cellIdx>=0) {
                                ccpotgrad(rp,cellIdx,1,cellsinfo,commdata);
                        } else {
                                for(s=0; s<8; s++) {
                                        newIdx=cellsinfo->octree[idx].child[s];
                                        if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                octheappush(&octh,newIdx);
                                }
                        }
                }
                octheapfree(&octh);
        }
        return;
}
