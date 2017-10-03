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

/*! \file octree.c
 *  \brief contains functions that build octree and functions
 *         that provide a convenient way of operating on octree.
 */

/*!
 * This function creates an empty octree node in a given position.
 */
void octEmptyNode(int64_t father,int level,unsigned int xbit,unsigned int ybit,unsigned int zbit)
{
  int i;
  unsigned int bit;
  unsigned int len;
  if(father==-1) {
    octree[octSize].xcode=0;
    octree[octSize].ycode=0;
    octree[octSize].zcode=0;
    octree[octSize].xlimit=1<<level;
    octree[octSize].ylimit=1<<level;
    octree[octSize].zlimit=1<<level;
  } else {
    bit=1<<level;
    octree[octSize].xcode=octree[father].xcode|xbit;
    octree[octSize].ycode=octree[father].ycode|ybit;
    octree[octSize].zcode=octree[father].zcode|zbit;
    len=(octree[father].xlimit-octree[father].xcode)/2;
    octree[octSize].xlimit=octree[octSize].xcode+len;
    octree[octSize].ylimit=octree[octSize].ycode+len;
    octree[octSize].zlimit=octree[octSize].zcode+len;
  }
  octree[octSize].father=father;
  for(i=0; i<8; i++)
    octree[octSize].child[i]=-1;
  octree[octSize].data=-1;
  octree[octSize].level=level;
  octSize++;
}

/*!
 * This function inserts a given cell c into a proper position in octree.
 */
void octInsertCell(int64_t c)
{
  int64_t octIdx=0;
  int64_t father,child;
  unsigned int level;
  unsigned int childBranchBit,childIndex;
  unsigned int xbit,ybit,zbit;
  int inserted=0;

  level=ROOT_LEVEL-1;
  while(!inserted) {
    if(octree[octIdx].data==-1) { /* insert cell into node */
      octree[octIdx].data=c;
      inserted=1;
    } else {
      if(octree[octIdx].data==-2) { /* this is not a leaf */
        childBranchBit = 1 << (level);
        xbit=(locCode[c].x) & childBranchBit;
        ybit=(locCode[c].y) & childBranchBit;
        zbit=(locCode[c].z) & childBranchBit;
        childIndex = (   ((xbit) >> (level))
                         + ((ybit) >> (level-1))
                         + ((zbit) >> (level-2)) );
        level--;
        father=octIdx;
        child = octree[octIdx].child[childIndex];
        if(child==-1) {
          octree[octIdx].child[childIndex]=octSize;
          child=octSize;
          octEmptyNode(father,level,xbit,ybit,zbit);
        }
        octIdx=child;
      } else { /* occupied leaf */
        int64_t oc;
        oc=octree[octIdx].data;
        childBranchBit = 1 << (level);
        xbit=(locCode[oc].x) & childBranchBit;
        ybit=(locCode[oc].y) & childBranchBit;
        zbit=(locCode[oc].z) & childBranchBit;
        childIndex = (   ((xbit) >> (level))
                         + ((ybit) >> (level-1))
                         + ((zbit) >> (level-2)) );
        level--;
        father=octIdx;
        octree[octIdx].child[childIndex]=octSize;
        child=octSize;
        octEmptyNode(father,level,xbit,ybit,zbit);
        octree[child].data=octree[father].data;
        octree[father].data=-2;
        octIdx=father;
        level++;
      }
    }
  }
}

/*!
 * This is a driving function for creating an octree.
 */
void octBuild()
{
  int i,c;
  int64_t octMaxSize;
  struct doubleVector3d bMin,bMax;
  double epsilon=0.01;

  if(lnc==0)
    return;

#ifdef DEBUG
  printf("DEBUG: octBuild() start\n");
#endif

  bMin.x= DBL_MAX;
  bMin.y= DBL_MAX;
  bMin.z= DBL_MAX;
  bMax.x=-DBL_MAX;
  bMax.y=-DBL_MAX;
  bMax.z=-DBL_MAX;
  for(i=0; i<lnc; i++) {
    double x,y,z;
    double e;
    x=cells[i].x;
    y=cells[i].y;
    z=cells[i].z;
    e=h+epsilon;
    bMin.x=(x-e<bMin.x?x-e:bMin.x);
    bMax.x=(x+e>bMax.x?x+e:bMax.x);
    bMin.y=(y-e<bMin.y?y-e:bMin.y);
    bMax.y=(y+e>bMax.y?y+e:bMax.y);
    bMin.z=(z-e<bMin.z?z-e:bMin.z);
    bMax.z=(z+e>bMax.z?z+e:bMax.z);
  }
  affScale=bMax.x-bMin.x;
  affScale=(affScale>bMax.y-bMin.y?affScale:bMax.y-bMin.y);
  affScale=(affScale>bMax.z-bMin.z?affScale:bMax.z-bMin.z);
  affShift.x=bMin.x;
  affShift.y=bMin.y;
  affShift.z=bMin.z;
  /* each cell coordinate will be shifted by affShift and scaled by affScale */

  locCode=(struct uintVector3d*) malloc(sizeof(struct uintVector3d)*lnc);
  //#pragma omp parallel for
  for(c=0; c<lnc; c++) {
    locCode[c].x=(unsigned int)( ((cells[c].x-affShift.x)/affScale) * MAXVALUE );
    locCode[c].y=(unsigned int)( ((cells[c].y-affShift.y)/affScale) * MAXVALUE );
    locCode[c].z=(unsigned int)( ((cells[c].z-affShift.z)/affScale) * MAXVALUE );
  }
  /* memory space required to store octree (to be reviewed again!) */
  octMaxSize=lnc*64;
#ifdef __MIC__
  octree=_mm_malloc(sizeof(octNode)*octMaxSize,64);
#else
  octree=(octNode*) malloc(sizeof(octNode)*octMaxSize);
#endif
  root=0;
  octSize=0;
  octEmptyNode(-1,ROOT_LEVEL,0,0,0);

  for(c=0; c<lnc; c++) {
    octInsertCell(c);
  }

  /* a security check for space available in the octree buffer should be implemented */

#ifdef DEBUG
  printf("DEBUG: octBuild() end\n");
#endif

}

/*!
 * This function computes binary code for a given cell.
 */
void octComputeCode(int64_t c,struct uintVector3d *code)
{
  code[0].x=(unsigned int)( ((cells[c].x-affShift.x)/affScale) * MAXVALUE );
  code[0].y=(unsigned int)( ((cells[c].y-affShift.y)/affScale) * MAXVALUE );
  code[0].z=(unsigned int)( ((cells[c].z-affShift.z)/affScale) * MAXVALUE );
}

/*!
 * This function computes bounding box of a remote cell to be used for neighbour searching.
 * Difference between remote and local versions is the box croping which is applied
 * to some remote cells.
 */
MIC_ATTR void octComputeBoxR(int64_t c,struct uintVector3d *minLocCode,struct uintVector3d *maxLocCode)
{
  struct doubleVector3d minCor,maxCor;
  /* compute corners */
  minCor.x=(recvData[c].x-h-affShift.x)/affScale;
  minCor.y=(recvData[c].y-h-affShift.y)/affScale;
  minCor.z=(recvData[c].z-h-affShift.z)/affScale;
  maxCor.x=(recvData[c].x+h-affShift.x)/affScale;
  maxCor.y=(recvData[c].y+h-affShift.y)/affScale;
  maxCor.z=(recvData[c].z+h-affShift.z)/affScale;
  /* for remote cells - box crop */
  minCor.x=(minCor.x<0.0?0.0:minCor.x);
  minCor.y=(minCor.y<0.0?0.0:minCor.y);
  minCor.z=(minCor.z<0.0?0.0:minCor.z);
  maxCor.x=(maxCor.x>1.0?1.0:maxCor.x);
  maxCor.y=(maxCor.y>1.0?1.0:maxCor.y);
  maxCor.z=(maxCor.z>1.0?1.0:maxCor.z);
  /* compute location codes of corners */
  minLocCode[0].x=(unsigned int)(minCor.x*MAXVALUE);
  minLocCode[0].y=(unsigned int)(minCor.y*MAXVALUE);
  minLocCode[0].z=(unsigned int)(minCor.z*MAXVALUE);
  maxLocCode[0].x=(unsigned int)(maxCor.x*MAXVALUE);
  maxLocCode[0].y=(unsigned int)(maxCor.y*MAXVALUE);
  maxLocCode[0].z=(unsigned int)(maxCor.z*MAXVALUE);
}

/*!
 * This function computes bounding box of a local cell to be used for neighbour searching.
 */
MIC_ATTR void octComputeBox(int64_t c,struct uintVector3d *minLocCode,struct uintVector3d *maxLocCode)
{
  struct doubleVector3d minCor,maxCor;
  /* compute corners */
  minCor.x=(cells[c].x-h-affShift.x)/affScale;
  minCor.y=(cells[c].y-h-affShift.y)/affScale;
  minCor.z=(cells[c].z-h-affShift.z)/affScale;
  maxCor.x=(cells[c].x+h-affShift.x)/affScale;
  maxCor.y=(cells[c].y+h-affShift.y)/affScale;
  maxCor.z=(cells[c].z+h-affShift.z)/affScale;
  /* compute location codes of corners */
  minLocCode[0].x=(unsigned int)(minCor.x*MAXVALUE);
  minLocCode[0].y=(unsigned int)(minCor.y*MAXVALUE);
  minLocCode[0].z=(unsigned int)(minCor.z*MAXVALUE);
  maxLocCode[0].x=(unsigned int)(maxCor.x*MAXVALUE);
  maxLocCode[0].y=(unsigned int)(maxCor.y*MAXVALUE);
  maxLocCode[0].z=(unsigned int)(maxCor.z*MAXVALUE);
}

/*!
 * This function traverses the tree to a given level and returns node index.
 */
MIC_ATTR int octTraverseToLevel(unsigned int level,unsigned int xloc,unsigned int yloc,unsigned int zloc,unsigned int lmin)
{
  int cellIdx=0;
  int parent=0;
  unsigned int n;
  unsigned int childBranchBit,childIndex;
  n=(level)-(lmin)+1;
  while(n--) {
    childBranchBit = 1 << (level);
    childIndex = (   (((xloc) & childBranchBit) >> (level))
                     + (((yloc) & childBranchBit) >> (level-1))
                     + (((zloc) & childBranchBit) >> ((level-2))) );
    level--;
    parent=cellIdx;
    cellIdx=octree[cellIdx].child[childIndex];
    if(octree[cellIdx].data != -2) break; /* a leaf */
  }
  if(cellIdx==-1) cellIdx=parent;
  return cellIdx;
}

/*!
 * This function determines the level of the smallest cell containing a given region.
 */
MIC_ATTR int octLocateRegion(struct uintVector3d minLocCode,struct uintVector3d maxLocCode)
{
  unsigned int l1,l2,lmin;
  unsigned int level;
  int cellIdx;
  struct uintVector3d locDiff;
  /* XOR of location codes */
  locDiff.x=minLocCode.x ^ maxLocCode.x;
  locDiff.y=minLocCode.y ^ maxLocCode.y;
  locDiff.z=minLocCode.z ^ maxLocCode.z;
  /* Determining the level of the smallest cell containing the region */
  l1=ROOT_LEVEL;
  l2=ROOT_LEVEL;
  lmin=ROOT_LEVEL;
  while(!(locDiff.x & (1<<l1)) && l1) l1--;
  while(!(locDiff.y & (1<<l2)) && (l2>l1)) l2--;
  while(!(locDiff.z & (1<<lmin)) && (lmin>l2)) lmin--;
  lmin++;
  level=ROOT_LEVEL-1;
  cellIdx=octTraverseToLevel(level,minLocCode.x,minLocCode.y,minLocCode.z,lmin);
  return cellIdx;
}

/*!
 * This function deallocates memory used by octree.
 */
void octFree()
{
  if(lnc==0) return;
#ifdef __MIC__
  _mm_free(octree);
#else
  free(octree);
#endif
  free(locCode);
}

#pragma offload_attribute(push,target(mic))

/*!
 * This function initializes heap for tree traversal.
 */
void octHeapInit(octHeap *ttHeap)
{
  ttHeap->size=64;
  ttHeap->count=0;
  ttHeap->data=malloc(sizeof(int)*ttHeap->size);
}

/*!
 * This function adds an element to the heap.
 */
void octHeapPush(octHeap *ttHeap,int idx)
{
  if(ttHeap->count==ttHeap->size) {
    ttHeap->size+=64;
    ttHeap->data=realloc(ttHeap->data,sizeof(int)*ttHeap->size);
    printf("realloc again\n");
  }
  ttHeap->data[ttHeap->count]=idx;
  ttHeap->count+=1;
}

/*!
 * This function removes an element from the heap.
 */
int octHeapPop(octHeap *ttHeap)
{
  ttHeap->count-=1;
  return ttHeap->data[ttHeap->count];
}

/*!
 * This functione deallocates memory used by the heap.
 */
void octHeapFree(octHeap *ttHeap)
{
  free(ttHeap->data);
}

#pragma offload_attribute(pop)
