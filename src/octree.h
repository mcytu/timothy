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

MIC_ATTR void octHeapInit(octHeap *ttHeap);
MIC_ATTR void octComputeBox(int64_t c,struct uintVector3d *minLocCode,struct uintVector3d *maxLocCode);
MIC_ATTR int octLocateRegion(struct uintVector3d minLocCode,struct uintVector3d maxLocCode);
MIC_ATTR void octHeapPush(octHeap *ttHeap,int idx);
MIC_ATTR int octHeapPop(octHeap *ttHeap);
MIC_ATTR void octHeapFree(octHeap *ttHeap);
MIC_ATTR static inline int octNodeIntersection(int idx,struct uintVector3d minLocCode,struct uintVector3d maxLocCode);
MIC_ATTR void octComputeBoxR(int64_t c,struct uintVector3d *minLocCode,struct uintVector3d *maxLocCode);
