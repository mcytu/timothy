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
#include "potential.h"
#include "fields.h"
#include "inline.h"

/*! \file compute.c
 *  \brief contains main computational function called in each time step of the simulation
 */

/*!
 * This function calls all important simulation steps (cellular dynamics and global fields computations).
 */
int computestep(system_t system, cellsinfo_t *cellsinfo, celltype_t* celltype,commdata_t *commdata)
{
        int p;
        double sf;

        /* 0. Initialization */

        //initCellsToGridExchange();

        /* initiate asynchronous data transfers between processors */
        cellscomminit(system,*cellsinfo,commdata);

        /* 1. Compute potential for local cells */

        /* compute potential for local cells */
        computepotential(cellsinfo,celltype,*commdata);

        /* 2. Solve global fields */

        //if(step>0) {
        //        waitCellsToGridExchange();
        //        fieldsSolve();
        //}


        /* wait for data transfers to finish */
        cellscommwait(system,*cellsinfo,commdata);

        /* 3. Compute potential for remote cells */
        computeremotepotential(cellsinfo,celltype,*commdata);

        /* 4. Add chemotactic term to potential */

        /* add chemotaxis term to potential */

        /* 5. Compute gradient of the potential for local cells */
        /* initiate transfer of the density and potential data from remote cells */
        dpcomminit(system,*cellsinfo,commdata);
        /* compute gradient of the potential for local cells */
        computegradient(cellsinfo,celltype,*commdata);

        /* 6. Interpolate global fields and compute gradient */

        /* interpolate data */
        //interpolateFieldsToCells();
        /* compute gradient of global fields */
        //fieldGradient();

        /* 7. Compute gradient of the potential for remote cells */
        /* wait for density and potential data from remote cells */
        dpcommwait(system,*cellsinfo,commdata);

        /* compute gradient of the potential for remote cells */
        computeremotegradient(cellsinfo,celltype,*commdata);

        return 0;
}
