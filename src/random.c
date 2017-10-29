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
#include <sprng.h>

#include "global.h"

/*!\file random.c
 * \brief contains functions which handle the RNG
 */

/*!
 * This function initializes the RNG.
 */
void randomstreaminit(system_t *system,settings_t *settings)
{

        settings->rseed=time(NULL) + system->rank;
        //      int gtype; // RNG type: DEFAULT_RNG_TYPE, SPRNG_LFG, SPRNG_LCG, SPRNG_LCG64, SPRNG_CMRG, SPRNG_MLFG, SPRNG_PMLCG
        //      gtype = DEFAULT_RNG_TYPE;
        //      stream=init_sprng(gtype,system->rank, system->size, SEED, SPRNG_DEFAULT);
        //settings->rstream = init_sprng(gtype, system->rank, system->size, SEED, SPRNG_DEFAULT);
}

/*!
 * This function deallocates memory used by the RNG.
 */
void randomstreamfree(settings_t *settings)
{
        //free_sprng(settings->rstream);
}
