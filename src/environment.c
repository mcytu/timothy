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

#include <float.h>
#include <inttypes.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "environment.h"
#include "mpi.h"

/* zmienne lokalne w pliku environment */
HYPRE_SStructGrid envGrid;
HYPRE_SStructGraph envGraph;
HYPRE_SStructStencil *envStencil;
HYPRE_SStructMatrix A;
HYPRE_SStructVector b;
HYPRE_SStructVector x;
HYPRE_ParCSRMatrix parA;
HYPRE_ParVector parb;
HYPRE_ParVector parx;
HYPRE_Solver envSolver;
HYPRE_Solver envPrecond;
int envObjectType;
HYPRE_Int envLower[3], envUpper[3];
HYPRE_Int bcLower[3];
HYPRE_Int bcUpper[3];
double *dt;
int envIter;
double *envZ;
HYPRE_SStructVariable *envVarTypes;
int numberOfIters;
/* koniec */

/* do polaczenia z fields */
double *fieldDt;
double **fieldAddr;
double **envField;
double *fieldConsumption; /* units - mol (cell s)^-1 */
double *fieldProduction;  /* units - mol (cell s)^-1 */
/* koniec */

/* tych zmiennych nie odnalazlem w strukturach */
double csize=0.2735;
double csizeInUnits=10.0;
int bvsim=0;
/* koniec */


void envalloc(systeminfo_t systeminfo,settings_t settings,grid_t grid,environment_t **environment,solverdata_t *solverdata) {
        int i;
#ifdef HYPRE
        solverdata->z=(double*)calloc(settings.numberoffields,sizeof(double));
        solverdata->vartypes=(HYPRE_SStructVariable*)calloc(settings.numberoffields,sizeof(HYPRE_SStructVariable));
        solverdata->stencil=(HYPRE_SStructStencil*)calloc(settings.numberoffields,sizeof(HYPRE_SStructStencil));
#endif
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].data=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->data", __FILE__, __LINE__);
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].production=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->production", __FILE__, __LINE__);
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].gradient=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z*3,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->gradient", __FILE__, __LINE__);
        return;
}

void envinit(systeminfo_t systeminfo,settings_t settings,grid_t grid,environment_t **environment) {
        int i,j,k,f;
        int yz=grid.localsize.y * grid.localsize.z;
        for(f=0; f<settings.numberoffields; f++) {
                for (i = 0; i < grid.localsize.x; i++)
                        for (j = 0; j < grid.localsize.y; j++)
                                for (k = 0; k < grid.localsize.z; k++) {
                                        (*environment)[f].data[yz * i + grid.localsize.z * j + k] =
                                                (*environment)[f].initialconditionmean + ((double)(rand())/RAND_MAX)*(*environment)[f].initialconditionvariance;
                                }
        }
        return;
}



/*!
 * This function sets boundary conditions for domain faces.
 */
void envSetBoundary(int coord, int boundary)
{
        if (coord == 0 && boundary == -1) {
                bcLower[0] = envLower[0];
                bcUpper[0] = envLower[0];
                bcLower[1] = envLower[1];
                bcUpper[1] = envUpper[1];
                bcLower[2] = envLower[2];
                bcUpper[2] = envUpper[2];
        }
        if (coord == 0 && boundary == 1) {
                bcLower[0] = envUpper[0];
                bcUpper[0] = envUpper[0];
                bcLower[1] = envLower[1];
                bcUpper[1] = envUpper[1];
                bcLower[2] = envLower[2];
                bcUpper[2] = envUpper[2];
        }
        if (coord == 1 && boundary == -1) {
                bcLower[0] = envLower[0];
                bcUpper[0] = envUpper[0];
                bcLower[1] = envLower[1];
                bcUpper[1] = envLower[1];
                bcLower[2] = envLower[2];
                bcUpper[2] = envUpper[2];
        }
        if (coord == 1 && boundary == 1) {
                bcLower[0] = envLower[0];
                bcUpper[0] = envUpper[0];
                bcLower[1] = envUpper[1];
                bcUpper[1] = envUpper[1];
                bcLower[2] = envLower[2];
                bcUpper[2] = envUpper[2];
        }
        if (coord == 2 && boundary == -1) {
                bcLower[0] = envLower[0];
                bcUpper[0] = envUpper[0];
                bcLower[1] = envLower[1];
                bcUpper[1] = envUpper[1];
                bcLower[2] = envLower[2];
                bcUpper[2] = envLower[2];
        }
        if (coord == 2 && boundary == 1) {
                bcLower[0] = envLower[0];
                bcUpper[0] = envUpper[0];
                bcLower[1] = envLower[1];
                bcUpper[1] = envUpper[1];
                bcLower[2] = envUpper[2];
                bcUpper[2] = envUpper[2];
        }
}

void envInit(int nnutrients) {
        envZ=(double*)malloc(nnutrients*sizeof(double));
        envStencil=(HYPRE_SStructStencil*)malloc(nnutrients*sizeof(HYPRE_SStructStencil));
        envVarTypes=(HYPRE_SStructVariable*)malloc(nnutrients*sizeof(HYPRE_SStructVariable));
}

/*!
 * This function initializes grid, stencil and matrix for a given envical field.
 */
void envinitsystem(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment)
{
        int i, j, k, c;
        int entry;
        int var;
        HYPRE_Int offsets[7][3] = {
                {0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                {0, 1, 0}, { 0, 0,-1}, {0, 0, 1}
        };

        double gridResolutionInUnits; /* grid resolution in centimeters */

        numberOfIters = 1;
        envIter=0;
        for(var=0; var<nnutrients; var++) {
                envVarTypes[var]=HYPRE_SSTRUCT_VARIABLE_NODE;
                dt[var]=fieldDt[var];
        }

        /* 1. INIT GRID */

        /* create an empty 3D grid object */
        HYPRE_SStructGridCreate(systeminfo.MPI_CART_COMM, 3, 1, &envGrid);

        /* set this process box */
        envLower[0] = grid->loweridx[systeminfo.rank].x;
        envLower[1] = grid->loweridx[systeminfo.rank].y;
        envLower[2] = grid->loweridx[systeminfo.rank].z;

        envUpper[0] = grid->upperidx[systeminfo.rank].x;
        envUpper[1] = grid->upperidx[systeminfo.rank].y;
        envUpper[2] = grid->upperidx[systeminfo.rank].z;

        /* add a new box to the grid */
        HYPRE_SStructGridSetExtents(envGrid, 0, envLower, envUpper);

        HYPRE_SStructGridSetVariables(envGrid, 0, nnutrients, envVarTypes);
        HYPRE_SStructGridAssemble(envGrid);

        //  2. INIT STENCIL
        //HYPRE_SStructStencilCreate(3, 7, &envStencil);
        for(var = 0; var < settings.numberoffields; var++) {
                HYPRE_SStructStencilCreate(3, 7, &envStencil[var]);
                for(entry = 0; entry < 7; entry++)
                        HYPRE_SStructStencilSetEntry(envStencil[var], entry, offsets[entry],var);

        }

        // 3. SET UP THE GRAPH
        // assumption - all stencils are the same
        envObjectType = HYPRE_PARCSR;
        HYPRE_SStructGraphCreate(systeminfo.MPI_CART_COMM, envGrid, &envGraph);
        HYPRE_SStructGraphSetObjectType(envGraph, envObjectType);
        for(var=0; var<settings.numberoffields; var++)
                HYPRE_SStructGraphSetStencil(envGraph, 0, var, envStencil[var]);
        HYPRE_SStructGraphAssemble(envGraph);

        // 4. SET UP MATRIX
        long long nentries = 7;
        long long nvalues;
        double *values;
        HYPRE_Int stencil_indices[7];

        nvalues = nentries * grid->localsize.x * grid->localsize.y * grid->localsize.z;
        // create an empty matrix object
        HYPRE_SStructMatrixCreate(systeminfo.MPI_CART_COMM, envGraph, &A);
        HYPRE_SStructMatrixSetObjectType(A, envObjectType);
        // indicate that the matrix coefficients are ready to be set
        HYPRE_SStructMatrixInitialize(A);

        values = calloc(nvalues, sizeof(double));

        for (j = 0; j < nentries; j++)
                stencil_indices[j] = j;

        gridResolutionInUnits = (grid->resolution / csize) * csizeInUnits * 0.0001;

        for(var=0; var<settings.numberoffields; var++) {

                envZ[var] =
                        nutrient[var].diffusion_coefficient * dt[var] / (gridResolutionInUnits *
                                                                         gridResolutionInUnits);

                // set the standard stencil at each grid point,
                //  we will fix the boundaries later
                for (i = 0; i < nvalues; i += nentries) {
                        values[i] = 1 + 6.0 * envZ[var] + nutrient[var].lambda_delay * dt[var];
                        for (j = 1; j < nentries; j++)
                                values[i + j] = -envZ[var];
                }

                HYPRE_SStructMatrixSetBoxValues(A, 0, envLower, envUpper, var,
                                                nentries, stencil_indices, values);
        }

        free(values);
}

/*!
 * This function computes cell production/consumption function based on
 * the interpolated cell density field.
 */
void envCellPC(struct state simstate,double *envPC, int var, struct gridData grid)
{
        int ch __attribute__((unused)), i, j, k;

        if (simstate.step == 0)
                return;

        int idx = 0;
        for (k = 0; k < grid.local_size.z; k++)
                for (j = 0; j < grid.local_size.y; j++)
                        for (i = 0; i < grid.local_size.x; i++, idx++) {
                                envPC[idx] = -fieldConsumption[var] * 0.1 * dt[var];
                                /*envPC[idx] = -fieldConsumption[nch] * tissueField[gridSize.z * gridSize.y * i + gridSize.z * j + k] * dt[nch];
                                   if(bvsim) envPC[idx]+=fieldProduction[nch] * vesselField[gridSize.z * gridSize.y * i + gridSize.z * j +k] * dt[nch] ;	*/
                                //*(cellVolume/boxVolume);//*(1.0/cellVolume);//*dt[nch];//*dt[nch];
                        }
}

/*!
 * This function initializes boundary conditions for a given envical field.
 */
void envInitBC(struct state simstate,int nnutrients,struct environment *nutrient,struct gridData grid)
{
        int i, j, k;
        int mi;
        int var;
        int nentries = 1;
        HYPRE_Int stencil_indices[1];
        long long nvalues = grid.local_size.x * grid.local_size.y * grid.local_size.z;
        double *values, *bvalues;
        double *envPC;

        envPC = (double *) calloc(nvalues, sizeof(double));
        values = calloc(nvalues, sizeof(double));
        bvalues = calloc(nvalues, sizeof(double));


        /* 5. SETUP STRUCT VECTORS FOR B AND X */

        /* create an empty vector object */
        HYPRE_SStructVectorCreate(simstate.MPI_CART_COMM, envGrid, &b);
        HYPRE_SStructVectorCreate(simstate.MPI_CART_COMM, envGrid, &x);

        /* as with the matrix, set the appropriate object type for the vectors */
        HYPRE_SStructVectorSetObjectType(b, envObjectType);
        HYPRE_SStructVectorSetObjectType(x, envObjectType);

        /* indicate that the vector coefficients are ready to be set */
        HYPRE_SStructVectorInitialize(b);
        HYPRE_SStructVectorInitialize(x);

        for(var=0; var<nnutrients; var++) {

                envCellPC(simstate,envPC,var,grid);

                /* set the values */
                mi = 0;
                for (k = envLower[2]; k <= envUpper[2]; k++)
                        for (j = envLower[1]; j <= envUpper[1]; j++)
                                for (i = envLower[0]; i <= envUpper[0]; i++) {
                                        values[mi] = nutrient[var].initial_condition_mean;
                                        mi++;
                                }

                HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper, var,
                                                values);

                mi = 0;
                for (k = envLower[2]; k <= envUpper[2]; k++)
                        for (j = envLower[1]; j <= envUpper[1]; j++)
                                for (i = envLower[0]; i <= envUpper[0]; i++) {
                                        values[mi] = nutrient[var].initial_condition_mean;
                                        mi++;
                                }

                HYPRE_SStructVectorSetBoxValues(x, 0, envLower, envUpper, var,
                                                values);

                /* incorporate boundary conditions; Dirichlet on 6 faces */

                for (i = 0; i < nvalues; i++)
                        values[i] = envZ[var];
                for (i = 0; i < nvalues; i++)
                        bvalues[i] = envZ[var] * nutrient[var].boundary_condition;

                if (simstate.MPIcoords[simstate.MPIrank][0] == 0) {
                        nvalues = nentries * grid.local_size.y * grid.local_size.z;
                        envSetBoundary(0, -1);
                        stencil_indices[0] = 1;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }
                if (simstate.MPIcoords[simstate.MPIrank][0] == simstate.MPIdim[0] - 1) {
                        nvalues = nentries * grid.local_size.y * grid.local_size.z;
                        envSetBoundary(0, 1);
                        stencil_indices[0] = 2;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }
                if (simstate.MPIcoords[simstate.MPIrank][1] == 0) {
                        nvalues = nentries * grid.local_size.x * grid.local_size.z;
                        envSetBoundary(1, -1);
                        stencil_indices[0] = 3;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }
                if (simstate.MPIcoords[simstate.MPIrank][1] == simstate.MPIdim[1] - 1) {
                        nvalues = nentries * grid.local_size.x * grid.local_size.z;
                        envSetBoundary(1, 1);
                        stencil_indices[0] = 4;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }
                if (simstate.MPIcoords[simstate.MPIrank][2] == 0) {
                        nvalues = nentries * grid.local_size.x * grid.local_size.y;
                        envSetBoundary(2, -1);
                        stencil_indices[0] = 5;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }
                if (simstate.MPIcoords[simstate.MPIrank][2] == simstate.MPIdim[2] - 1) {
                        nvalues = nentries * grid.local_size.x * grid.local_size.y;
                        envSetBoundary(2, 1);
                        stencil_indices[0] = 6;
                        HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                                          nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                                          bvalues);
                }

                /* add production consumption function to the right side */
                HYPRE_SStructVectorAddToBoxValues(b, 0, envLower, envUpper, var,
                                                  envPC);

        }

        free(envPC);
        free(values);
        free(bvalues);
        /* stdout brought back */
}

/*!
 * This function initializes Hypre for solving a given envical field.
 */
void envInitSolver(struct state simstate)
{

        HYPRE_SStructMatrixAssemble(A);
        /* This is a collective call finalizing the vector assembly.
           The vector is now ``ready to be used'' */
        HYPRE_SStructVectorAssemble(b);
        HYPRE_SStructVectorAssemble(x);

        HYPRE_SStructMatrixGetObject(A, (void **) &parA);
        HYPRE_SStructVectorGetObject(b, (void **) &parb);
        HYPRE_SStructVectorGetObject(x, (void **) &parx);

        HYPRE_ParCSRPCGCreate(simstate.MPI_CART_COMM, &envSolver);
        HYPRE_ParCSRPCGSetTol(envSolver, 1.0e-12);
        HYPRE_ParCSRPCGSetPrintLevel(envSolver, 2);
        HYPRE_ParCSRPCGSetMaxIter(envSolver, 50);

        HYPRE_BoomerAMGCreate(&envPrecond);
        HYPRE_BoomerAMGSetMaxIter(envPrecond, 1);
        HYPRE_BoomerAMGSetTol(envPrecond, 0.0);
        HYPRE_BoomerAMGSetPrintLevel(envPrecond, 2);
        HYPRE_BoomerAMGSetCoarsenType(envPrecond, 6);
        HYPRE_BoomerAMGSetRelaxType(envPrecond, 6);
        HYPRE_BoomerAMGSetNumSweeps(envPrecond, 1);

        HYPRE_ParCSRPCGSetPrecond(envSolver, HYPRE_BoomerAMGSolve,
                                  HYPRE_BoomerAMGSetup, envPrecond);
        HYPRE_ParCSRPCGSetup(envSolver, parA, parb,
                             parx);

}

/*!
 * This is a driving function for solving next time step
 * of a given envical field.
 */
void envSolve(struct state simstate,int nnutrients,struct environment *nutrient,struct gridData grid)
{
        int i, j, k;
        int idx;
        int var;
        double *values;
        int stepIter = 0;
        long long nvalues = grid.local_size.x * grid.local_size.y * grid.local_size.z;
        double *envPC;
        if (simstate.MPIrank == 0 ) {
                printf("Solving field.");
                fflush(stdout);
        }

        values = (double *) calloc(nvalues, sizeof(double));
        envPC = (double *) calloc(nvalues, sizeof(double));

        while (stepIter < numberOfIters) {
                if (envIter > 0) {
                        /* update right hand side */
                        for(var=0; var<nnutrients; var++) {

                                envCellPC(simstate,envPC,var,grid);
                                HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper,
                                                                var, values);
                                HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper,
                                                                var, values);
                                for (i = 0; i < nvalues; i++)
                                        values[i] = envZ[var] * nutrient[var].boundary_condition;
                                if (simstate.MPIcoords[simstate.MPIrank][0] == 0) {
                                        envSetBoundary(0, -1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                if (simstate.MPIcoords[simstate.MPIrank][0] == simstate.MPIdim[0] - 1) {
                                        envSetBoundary(0, 1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                if (simstate.MPIcoords[simstate.MPIrank][1] == 0) {
                                        envSetBoundary(1, -1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                if (simstate.MPIcoords[simstate.MPIrank][1] == simstate.MPIdim[1] - 1) {
                                        envSetBoundary(1, 1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                if (simstate.MPIcoords[simstate.MPIrank][2] == 0) {
                                        envSetBoundary(2, -1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                if (simstate.MPIcoords[simstate.MPIrank][2] == simstate.MPIdim[2] - 1) {
                                        envSetBoundary(2, 1);
                                        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                                                          var, values);
                                }
                                HYPRE_SStructVectorAddToBoxValues(b, 0, envLower,
                                                                  envUpper, var, envPC);
                                HYPRE_SStructVectorAssemble(b);
                                HYPRE_SStructVectorAssemble(x);
                        }

                }

                HYPRE_ParCSRPCGSolve(envSolver, parA, parb,
                                     parx);

                for(var=0; var<nnutrients; var++) {
                        HYPRE_SStructVectorGather(x);
                        HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, var,
                                                        values);
                        /*         idx = 0;
                           for (k = 0; k < gridSize.z; k++)
                             for (j = 0; j < gridSize.y; j++)
                               for (i = 0; i < gridSize.x; i++, idx++) {
                                 //envField[nch][gridSize.y * gridSize.z * i + gridSize.z * j +
                                 //               k] = values[idx];
                                 printf("[%d,%d,%d] %.12f\n",i,j,k,values[idx]);
                               }
                         */
                }

                /* copy solution to field buffer */
/*   HYPRE_SStructVectorGather(x);
    HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, 0,
                                    values);
    idx = 0;
    for (k = 0; k < gridSize.z; k++)
      for (j = 0; j < gridSize.y; j++)
        for (i = 0; i < gridSize.x; i++, idx++) {
          envField[nch][gridSize.y * gridSize.z * i + gridSize.z * j +
                         k] = values[idx];
        }
 */
                envIter++;
                stepIter++;
        }

        free(values);
        free(envPC);
}

/*void allocateFields(int nnutrients,struct environment **nutrient,struct gridData grid) {

        int i,nf;

        dt=(double*)malloc(nnutrients*sizeof(double));
        fieldDt=(double*)malloc(nnutrients*sizeof(double));
        fieldAddr=(double**)malloc(nnutrients*sizeof(double*));
        envField=(double**)malloc(nnutrients*sizeof(double*));

        fieldConsumption=(double*)malloc(nnutrients*sizeof(double));
        fieldProduction=(double*)malloc(nnutrients*sizeof(double));
 * nutrient=(struct environment*)malloc(nnutrients*sizeof(struct environment));

        for(nf=0; nf<nnutrients; nf++) {
                fieldAddr[nf] =
                        (double *) calloc(grid.local_size.x * grid.local_size.y * grid.local_size.z,
                                          sizeof(double));
                envField[nf] = (double *) fieldAddr[nf];
                for (i = 0; i < grid.local_size.x * grid.local_size.y * grid.local_size.z; i++)
                        envField[nf][i] = 0.0;
        }

        printf("Fields allocation finished\n");
   }*/

/* ustawia testowe parametry dla fields */
/*void initFields(int nnutrients,struct environment *nutrient,struct gridData grid) {

        int i,j,k;
        int var;

        for(var=0; var<nnutrients; var++) {
                nutrient[var].diffusion_coefficient=1.82e-5;
                nutrient[var].lambda_delay=0.0;
                nutrient[var].boundary_condition=((double)var+1.0)*0.1575e-6;
                nutrient[var].initial_condition_mean=((double)var+1.0)*0.1575e-6;
                nutrient[var].initial_condition_variance=0.0;
                fieldConsumption[var]=8.3e-17;
                fieldProduction[var]=8.3e-1;
                fieldDt[var]=4000;
                for (k = 0; k < grid.local_size.z; k++)
                        for (j = 0; j < grid.local_size.y; j++)
                                for (i = 0; i < grid.local_size.x; i++) {
                                        envField[var][grid.local_size.y * grid.local_size.z * i + grid.local_size.z * j + k] =
                                                nutrient[var].initial_condition_mean;
                                }
        }

   }*/

/* to jest glowna funkcja "biblioteczna" */
void envsolve(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment) {

        envinitsystem(systeminfo,settings,grid,environment);
        //envinitboundary(simstate,nnutrients,nutrient,grid);
        //envinitsolver(simstate);
        //envSolve(simstate,nnutrients,nutrient,grid);
        return;
}

/* zmienne globalne z timothy-ego */
struct gridData grid;
struct environment *nutrient;
struct settings simsetup;
struct state simstate;
int nnutrients;
