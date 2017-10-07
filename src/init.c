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
#include <string.h>
#include <math.h>

#include "global.h"
#include "inline.h"

/*! \file init.c
 *  \brief contains initialization functions
 */

void getsysteminfo(system_t* system) {
  checkEndiannes();
  getLocalRankAndSize(system->rank, system->size, &(system->noderank), &(system->nodesize));
  system->memperproc = getMemoryPerProcess(system->nodesize);
  if (!POWER_OF_TWO(system->size))
    stopRun(101, NULL, __FILE__, __LINE__);
}

void initialsettings(settings_t* settings){
  settings->maxcells=0;
  settings->numberofsteps=0;
  settings->secondsperstep=0;
  settings->numberofcelltypes=0;
  settings->numberoffields=0;
  settings->dimension=3;
  settings->restart=0;
  strcpy(settings->rstfilename,"restart.bin");
  strcpy(settings->outdir,"results");
  settings->visoutstep=0;
  settings->statoutstep=0;
  settings->rstoutstep=0;
  settings->maxspeed=0;
}

void initialisation(int argc, char **argv, system_t system, settings_t* settings,celltype_t* celltype) {
  if (argc < 2 || argc >2) {
    if(system.rank==0) { printf("usage: timothy <parameter file>\n"); fflush(stdout); }
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  initialsettings(settings);
  readparamfile(argc,argv,system,settings);
  if (settings->numberofcelltypes<1) {
    if(system.rank==0) { printf("0 cell types. aborting.\n"); fflush(stdout); }
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  celltype=(celltype_t*)malloc((settings->numberofcelltypes)*sizeof(celltype_t));
//  readcellsfile(settings,celltype);
  //readenvfile();

//  if (strcmp(argv[1], "-h") == 0)
//    printHelp();

}

void printinfo(system_t system) {
  if (system.rank == 0) {
    printf("\ntimothy, tissue modelling framework\n");
    printf("http://timothy.icm.edu.pl\n");
    printf("version %s\n\n", VERSION);
    printf("number of MPI processes: %d\n",system.size);
    printf("processes per node: %d\n",system.nodesize);
    printf("threads per process: %d\n",system.nthreads);
    printf("system: ");
    if (endian)
      printf("%s, little endian\n", CPUARCH);
    else
      printf("%s, big endian\n", CPUARCH);
    printf("\n");
    fflush(stdout);
  }
}

/*!
 * This function sets default values for the simulation.
*/
void defaultValues()
{
  rst = 0;
  rstReset = 0;
  nhs = -1;
  tgs = 0;

  bvsim=0;
  bnsim=0;
  scsim=0;

  statOutStep = 1;
  rstOutStep = 1;
  vtkOutStep = 1;

  povout = 0;
  vtkout = 0;
  vnfout = 0;

  csizeInUnits = 10.0;		/* eukaryotic cell size is usually between 10-30 micrometers */
  cellVolume = (4.0 / 3.0) * M_PI * pow(csizeInUnits * 0.0001, 3.0);

}

void scInit()
{
  nscstages=5;
  sctprob=(double*)malloc((nscstages-1)*sizeof(double));
  nscinst=(int64_t*)calloc(nscstages,sizeof(int64_t));
  gnscinst=(int64_t*)calloc(nscstages,sizeof(int64_t));
  sctprob[0]=0.8;
  sctprob[1]=0.8;
  sctprob[2]=0.8;
  sctprob[3]=0.8;
  //sctprob[4]=0.8;
  localbc=0;
  globalbc=0;
}


/*!
 * This functions checks the consistence of key simulation parameters.
*/
void checkParameters()
{
  /* under construction */
  if (strcmp(cOutType, "VTK") && strcmp(cOutType, "POV")
      && strcmp(cOutType, "NON"))
    stopRun(116, "COUTTYPE", __FILE__, __LINE__);
  if (strcmp(fOutType, "VNF") && strcmp(fOutType, "NON"))
    stopRun(116, "FOUTTYPE", __FILE__, __LINE__);
}

/*!
 * This function intializes the simulation by setting all most important simulation parameters.
*/
void simulationInit(int argc, char **argv)
{
  int i;
  int periods[3];
  int reorder;

  /* set default values */
  defaultValues();

  /* stem cells parameters init */
  scInit();

  /* initialize random number generator */
  //randomStreamInit();

  /* read parameters file and restart file (if present) */
  //readParams(argc, argv);
  /* generate random cells if not a restart simulation */
  if (!rst) {
    simStart = 0;
    /* calculating number of cells per process */
    lnc = nc / MPIsize;
    if (MPIrank < nc % MPIsize)
      lnc++;
    /* allocating tables */
    cellsAllocate();
    if(bvsim) initVessel();
    if(bnsim) initBone();
    /* total number of cells - initialization */
//    for (i = 0; i < MPIsize; i++)
//      tlnc[i] = lnc;		//ZLE!!!!!
    /* cell cycle init */
    cellsCycleInit();
    /* random cell placement */
    //initVessel();
    cellsRandomInit();
    //initVessel();
    /* decomposition - initialization */
    decompositionInit(argc, argv, MPI_COMM_WORLD);
  }

  /* maximum distance cell can travel in 1 sec */
  maxSpeedInUnits = (maxSpeed * csize) / (24.0 * 60.0 * 60.0);
  /* at least one global fields iteration per step */
  gfDt = (gfDt > secondsPerStep ? secondsPerStep : gfDt);
  /* global fields iterations per step */
  gfIterPerStep = (int) (secondsPerStep / gfDt);

  if (sdim == 3)
    tnc = 8;
  if (sdim == 2)
    tnc = 4;

  printf("h=%f\n",h);
  printf("h3=%f %f\n",h3,h*h*h);

  /* density critical levels (very important parameters) */
  if (sdim == 3) {
    densityCriticalLevel1 = 6 * h3 * sph_kernel(1.6 * csize);	//1.8  //1.4 //1.75
    densityCriticalLevel2 = 6 * h3 * sph_kernel(1.1 * csize);	//1.1 //1.4
    printf("%f %f\n",densityCriticalLevel1,densityCriticalLevel2);
  }
  if (sdim == 2) {
    densityCriticalLevel1 = 4 * h2 * sph_kernel(1.4 * csize);	//1.4 //1.75
    densityCriticalLevel2 = 4 * h2 * sph_kernel(1.15 * csize);	//1.1 //1.4
  }

  /* checking the consistency */
  checkParameters();

  if (!strcmp(cOutType, "POV"))
    povout = 1;
  if (!strcmp(cOutType, "VTK"))
    vtkout = 1;
  if (!strcmp(fOutType, "VNF"))
    vnfout = 1;

  /* organizing processes in a Cartesian grid for global fields computations */
  MPI_Dims_create(MPIsize, sdim, MPIdim);
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, sdim, MPIdim, periods, reorder,
                  &MPI_CART_COMM);
  MPIcoords = (int **) malloc(MPIsize * sizeof(int *));
  for (i = 0; i < MPIsize; i++) {
    MPIcoords[i] = (int *) malloc(3 * sizeof(int));
    MPI_Cart_coords(MPI_CART_COMM, i, sdim, MPIcoords[i]);
  }
  /* compute grid size */
  computeGridSize();

  /* allocate grid data */
  allocateGrid();
  /* initialize global fields (might take some time) */
  fieldsInit();

  /* define colormaps for PovRay outputs */
  defineColormaps();

  /* define vessels - TBD */
  //initVessel();
}
