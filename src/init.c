/* **************************************************************************
 * Timothy - Tissue Modelling Framework
 * Copyright (C) 2014-2018 Maciej Cytowski
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
    terminate(system,"number of processes must be power of two", __FILE__, __LINE__);
  return;
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
  settings->gfdt=0;
  settings->gfh=0;
  settings->simulationstart=0;
  return;
}

void initialcelltype(int numberofcelltypes,int numberoffields,celltype_t* celltype){
  int i,j;
  for(i=0;i<numberofcelltypes;i++) {
    sprintf(celltype->name,"celltype%d",i);
  	celltype[i].g1=CELLTYPE_G1_DEFAULT;
  	celltype[i].s=CELLTYPE_S_DEFAULT;
    celltype[i].g2=CELLTYPE_G2_DEFAULT;
    celltype[i].m=CELLTYPE_M_DEFAULT;
    celltype[i].v=CELLTYPE_V_DEFAULT;
    celltype[i].rd=CELLTYPE_RD_DEFAULT;
    celltype[i].criticaldensity=CELLTYPE_CDENS_DEFAULT;
  }
  for(i=0;i<numberofcelltypes;i++) {
    for(j=0;j<numberoffields;j++) {
      celltype[i].production[j]=CELLTYPE_PROD_DEFAULT;
      celltype[i].consumption[j]=CELLTYPE_CONS_DEFAULT;
      celltype[i].criticallevel1[j]=CELLTYPE_CL1_DEFAULT;
      celltype[i].criticallevel2[j]=CELLTYPE_CL2_DEFAULT;
    }
  }
  return;
}

void initialfields(int numberoffields,environment_t* environment) {
  int i,j;
  for(i=0;i<numberoffields;i++) {
    sprintf(environment[i].name,"environment%d",i);
    environment[i].diffusioncoefficient=ENVIRONMENT_DC_DEFAULT;
    environment[i].boundarycondition=ENVIRONMENT_BC_DEFAULT;
    environment[i].initialconditionmean=ENVIRONMENT_ICMEAN_DEFAULT;
    environment[i].initialconditionvariance=ENVIRONMENT_ICVAR_DEFAULT;
    environment[i].lambdadelay=ENVIRONMENT_LAMBDA_DEFAULT;
  }
  return;
}

void initialisation(int argc, char **argv, system_t *system, settings_t* settings,celltype_t* celltype,environment_t* environment) {
  int i;
  int periods[3];
  int reorder;

  if (argc < 2 || argc >2) {
    if(system->rank==0) { printf("usage: timothy <parameter file>\n"); fflush(stdout); }
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  initialsettings(settings);
  readparamfile(argc,argv,*system,settings);
  if (settings->numberofcelltypes<1)
    terminate(system,"no cell types specified", __FILE__, __LINE__);

  if(!(celltype=(celltype_t*)malloc((settings->numberofcelltypes)*sizeof(celltype_t))))
    terminate(system,"cannot allocate celltype", __FILE__, __LINE__);
  if(!(environment=(environment_t*)malloc((settings->numberoffields)*sizeof(environment_t))))
    terminate(system,"cannot allocate environment", __FILE__, __LINE__);

  for(i=0;i<settings->numberofcelltypes;i++) {
    int size=settings->numberoffields*sizeof(float);
    celltype[i].production=(float*)malloc(size);
    celltype[i].consumption=(float*)malloc(size);
    celltype[i].criticallevel1=(float*)malloc(size);
    celltype[i].criticallevel2=(float*)malloc(size);
  }
  initialcelltype(settings->numberofcelltypes,settings->numberoffields,celltype);
  readcellsfile(*system,settings,celltype);
  initialfields(settings->numberoffields,environment);
  readenvfile(*system,settings,environment);

  /* organizing processes in a Cartesian grid for global fields computations */
  MPI_Dims_create(system->size, settings->dimension, system->dim);
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, settings->dimension, system->dim, periods, reorder,
                  &MPI_CART_COMM);

  if(!(system->coords = (int **) malloc(system->size * sizeof(int *))))
    terminate(system,"cannot allocate system->coords", __FILE__, __LINE__);
  for (i = 0; i < system->size; i++) {
    if(!(system->coords[i] = (int *) malloc(3 * sizeof(int))))
      terminate(system,"cannot allocate system->coords[i]", __FILE__, __LINE__);
    MPI_Cart_coords(MPI_CART_COMM, i, settings->dimension, system->coords[i]);
  }

  /* allocating tables */
  //cellsallocate();
  /* cell cycle init */
  //cellsCycleInit();
  /* random cell placement */
  //cellsRandomInit();
  /* decomposition - initialization */
  //decompositionInit(argc, argv, MPI_COMM_WORLD);

  if(settings->restart) {
    // read restart file
  }
//  if (strcmp(argv[1], "-h") == 0)
//    printHelp();

}

void initcount(cellcount_t *cellcount) {
  cellcount->n=0;
  cellcount->g0phase=0;
  cellcount->g1phase=0;
  cellcount->sphase=0;
  cellcount->g2phase=0;
  cellcount->mphase=0;
  cellcount->necroticphase=0;
  return;
}

void allocatecells(system_t system,settings_t settings,cellsinfo_t *cellsinfo) {
  int i;
  int maxlocalcells;
  maxlocalcells=settings.maxcells / system.size;
  if (system.rank < maxlocalcells % system.size)
    maxlocalcells++;
  if(!(cellsinfo->cells=(celldata_t*)malloc(maxlocalcells*sizeof(celldata_t))))
    terminate(system,"cannot allocate cellsinfo->cells", __FILE__, __LINE__);
  if(!(cellsinfo->forces=(double3dv_t*)malloc(maxlocalcells*sizeof(celldata_t))))
    terminate(system,"cannot allocate cellsinfo->forces", __FILE__, __LINE__);
  if(!(cellsinfo->cellsperproc=(uint64_t*)malloc(sizeof(uint64_t)*system.size)))
    terminate(system,"cannot allocate cellsinfo->cellsperproc", __FILE__, __LINE__);
  if(!(cellsinfo->typecount=(cellcount_t*)malloc(sizeof(cellcount_t)*settings.numberoffields)))
    terminate(system,"cannot allocate cellsinfo->typecount", __FILE__, __LINE__);
  initcount(&(cellsinfo->localcount));
  initcount(&(cellsinfo->globalcount));
  for(i=0;i<settings.numberoffields;i++)
    initcount(&(cellsinfo->typecount[i]));
  for(i=0;i<system.size;i++)
    cellsinfo->cellsperproc[i]=0;
  return;
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

/*!
 * This function intializes the simulation by setting all most important simulation parameters.
*/
void simulationInit(int argc, char **argv)
{
  int i;
  int periods[3];
  int reorder;

  /* set default values */
  //defaultValues();

  /* stem cells parameters init */
  //scInit();

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
  //checkParameters();

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
