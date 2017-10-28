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
#include <mpi.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#define _GNU_SOURCE
#include <fcntl.h>
#include <sys/stat.h>

#include "global.h"
#include "io.h"
#include "fields.h"

/*! \file io.c
*  \brief contains I/O functions
*/

#define INPUT_FILE_LINE_LENGTH 256

void readRstFile(int argc, char **argv);


/*!
* This function initializes all parameters.
*/

//void initParams(int argc, char **argv)
/*{

int nr;

nr = 0;

strcpy(params[nr], "SIZEX");
strcpy(desc[nr], "box size (X coordinate)");
req[nr] = 1;
addr[nr] = &nx;
type[nr] = INT;
nr++;

strcpy(params[nr], "SIZEY");
strcpy(desc[nr], "box size (Y coordinate)");
req[nr] = 1;
addr[nr] = &ny;
type[nr] = INT;
nr++;

strcpy(params[nr], "SIZEZ");
strcpy(desc[nr], "box size (Z coordinate)");
req[nr] = 1;
addr[nr] = &nz;
type[nr] = INT;
nr++;

strcpy(params[nr], "RSTFILE");
strcpy(desc[nr], "restart file name");
req[nr] = 0;
addr[nr] = rstFileName;
type[nr] = STRING;
nr++;

strcpy(params[nr], "NC");
strcpy(desc[nr], "number of cells");
req[nr] = 1;
addr[nr] = &nc;
type[nr] = LONG;
nr++;

strcpy(params[nr], "RNG");
strcpy(desc[nr], "random number generator type");
req[nr] = 1;
addr[nr] = rng;
type[nr] = STRING;
nr++;

strcpy(params[nr], "H");
strcpy(desc[nr], "SPH kernel function diameter");
req[nr] = 0;
addr[nr] = &h;
type[nr] = REAL;
nr++;

strcpy(params[nr], "NSTEPS");
strcpy(desc[nr], "number of simulation steps");
req[nr] = 1;
addr[nr] = &nsteps;
type[nr] = INT;
nr++;

strcpy(params[nr], "OUTDIR");
strcpy(desc[nr], "output directory");
req[nr] = 0;
addr[nr] = &outdir;
type[nr] = STRING;
nr++;

strcpy(params[nr], "G1");
strcpy(desc[nr],
"mean duration of G1 phase (in hours) - healthy tissue");
req[nr] = 1;
addr[nr] = &g1;
type[nr] = REAL;
nr++;

strcpy(params[nr], "S");
strcpy(desc[nr], "mean duration of S phase (in hours) - healthy tissue");
req[nr] = 1;
addr[nr] = &s;
type[nr] = REAL;
nr++;

strcpy(params[nr], "G2");
strcpy(desc[nr],
"mean duration of G2 phase (in hours) - healthy tissue");
req[nr] = 1;
addr[nr] = &g2;
type[nr] = REAL;
nr++;

strcpy(params[nr], "M");
strcpy(desc[nr], "mean duration of M phase (in hours) - healthy tissue");
req[nr] = 1;
addr[nr] = &m;
type[nr] = REAL;
nr++;

strcpy(params[nr], "SECPERSTEP");
strcpy(desc[nr], "time step");
req[nr] = 1;
addr[nr] = &secondsPerStep;
type[nr] = REAL;
nr++;

strcpy(params[nr], "DIM");
strcpy(desc[nr], "dimensionality of the system (2D/3D)");
req[nr] = 1;
addr[nr] = &sdim;
type[nr] = INT;
nr++;

strcpy(params[nr], "MITRAND");
strcpy(desc[nr], "mitosis random direction (1/0)");
req[nr] = 1;
addr[nr] = &mitrand;
type[nr] = INT;
nr++;

strcpy(params[nr], "V");
strcpy(desc[nr], "variability of duration of cell cycles, 0<V<1");
req[nr] = 1;
addr[nr] = &v;
type[nr] = REAL;
nr++;

strcpy(params[nr], "RD");
strcpy(desc[nr],
"radnom death - probability (for each cell) of being marked for dying. 0.0<=RD<1.0");
req[nr] = 1;
addr[nr] = &rd;
type[nr] = REAL;
nr++;

strcpy(params[nr], "RSTRESET");
strcpy(desc[nr], "reset simulation parameters of restart file");
req[nr] = 1;
addr[nr] = &rstReset;
type[nr] = INT;
nr++;

strcpy(params[nr], "NHS");
strcpy(desc[nr],
"number of cells needed to activate random dying (for homeostasis)");
req[nr] = 0;
addr[nr] = &nhs;
type[nr] = LONG;
nr++;

strcpy(params[nr], "TGS");
strcpy(desc[nr], "switches on tumor growth simulation");
req[nr] = 1;
addr[nr] = &tgs;
type[nr] = INT;
nr++;

strcpy(params[nr], "STATOUTSTEP");
strcpy(desc[nr], "every how many steps statistics are printed");
req[nr] = 1;
addr[nr] = &statOutStep;
type[nr] = INT;
nr++;

strcpy(params[nr], "RSTOUTSTEP");
strcpy(desc[nr], "every how many steps restart file is printed");
req[nr] = 1;
addr[nr] = &rstOutStep;
type[nr] = INT;
nr++;

strcpy(params[nr], "VISOUTSTEP");
strcpy(desc[nr], "every how many steps VTK file is printed");
req[nr] = 1;
addr[nr] = &vtkOutStep;
type[nr] = INT;
nr++;

strcpy(params[nr], "CG1");
strcpy(desc[nr], "mean duration of G1 phase (in hours) - cancer cells");
req[nr] = 1;
addr[nr] = &cg1;
type[nr] = REAL;
nr++;

strcpy(params[nr], "CS");
strcpy(desc[nr], "mean duration of S phase (in hours) - cancer cells");
req[nr] = 1;
addr[nr] = &cs;
type[nr] = REAL;
nr++;

strcpy(params[nr], "CG2");
strcpy(desc[nr], "mean duration of G2 phase (in hours) - cancer cells");
req[nr] = 1;
addr[nr] = &cg2;
type[nr] = REAL;
nr++;

strcpy(params[nr], "CM");
strcpy(desc[nr], "mean duration of M phase (in hours) - cancer cells");
req[nr] = 1;
addr[nr] = &cm;
type[nr] = REAL;
nr++;

strcpy(params[nr], "MAXSPEED");
strcpy(desc[nr],
"maximal displacement of cells in one time step (0.0<MAXMOVE<1.0)");
req[nr] = 1;
addr[nr] = &maxSpeed;
type[nr] = REAL;
nr++;

strcpy(params[nr], "GFLOGDIR");
strcpy(desc[nr], "log directory");
req[nr] = 0;
addr[nr] = &logdir;
type[nr] = STRING;
nr++;

strcpy(params[nr], "GFDT");
strcpy(desc[nr],
"the length of time step for solving global fields (unit: seconds)");
req[nr] = 1;
addr[nr] = &gfDt;
type[nr] = REAL;
nr++;

strcpy(params[nr], "GFH");
strcpy(desc[nr],
"grid resolution for solving global fields (unit: cell size)");
req[nr] = 1;
addr[nr] = &gfH;
type[nr] = REAL;
nr++;

strcpy(params[nr], "GFIELDS");
strcpy(desc[nr],
"do we use global fields in the simulation? (0 - no, 1 - yes)");
req[nr] = 1;
addr[nr] = &gfields;
type[nr] = INT;
nr++;

strcpy(params[nr], "OXYGEN");
strcpy(desc[nr],
"do we use oxygen field in the simulation? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &oxygen;
type[nr] = INT;
nr++;

strcpy(params[nr], "OXYGENDC");
strcpy(desc[nr], "oxygen field diffusion coefficient");
req[nr] = 0;
addr[nr] = &fieldDiffCoef[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENBC");
strcpy(desc[nr],
"oxygen field boundary condition (Dirichlet), mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldBC[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENICMEAN");
strcpy(desc[nr], "oxygen field initial condition mean, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICMean[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENICVAR");
strcpy(desc[nr], "oxygen field initial condition variance, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICVar[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENCONS");
strcpy(desc[nr], "oxygen field consumption, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldConsumption[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENPROD");
strcpy(desc[nr], "oxygen field production, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldProduction[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENLAMBDA");
strcpy(desc[nr], "oxygen field lambda");
req[nr] = 0;
addr[nr] = &fieldLambda[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENCL1");
strcpy(desc[nr], "oxygen field critical level 1, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel1[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "OXYGENCL2");
strcpy(desc[nr], "oxygen field critical level 2, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel2[OXYG];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSE");
strcpy(desc[nr],
"do we use glucose field in the simulation? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &glucose;
type[nr] = INT;
nr++;

strcpy(params[nr], "GLUCOSEDC");
strcpy(desc[nr], "glucose field diffusion coefficient");
req[nr] = 0;
addr[nr] = &fieldDiffCoef[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSEBC");
strcpy(desc[nr],
"glucose field boundary condition (Dirichlet), mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldBC[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSEICMEAN");
strcpy(desc[nr], "glucose field initial condition mean, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICMean[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSEICVAR");
strcpy(desc[nr], "glucose field initial condition variance, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICVar[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSECONS");
strcpy(desc[nr], "glucose field consumption, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldConsumption[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSEPROD");
strcpy(desc[nr], "glucose field production, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldProduction[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSELAMBDA");
strcpy(desc[nr], "glucose field lambda");
req[nr] = 0;
addr[nr] = &fieldLambda[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSECL1");
strcpy(desc[nr], "glucose field critical level 1, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel1[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "GLUCOSECL2");
strcpy(desc[nr], "glucose field critical level 2, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel2[GLUC];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENION");
strcpy(desc[nr],
"do we use hydrogen ion field in the simulation? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &hydrogenIon;
type[nr] = INT;
nr++;

strcpy(params[nr], "HYDROGENIONDC");
strcpy(desc[nr], "hydrogen ion field diffusion coefficient");
req[nr] = 0;
addr[nr] = &fieldDiffCoef[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONBC");
strcpy(desc[nr],
"hydrogen ion field boundary condition (Dirichlet), mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldBC[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONICMEAN");
strcpy(desc[nr], "hydrogen ion field initial condition mean, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICMean[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONICVAR");
strcpy(desc[nr],
"hydrogen ion field initial condition variance, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldICVar[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONCONS");
strcpy(desc[nr], "hydrogen ion field consumption, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldConsumption[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONPROD");
strcpy(desc[nr], "hydrogen ion field production, mol/(cell s)");
req[nr] = 0;
addr[nr] = &fieldProduction[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONLAMBDA");
strcpy(desc[nr], "hydrogen ion field lambda");
req[nr] = 0;
addr[nr] = &fieldLambda[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONCL1");
strcpy(desc[nr], "hydrogen ion field critical level 1, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel1[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "HYDROGENIONCL2");
strcpy(desc[nr], "hydrogen ion field critical level 2, mol/cm^3");
req[nr] = 0;
addr[nr] = &fieldCriticalLevel2[HYDR];
type[nr] = DOUBLE;
nr++;

strcpy(params[nr], "TEMPERATURE");
strcpy(desc[nr],
"do we use temperature field in the simulation? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &temperature;
type[nr] = INT;
nr++;

strcpy(params[nr], "MAXCELLS");
strcpy(desc[nr], "maximum number of cells in the simulation");
req[nr] = 1;
addr[nr] = &maxCells;
type[nr] = LONG;
nr++;

strcpy(params[nr], "COUTTYPE");
strcpy(desc[nr], "type of cellular data output (VTK or POV)");
req[nr] = 1;
addr[nr] = cOutType;
type[nr] = STRING;
nr++;

strcpy(params[nr], "FOUTTYPE");
strcpy(desc[nr], "type of fields data output (VNF)");
req[nr] = 1;
addr[nr] = fOutType;
type[nr] = STRING;
nr++;

strcpy(params[nr], "SCSIM");
strcpy(desc[nr],
"is this stem cell simulation? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &scsim;
type[nr] = INT;
nr++;

strcpy(params[nr], "BVSIM");
strcpy(desc[nr],
"do we simulate blood vessels? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &bvsim;
type[nr] = INT;
nr++;

strcpy(params[nr], "BNSIM");
strcpy(desc[nr],
"do we simulate bones? (0 - no, 1 - yes)");
req[nr] = 0;
addr[nr] = &bnsim;
type[nr] = INT;
nr++;

}*/

/*!
* This function reads the parameter file.
*/

/*
  strcpy(params[nr], "GFLOGDIR");
  strcpy(desc[nr], "log directory");
  req[nr] = 0;
  addr[nr] = &logdir;
  type[nr] = STRING;
  nr++;

  strcpy(params[nr], "GFDT");
  strcpy(desc[nr],
  "the length of time step for solving global fields (unit: seconds)");
  req[nr] = 1;
  addr[nr] = &gfDt;
  type[nr] = REAL;
  nr++;

  strcpy(params[nr], "GFH");
  strcpy(desc[nr],
  "grid resolution for solving global fields (unit: cell size)");
  req[nr] = 1;
  addr[nr] = &gfH;
  type[nr] = REAL;
  nr++;

  strcpy(params[nr], "GFIELDS");
  strcpy(desc[nr],
  "do we use global fields in the simulation? (0 - no, 1 - yes)");
  req[nr] = 1;
  addr[nr] = &gfields;
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "OXYGEN");
  strcpy(desc[nr],
  "do we use oxygen field in the simulation? (0 - no, 1 - yes)");
  req[nr] = 0;
  addr[nr] = &oxygen;
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "OXYGENDC");
  strcpy(desc[nr], "oxygen field diffusion coefficient");
  req[nr] = 0;
  addr[nr] = &fieldDiffCoef[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENBC");
  strcpy(desc[nr],
  "oxygen field boundary condition (Dirichlet), mol/cm^3");
  req[nr] = 0;
  addr[nr] = &fieldBC[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENICMEAN");
  strcpy(desc[nr], "oxygen field initial condition mean, mol/cm^3");
  req[nr] = 0;
  addr[nr] = &fieldICMean[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENICVAR");
  strcpy(desc[nr], "oxygen field initial condition variance, mol/cm^3");
  req[nr] = 0;
  addr[nr] = &fieldICVar[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENCONS");
  strcpy(desc[nr], "oxygen field consumption, mol/(cell s)");
  req[nr] = 0;
  addr[nr] = &fieldConsumption[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENPROD");
  strcpy(desc[nr], "oxygen field production, mol/(cell s)");
  req[nr] = 0;
  addr[nr] = &fieldProduction[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENLAMBDA");
  strcpy(desc[nr], "oxygen field lambda");
  req[nr] = 0;
  addr[nr] = &fieldLambda[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENCL1");
  strcpy(desc[nr], "oxygen field critical level 1, mol/cm^3");
  req[nr] = 0;
  addr[nr] = &fieldCriticalLevel1[OXYG];
  type[nr] = DOUBLE;
  nr++;

  strcpy(params[nr], "OXYGENCL2");
  strcpy(desc[nr], "oxygen field critical level 2, mol/cm^3");
  req[nr] = 0;
  addr[nr] = &fieldCriticalLevel2[OXYG];
  type[nr] = DOUBLE;
  nr++; */

void readenvfile(system_t system,settings_t* settings,environment_t* environment) {
  #define NENVPAR 6
  char envfile[FNLEN];
  char buf[400], buf1[100], buf2[100], buf3[100];
  FILE *fhandle;
  char params[NENVPAR][64];
  void *addr[NENVPAR];
  int req[NENVPAR];
  int set[NENVPAR];
  int type[NENVPAR];
  int i,j;
  int nr;
  int nameset=0;

  sprintf(envfile,"environment.inp");

  if (system.rank == 0)
  printf("reading environment file: %s\n", envfile);
  fflush(stdout);

  fhandle = fopen(envfile, "r");
  if (fhandle == NULL)
    terminate(system, "could not open environment file", __FILE__, __LINE__);

  for(i=0;i<settings->numberoffields;i++){

    nr=0;

    strcpy(params[nr], "ENVNAME");
    req[nr] = 0;
    addr[nr] = &(environment[i].name);
    type[nr] = STRING;
    nr++;

    strcpy(params[nr], "DC");
    req[nr] = 0;
    addr[nr] = &(environment[i].diffusioncoefficient);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "BC");
    req[nr] = 0;
    addr[nr] = &(environment[i].boundarycondition);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "ICMEAN");
    req[nr] = 0;
    addr[nr] = &(environment[i].initialconditionmean);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "ICVAR");
    req[nr] = 0;
    addr[nr] = &(environment[i].initialconditionvariance);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "LAMBDA");
    req[nr] = 0;
    addr[nr] = &(environment[i].lambdadelay);
    type[nr] = REAL;
    nr++;

    while(!feof(fhandle)) {
      char *ret;
      ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

      if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
      continue;

      if (feof(fhandle))
      break;

      if (buf1[0] == '#')
      continue;

      for (j = 1; j < NENVPAR; j++) {
        if (strcmp(buf1, params[j]) == 0) {
          switch (type[j]) {
            case REAL:
            *((float *) addr[j]) = atof(buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: environment %d read %s = %f\n", i, params[j], *((float *) addr[j]));
              fflush(stdout);
            }
            #endif
            break;
            case STRING:
            strcpy((char *) addr[j], buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: environment %d read %s = %s\n", i, params[j], buf2);
              fflush(stdout);
            }
            #endif
            break;
          }
        }
      }

      if (strcmp(buf1, params[0]) == 0) {
        if(nameset) {
          if(i+1==settings->numberoffields) break;
          else {
            strcpy((char *) addr[i+1], buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: read %s = %s\n", params[i], buf2);
              fflush(stdout);
            }
            #endif
            break;
          }
        } else {
          nameset++;
          strcpy((char *) addr[i], buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %s\n", params[i], buf2);
            fflush(stdout);
          }
          #endif
        }
      }
    }
  }
  fclose(fhandle);
  return;
}

void readcellsfile(system_t system, settings_t* settings, celltype_t* celltype) {
  #define NCELLPAR 13
  char cellsfile[FNLEN];
  char buf[400], buf1[100], buf2[100], buf3[100], buf4[100],buf5[100];
  FILE *fhandle;
  char params[NCELLPAR][64];
  void *addr[NCELLPAR];
  int req[NCELLPAR];
  int set[NCELLPAR];
  int type[NCELLPAR];
  int i,j,k;
  int nr;
  int nbytes,bytesread;
  int nameset=0;

  sprintf(cellsfile,"cells.inp");

  if (system.rank == 0)
  printf("reading cells file: %s\n", cellsfile);
  fflush(stdout);

  fhandle = fopen(cellsfile, "r");
  if (fhandle == NULL)
    terminate(system, "could not open cell type file", __FILE__, __LINE__);

  for(i=0;i<settings->numberofcelltypes;i++){

    nr=0;

    strcpy(params[nr], "CNAME");
    req[nr] = 0;
    addr[nr] = &(celltype[i].name);
    type[nr] = STRING;
    nr++;

    strcpy(params[nr], "G1");
    req[nr] = 0;
    addr[nr] = &(celltype[i].g1);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "S");
    req[nr] = 0;
    addr[nr] = &(celltype[i].s);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "G2");
    req[nr] = 0;
    addr[nr] = &(celltype[i].g2);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "M");
    req[nr] = 0;
    addr[nr] = &(celltype[i].m);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "V");
    req[nr] = 0;
    addr[nr] = &(celltype[i].v);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "RD");
    req[nr] = 1;
    addr[nr] = &(celltype[i].rd);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "INPUT");
    req[nr] = 0;
    addr[nr] = &(celltype[i].inputfile);
    type[nr] = STRING;
    nr++;

    strcpy(params[nr], "CDENS");
    req[nr] = 0;
    addr[nr] = &(celltype[i].criticaldensity);
    type[nr] = REAL;
    nr++;

    strcpy(params[nr], "ENVPROD");
    req[nr] = 0;
    type[nr] = REALVECTOR;
    addr[nr] = &(celltype[i].production[0]);
    nr++;

    strcpy(params[nr], "ENVCONS");
    req[nr] = 0;
    type[nr] = REALVECTOR;
    addr[nr] = &(celltype[i].consumption[0]);
    nr++;

    strcpy(params[nr], "ENVCL1");
    req[nr] = 0;
    type[nr] = REALVECTOR;
    addr[nr] = &(celltype[i].criticallevel1[0]);
    nr++;

    strcpy(params[nr], "ENVCL2");
    req[nr] = 0;
    type[nr] = REALVECTOR;
    addr[nr] = &(celltype[i].criticallevel2[0]);
    nr++;

    while(!feof(fhandle)) {
      char *ret;
      ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

      if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
      continue;

      if (feof(fhandle))
      break;

      if (buf1[0] == '#')
      continue;

      for (j = 1; j < NCELLPAR; j++) {
        if (strcmp(buf1, params[j]) == 0) {
          switch (type[j]) {

            case REAL:
            *((float *) addr[j]) = atof(buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: cell type %d read %s = %f\n", i, params[j], *((float *) addr[j]));
              fflush(stdout);
            }
            #endif
            break;

            case STRING:
            strcpy((char *) addr[j], buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: cell type %d read %s = %s\n", i, params[j], buf2);
              fflush(stdout);
            }
            #endif
            break;

            case REALVECTOR:
            bytesread=0;
            sscanf(buf, "%s%n",buf4,&nbytes);
            bytesread+=nbytes;
            for(k=0;k<settings->numberoffields;k++) {
              int ret;
              ret=sscanf(buf+bytesread,"%s%n",buf5,&nbytes);
              if(ret<=0) break;
              bytesread+=nbytes;
              ((float *) addr[j])[k] = atof(buf5);
              #ifdef DEBUG
              if (system.rank == 0) {
                printf("debug: cell type %d read %s[%d] = %f\n", i, params[j],k, ((float *) addr[j])[k]);
                fflush(stdout);
              }
              #endif
            }
            break;
          }
        }
      }

      if (strcmp(buf1, params[0]) == 0) {
        if(nameset) {
          if(i+1==settings->numberofcelltypes) break;
          else {
            strcpy((char *) addr[i+1], buf2);
            #ifdef DEBUG
            if (system.rank == 0) {
              printf("debug: read %s = %s\n", params[i], buf2);
              fflush(stdout);
            }
            #endif
            break;
          }
        } else {
          nameset++;
          strcpy((char *) addr[i], buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %s\n", params[i], buf2);
            fflush(stdout);
          }
          #endif
        }
      }
    }
  }
  fclose(fhandle);
  return;
}


void readparamfile(int argc, char **argv, system_t system, settings_t* settings)
{
  #define NPAR 15
  char paramfile[FNLEN];
  char buf[400], buf1[100], buf2[100], buf3[100];
  FILE *fhandle;
  char params[NPAR][64];
  void *addr[NPAR];
  int req[NPAR];
  int set[NPAR];
  int type[NPAR];
  int i;
  int nr;

  for(i=0;i<NPAR;i++)
    set[i]=0;

  nr = 0;

  strcpy(params[nr], "GFDT");
  req[nr] = 1;
  addr[nr] = &(settings->gfdt);
  type[nr] = REAL;
  nr++;

  strcpy(params[nr], "GFH");
  req[nr] = 1;
  addr[nr] = &(settings->gfh);
  type[nr] = REAL;
  nr++;

  strcpy(params[nr], "MAXCELLS");
  req[nr] = 1;
  addr[nr] = &(settings->maxcells);
  type[nr] = LONG;
  nr++;

  strcpy(params[nr], "NSTEPS");
  req[nr] = 1;
  addr[nr] = &(settings->numberofsteps);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "SECPERSTEP");
  req[nr] = 1;
  addr[nr] = &(settings->secondsperstep);
  type[nr] = REAL;
  nr++;

  strcpy(params[nr], "NCELLTYPES");
  req[nr] = 1;
  addr[nr] = &(settings->numberofcelltypes);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "NFIELDS");
  req[nr] = 1;
  addr[nr] = &(settings->numberoffields);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "DIMENSIONS");
  req[nr] = 1;
  addr[nr] = &(settings->dimension);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "RESTART");
  req[nr] = 1;
  addr[nr] = &(settings->restart);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "RSTFILE");
  req[nr] = 0;
  addr[nr] = &(settings->rstfilename);
  type[nr] = STRING;
  nr++;

  strcpy(params[nr], "OUTDIR");
  req[nr] = 0;
  addr[nr] = &(settings->outdir);
  type[nr] = STRING;
  nr++;

  strcpy(params[nr], "VISOUTSTEP");
  req[nr] = 1;
  addr[nr] = &(settings->visoutstep);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "STATOUTSTEP");
  req[nr] = 1;
  addr[nr] = &(settings->statoutstep);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "RSTOUTSTEP");
  req[nr] = 1;
  addr[nr] = &(settings->rstoutstep);
  type[nr] = INT;
  nr++;

  strcpy(params[nr], "MAXSPEED");
  req[nr] = 1;
  addr[nr] = &(settings->maxspeed);
  type[nr] = REAL;
  nr++;

  strcpy(settings->outdir, "results");

  if (strlen(argv[1]) >= FNLEN)
    terminate(system, "parameter file name too long", __FILE__, __LINE__);

  sprintf(paramfile, "%s", argv[1]);

  if (system.rank == 0)
  printf("reading parameters file: %s\n", paramfile);

  fflush(stdout);

  fhandle = fopen(paramfile, "r");
  if (fhandle == NULL)
    terminate(system, "could not open parameter file", __FILE__, __LINE__);

  /* look for parameters in the file */
  while (!feof(fhandle)) {

    char *ret;
    ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

    if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
    continue;

    if (feof(fhandle))
    break;

    if (buf1[0] == '#')
    continue;

    for (i = 0; i < NPAR; i++) {
      if (strcmp(buf1, params[i]) == 0
      && strcmp(params[i], "RSTFILE") == 0) {
        strcpy((char *) addr[i], buf2);
        set[i] = 1;
        rst = 1;
      }
      if (strcmp(buf1, params[i]) == 0
      && strcmp(params[i], "MAXCELLS") == 0) {
        *((int64_t *) addr[i]) = atol(buf2);
      }
    }
  }

  /* read restart file if given */
  if (settings->restart == 1)
  readRstFile(argc, argv);

  /* rewind the file */
  rewind(fhandle);

  /* convert some of the parameters to data */
  while (!feof(fhandle)) {

    char *ret;
    ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

    if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
    continue;

    if (feof(fhandle))
    break;

    if (buf1[0] == '#')
    continue;

    for (i = 0; i < NPAR; i++) {

      if (strcmp(buf1, params[i]) == 0) {
        switch (type[i]) {
          case REAL:
          *((float *) addr[i]) = atof(buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %f\n", params[i], *((float *) addr[i]));
            fflush(stdout);
          }
          #endif
          break;
          case DOUBLE:
          *((double *) addr[i]) = atof(buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %f\n", params[i], *((double *) addr[i]));
            fflush(stdout);
          }
          #endif
          break;
          case STRING:
          strcpy((char *) addr[i], buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %s\n", params[i], buf2);
            fflush(stdout);
          }
          #endif
          break;
          case INT:
          *((int *) addr[i]) = atoi(buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %d\n", params[i], *((int *) addr[i]));
            fflush(stdout);
          }
          #endif
          break;
          case LONG:
          *((int64_t *) addr[i]) = atol(buf2);
          #ifdef DEBUG
          if (system.rank == 0) {
            printf("debug: read %s = %" PRId64 "\n", params[i], *((int64_t *) addr[i]));
            fflush(stdout);
          }
          #endif
          break;
        }
        set[i] = 1;
        break;
      }
    }
  }

  for (i = 0; i < NPAR; i++)
  if (req[i] == 1 && set[i] == 0) {
    char errmsg[128];
    sprintf(errmsg,"missing parameter %s",params[i]);
    terminate(system,errmsg, __FILE__, __LINE__);
  }

  if (settings->maxspeed <= 0.0 || settings->maxspeed >= 4.0)
    terminate(system,"maxspeed out of range", __FILE__, __LINE__);

  if (system.rank == 0) {
    struct stat s;
    int err;
    printf("output directory: %s/\n", settings->outdir);

    err = stat(settings->outdir, &s);
    if (err == -1) {
      printf("creating output directory.\n");
      mkdir(settings->outdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
  }
  fclose(fhandle);
  return;
}

/*!
* This function defines output fields.
* Edit this part in order to introduce new output fields.
*/
void ioDefineOutputFields()
{

  nfOut = 0;

  if (lnc > 0) {

    strcpy(nameOut[nfOut], "density");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = REAL;
    addrOut[nfOut] = &cells[0].density;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].density - (int64_t) & cells[0].density;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "size");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = REAL;
    addrOut[nfOut] = &cells[0].size;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].size - (int64_t) & cells[0].size;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "rank");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &MPIrank;
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "phase");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].phase;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].phase - (int64_t) & cells[0].phase;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "tumor");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].tumor;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].tumor - (int64_t) & cells[0].tumor;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "halo");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].halo;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].halo - (int64_t) & cells[0].halo;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "velocity");
    dimOut[nfOut] = VECTOR;
    typeOut[nfOut] = REAL;
    addrOut[nfOut] = &velocity[0];
    if (lnc > 1)
    jumpOut[nfOut] = (int64_t) & velocity[1] - (int64_t) & velocity[0];
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "age");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].age;
    if (lnc > 1)
    jumpOut[nfOut] = (int64_t) & cells[1].age - (int64_t) & cells[0].age;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "scalarField");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = REAL;
    addrOut[nfOut] = &cells[0].scalarField;
    if (lnc > 1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].scalarField -
    (int64_t) & cells[0].scalarField;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "ctype");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].ctype;
    if(lnc>1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].ctype -
    (int64_t) & cells[0].ctype;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

    strcpy(nameOut[nfOut], "scstage");
    dimOut[nfOut] = SCALAR;
    typeOut[nfOut] = INT;
    addrOut[nfOut] = &cells[0].scstage;
    if(lnc>1)
    jumpOut[nfOut] =
    (int64_t) & cells[1].scstage -
    (int64_t) & cells[0].scstage;
    else
    jumpOut[nfOut] = 0;

    nfOut++;

  }
}

/*!
* This function writes a VTK file with all cells for a given step.
*/
void ioWriteStepVTK(int step)
{

  int i, j;
  MPI_File fh;
  float *floatVectorField;
  float *floatScalarField;
  int *integerScalarField;
  int64_t nprev = 0;
  char fstname[256];
  char header[256];
  MPI_Offset offset, goffset;

  ioDefineOutputFields();

  floatVectorField = (float *) calloc(lnc * 3, sizeof(float));
  floatScalarField = (float *) calloc(lnc, sizeof(float));
  integerScalarField = (int *) calloc(lnc, sizeof(int));

  sprintf(fstname, "%s/step%08d.vtk", outdir, step);

  goffset = 0;
  MPI_File_open(MPI_COMM_WORLD, fstname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
    MPI_INFO_NULL, &fh);
    /* truncate the file */
    MPI_File_set_size(fh, 0);

    sprintf(header,
      "# vtk DataFile Version 2.0\nTimothy output\nBINARY\nDATASET UNSTRUCTURED_GRID\n");

      /* gather number of cells from each process */
      MPI_Allgather(&lnc, 1, MPI_LONG_LONG, tlnc, 1, MPI_LONG_LONG,
        MPI_COMM_WORLD);
        for (i = 0; i < MPIrank; i++)
        nprev += tlnc[i];


        /* write the VTK header */
        if (MPIrank == 0)
        MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
        MPI_STATUS_IGNORE);
        goffset += strlen(header);
        MPI_File_seek(fh, goffset, MPI_SEEK_SET);

        /* adding positions */
        memset(header, 0, 256);
        sprintf(header, "\nPOINTS %" PRId64 " float\n", nc);
        if (MPIrank == 0)
        MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
        MPI_STATUS_IGNORE);
        goffset += strlen(header);
        MPI_File_seek(fh, goffset, MPI_SEEK_SET);

        offset = nprev * sizeof(float) * 3;
        MPI_File_seek(fh, offset, MPI_SEEK_CUR);
        for (j = 0; j < lnc; j++) {
          floatVectorField[3 * j] = (float) (cells[j].x);
          floatVectorField[3 * j + 1] = (float) (cells[j].y);
          floatVectorField[3 * j + 2] = (float) (cells[j].z);
        }
        if (endian)
        swap_Nbyte((char *) floatVectorField, lnc * 3, sizeof(float));
        MPI_File_write(fh, floatVectorField, 3 * lnc, MPI_FLOAT,
          MPI_STATUS_IGNORE);
          goffset += nc * sizeof(float) * 3;
          MPI_File_seek(fh, goffset, MPI_SEEK_SET);

          /* adding cell types */
          memset(header, 0, 256);
          sprintf(header, "\nCELL_TYPES %" PRId64 "\n", nc);
          if (MPIrank == 0)
          MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
          MPI_STATUS_IGNORE);
          goffset += strlen(header);
          MPI_File_seek(fh, goffset, MPI_SEEK_SET);
          offset = nprev * sizeof(int);
          MPI_File_seek(fh, offset, MPI_SEEK_CUR);
          for (j = 0; j < lnc; j++)
          integerScalarField[j] = 1;
          if (endian)
          swap_Nbyte((char *) integerScalarField, lnc, sizeof(int));
          MPI_File_write(fh, integerScalarField, lnc, MPI_INT, MPI_STATUS_IGNORE);
          goffset += nc * sizeof(int);
          MPI_File_seek(fh, goffset, MPI_SEEK_SET);

          /* point data */
          memset(header, 0, 256);
          sprintf(header, "\nPOINT_DATA %" PRId64, nc);
          if (MPIrank == 0)
          MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
          MPI_STATUS_IGNORE);
          goffset += strlen(header);
          MPI_File_seek(fh, goffset, MPI_SEEK_SET);

          /* adding fields */
          for (i = 0; i < nfOut; i++) {
            memset(header, 0, 256);
            if (dimOut[i] == SCALAR) {
              switch (typeOut[i]) {
                case REAL:
                sprintf(header, "\nSCALARS %s float 1\nLOOKUP_TABLE default\n",
                nameOut[i]);
                if (MPIrank == 0)
                MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
                MPI_STATUS_IGNORE);
                goffset += strlen(header);
                MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                for (j = 0; j < lnc; j++)
                floatScalarField[j] =
                (float) (*((double *) (addrOut[i] + j * jumpOut[i])));
                offset = nprev * sizeof(float);
                MPI_File_seek(fh, offset, MPI_SEEK_CUR);
                if (endian)
                swap_Nbyte((char *) floatScalarField, lnc, sizeof(float));
                MPI_File_write(fh, floatScalarField, lnc, MPI_FLOAT,
                  MPI_STATUS_IGNORE);
                  goffset += nc * sizeof(float);
                  MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                  break;
                  case INT:
                  sprintf(header, "\nSCALARS %s integer 1\nLOOKUP_TABLE default\n",
                  nameOut[i]);
                  if (MPIrank == 0)
                  MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
                  MPI_STATUS_IGNORE);
                  goffset += strlen(header);
                  MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                  for (j = 0; j < lnc; j++)
                  integerScalarField[j] = *((int *) (addrOut[i] + j * jumpOut[i]));
                  offset = nprev * sizeof(int);
                  MPI_File_seek(fh, offset, MPI_SEEK_CUR);
                  if (endian)
                  swap_Nbyte((char *) integerScalarField, lnc, sizeof(int));
                  MPI_File_write(fh, integerScalarField, lnc, MPI_INT,
                    MPI_STATUS_IGNORE);
                    goffset += nc * sizeof(int);
                    MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                    break;
                  }
                }
                if (dimOut[i] == VECTOR) {
                  switch (typeOut[i]) {
                    case REAL:
                    sprintf(header, "\nVECTORS %s float\n", nameOut[i]);
                    if (MPIrank == 0)
                    MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
                    MPI_STATUS_IGNORE);
                    goffset += strlen(header);
                    MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                    for (j = 0; j < lnc; j++) {
                      floatVectorField[3 * j] =
                      (float) (*((double *) (addrOut[i] + j * jumpOut[i])));
                      floatVectorField[3 * j + 1] =
                      (float) (*
                        ((double *) (addrOut[i] + 1 * sizeof(double) +
                        j * jumpOut[i])));
                        floatVectorField[3 * j + 2] =
                        (float) (*
                          ((double *) (addrOut[i] + 2 * sizeof(double) +
                          j * jumpOut[i])));
                        }
                        offset = nprev * sizeof(float) * 3;
                        MPI_File_seek(fh, offset, MPI_SEEK_CUR);
                        if (endian)
                        swap_Nbyte((char *) floatVectorField, lnc * 3, sizeof(float));
                        MPI_File_write(fh, floatVectorField, lnc * 3, MPI_FLOAT,
                          MPI_STATUS_IGNORE);
                          goffset += nc * 3 * sizeof(float);
                          MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                          break;
                          case INT:
                          /* note: INTs are converted to FLOATs */
                          printf(header, "\nVECTORS %s float\n", nameOut[i]);
                          if (MPIrank == 0)
                          MPI_File_write(fh, &header, strlen(header), MPI_BYTE,
                          MPI_STATUS_IGNORE);
                          goffset += strlen(header);
                          MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                          for (j = 0; j < lnc; j++) {
                            floatVectorField[3 * j] =
                            (*((float *) (addrOut[i] + j * jumpOut[i])));
                            floatVectorField[3 * j + 1] =
                            (*
                              ((float *) (addrOut[i] + 1 * sizeof(int) +
                              j * jumpOut[i])));
                              floatVectorField[3 * j + 2] =
                              (*
                                ((float *) (addrOut[i] + 2 * sizeof(int) +
                                j * jumpOut[i])));
                              }
                              offset = nprev * sizeof(float) * 3;
                              MPI_File_seek(fh, offset, MPI_SEEK_CUR);
                              if (endian)
                              swap_Nbyte((char *) floatVectorField, lnc * 3, sizeof(float));
                              MPI_File_write(fh, floatVectorField, lnc * 3, MPI_FLOAT,
                                MPI_STATUS_IGNORE);
                                goffset += nc * 3 * sizeof(float);
                                MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                                break;
                              }
                            }
                          }

                          free(floatVectorField);
                          free(floatScalarField);
                          free(integerScalarField);

                          MPI_File_close(&fh);
                        }

                        /*!
                        * This function prints the actual step number and
                        * total number of cells in the system.
                        */
                        void printStepNum()
                        {
                          if (MPIrank == 0) {
                            printf
                            ("\n-------------------------------------------------------------------------\n");
                            printf(" Step %8d,%15s%8.4f%20s%14" PRId64 " ", step, "Time step = ",
                            secondsPerStep, "Number of cells = ", nc);
                            fflush(stdout);
                            printf
                            ("\n-------------------------------------------------------------------------\n\n");
                            printf(" Time: %8.4f\n\n", simTime);
                          }
                        }

                        /*!
                        * This function defines global parametes for the restart file.
                        * Defined parameters are used during reading and writing of the restart file.
                        * Edit this function in order to introduce new restart file parameters.
                        */
                        void ioDefineRstGlobalParams()
                        {

                          nRst = 0;

                          /* endiannnes */
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &one;
                          nRst++;
                          /* system dimensions */
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &sdim;
                          nRst++;
                          /* box sizes */
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &nx;
                          nRst++;
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &ny;
                          nRst++;
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &nz;
                          nRst++;
                          /* output directory */
                          typeRst[nRst] = CHAR;
                          sizeRst[nRst] = 128;
                          addrRst[nRst] = outdir;
                          nRst++;
                          /* RNG type */
                          typeRst[nRst] = CHAR;
                          sizeRst[nRst] = 3;
                          addrRst[nRst] = rng;
                          nRst++;
                          /* number of cells */
                          typeRst[nRst] = INT64_T;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &nc;
                          nRst++;
                          /* "Simulation started" flag */
                          typeRst[nRst] = INT;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &simStart;
                          nRst++;
                          /* simulation time step */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &secondsPerStep;
                          nRst++;
                          /* simulation time */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &simTime;
                          nRst++;
                          /* fraction */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &dummy;
                          nRst++;
                          /* cell cycle phases duration - healthy tissue */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &g1;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &s;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &g2;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &m;
                          nRst++;
                          /* cell cycle variability */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &v;
                          nRst++;
                          /* random death probability */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &rd;
                          nRst++;
                          /* neighborhood h parameter */
                          typeRst[nRst] = DOUBLE;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &h;
                          nRst++;
                          /* cell size */
                          typeRst[nRst] = DOUBLE;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &csize;
                          nRst++;
                          /* cell cycle phases duration - cancer cells */
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &cg1;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &cs;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &cg2;
                          nRst++;
                          typeRst[nRst] = REAL;
                          sizeRst[nRst] = 1;
                          addrRst[nRst] = &cm;
                          nRst++;

                        }

                        /*!
                        * This function saves a restart file of the simulation.
                        */
                        void saveRstFile()
                        {

                          int i;
                          MPI_File fh;
                          int64_t nprev = 0;
                          char fstname[256];
                          MPI_Offset goffset;
                          MPI_Offset offset;

                          sprintf(fstname, "%s/step%08d.rst", outdir, step);

                          goffset = 0;
                          MPI_File_open(MPI_COMM_WORLD, fstname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh);
                            /* truncate the file */
                            MPI_File_set_size(fh, 0);

                            /* gather number of cells from each process */
                            MPI_Allgather(&lnc, 1, MPI_INT64_T, tlnc, 1, MPI_INT64_T,
                              MPI_COMM_WORLD);
                              for (i = 0; i < MPIrank; i++)
                              nprev += tlnc[i];

                              /* write out the simulation global variables and parameters (single process) */
                              if (MPIrank == 0) {

                                ioDefineRstGlobalParams();
                                offset = 0;

                                /* write out to rst file */
                                for (i = 0; i < nRst; i++) {
                                  switch (typeRst[i]) {
                                    case REAL:
                                    MPI_File_write(fh, addrRst[i], sizeRst[i], MPI_FLOAT,
                                      MPI_STATUS_IGNORE);
                                      offset += sizeRst[i] * sizeof(float);
                                      break;
                                      case INT:
                                      MPI_File_write(fh, addrRst[i], sizeRst[i], MPI_INT,
                                        MPI_STATUS_IGNORE);
                                        offset += sizeRst[i] * sizeof(int);
                                        break;
                                        case CHAR:
                                        MPI_File_write(fh, addrRst[i], sizeRst[i], MPI_CHAR,
                                          MPI_STATUS_IGNORE);
                                          offset += sizeRst[i] * sizeof(char);
                                          break;
                                          case DOUBLE:
                                          MPI_File_write(fh, addrRst[i], sizeRst[i], MPI_DOUBLE,
                                            MPI_STATUS_IGNORE);
                                            offset += sizeRst[i] * sizeof(double);
                                            break;
                                            case INT64_T:
                                            MPI_File_write(fh, addrRst[i], sizeRst[i], MPI_INT64_T,
                                              MPI_STATUS_IGNORE);
                                              offset += sizeRst[i] * sizeof(int64_t);
                                              break;
                                            }
                                          }
                                        }

                                        /* distribute information on global offset in the rst file */
                                        MPI_Bcast(&offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
                                        /* move the global file offset */
                                        goffset = offset + nprev * sizeof(struct cellData);

                                        MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                                        /* write out cells' data */
                                        MPI_File_write(fh, cells, lnc * sizeof(struct cellData), MPI_BYTE,
                                        MPI_STATUS_IGNORE);
                                        /* close file */
                                        MPI_File_close(&fh);
                                      }

                                      /*!
                                      * This function reads the restart file.
                                      */
                                      void readRstFile(int argc, char **argv)
                                      {

                                        int i;
                                        MPI_File fh;
                                        MPI_Offset goffset;
                                        int64_t nk, nr;
                                        int64_t nprev = 0;
                                        int swap;
                                        int lcancer = 0;

                                        if (MPIrank == 0)
                                        printf("Reading restart file: %s\n", rstFileName);

                                        goffset = 0;

                                        /* open restart file */
                                        MPI_File_open(MPI_COMM_WORLD, rstFileName, MPI_MODE_RDONLY,
                                          MPI_INFO_NULL, &fh);

                                          /* each process defines global parameters to be read from restart file */
                                          ioDefineRstGlobalParams();

                                          goffset = 0;

                                          for (i = 0; i < nRst; i++) {
                                            switch (typeRst[i]) {
                                              case REAL:
                                              MPI_File_read(fh, addrRst[i], sizeRst[i], MPI_FLOAT,
                                                MPI_STATUS_IGNORE);
                                                goffset += sizeRst[i] * sizeof(float);
                                                break;
                                                case INT:
                                                MPI_File_read(fh, addrRst[i], sizeRst[i], MPI_INT,
                                                  MPI_STATUS_IGNORE);
                                                  goffset += sizeRst[i] * sizeof(int);
                                                  break;
                                                  case CHAR:
                                                  MPI_File_read(fh, addrRst[i], sizeRst[i], MPI_CHAR,
                                                    MPI_STATUS_IGNORE);
                                                    goffset += sizeRst[i] * sizeof(char);
                                                    break;
                                                    case DOUBLE:
                                                    MPI_File_read(fh, addrRst[i], sizeRst[i], MPI_DOUBLE,
                                                      MPI_STATUS_IGNORE);
                                                      goffset += sizeRst[i] * sizeof(double);
                                                      break;
                                                      case INT64_T:
                                                      MPI_File_read(fh, addrRst[i], sizeRst[i], MPI_INT64_T,
                                                        MPI_STATUS_IGNORE);
                                                        goffset += sizeRst[i] * sizeof(int64_t);
                                                        break;
                                                      }
                                                    }

                                                    /* check endian */
                                                    char *p = (char *) &one;
                                                    swap = 0;
                                                    if (p[0] == 1) {
                                                      if (MPIrank == 0)
                                                      printf("Restart file in little endian format.\n");
                                                      if (endian == 0) {
                                                        if (MPIrank == 0)
                                                        printf("Conversion to big endian.\n");
                                                        swap = 1;
                                                      }
                                                    } else {
                                                      if (MPIrank == 0)
                                                      printf("Restart file in big endian format.\n");
                                                      if (endian == 1) {
                                                        if (MPIrank == 0)
                                                        printf("Conversion to little endian.\n");
                                                        swap = 1;
                                                      }
                                                    }

                                                    if (swap) {
                                                      for (i = 0; i < nRst; i++) {
                                                        switch (typeRst[i]) {
                                                          case REAL:
                                                          swap_Nbyte((char *) addrRst[i], sizeRst[i], sizeof(float));
                                                          break;
                                                          case INT:
                                                          swap_Nbyte((char *) addrRst[i], sizeRst[i], sizeof(int));
                                                          break;
                                                          case CHAR:
                                                          swap_Nbyte((char *) addrRst[i], sizeRst[i], sizeof(char));
                                                          break;
                                                          case DOUBLE:
                                                          swap_Nbyte((char *) addrRst[i], sizeRst[i], sizeof(double));
                                                          break;
                                                          case INT64_T:
                                                          swap_Nbyte((char *) addrRst[i], sizeRst[i], sizeof(int64_t));
                                                          break;
                                                        }
                                                      }
                                                    }

                                                    /* set local number of cells to be read by each process */
                                                    nk = nc / MPIsize;
                                                    nr = nc % MPIsize;
                                                    lnc = (MPIrank < nr ? nk + 1 : nk);

                                                    if (nc > maxCells)
                                                    stopRun(115, NULL, __FILE__, __LINE__);

                                                    h = 2.5 * csize;

                                                    h2 = h * h;
                                                    h3 = h2 * h;
                                                    h4 = h3 * h;

                                                    cellsAllocate();

                                                    /* gather number of cells from each process */
                                                    MPI_Allgather(&lnc, 1, MPI_INT64_T, tlnc, 1, MPI_INT64_T,
                                                      MPI_COMM_WORLD);
                                                      for (i = 0; i < MPIrank; i++)
                                                      nprev += tlnc[i];

                                                      goffset += nprev * sizeof(struct cellData);

                                                      MPI_File_seek(fh, goffset, MPI_SEEK_SET);
                                                      /* read cells data */
                                                      MPI_File_read(fh, cells, lnc * sizeof(struct cellData), MPI_BYTE,
                                                      MPI_STATUS_IGNORE);

                                                      if (swap) {
                                                        for (i = 0; i < lnc; i++) {
                                                          swap_Nbyte((char *) (&cells[i].gid), 1, sizeof(ZOLTAN_ID_TYPE));
                                                          swap_Nbyte((char *) (&cells[i].lifetime), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].phase), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].phasetime), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].g1), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].s), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].g2), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].m), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].x), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].y), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].z), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].size), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].h), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].v), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].density), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].age), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].death), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].young), 1, sizeof(float));
                                                          swap_Nbyte((char *) (&cells[i].tumor), 1, sizeof(unsigned char));
                                                          swap_Nbyte((char *) (&cells[i].scalarField), 1, sizeof(double));
                                                          swap_Nbyte((char *) (&cells[i].ctype), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].scstage), 1, sizeof(int));
                                                          swap_Nbyte((char *) (&cells[i].halo), 1, sizeof(int));
                                                        }
                                                      }

                                                      lg0nc = 0;
                                                      lg1nc = 0;
                                                      lsnc = 0;
                                                      lg2nc = 0;
                                                      lmnc = 0;
                                                      lcnc = 0;
                                                      lnnc = 0;
                                                      lvc = 0;
                                                      cancer = 0;

                                                      /* count cell types and phases */
                                                      for (i = 0; i < lnc; i++) {

                                                        cells[i].gid =
                                                        (unsigned long long int) MPIrank *(unsigned long long int)
                                                        maxCellsPerProc + (unsigned long long int) i;

                                                        switch (cells[i].phase) {
                                                          case 0:			/* G0 phase */
                                                          lg0nc++;
                                                          break;
                                                          case 1:			/* G1 phase */
                                                          lg1nc++;
                                                          break;
                                                          case 2:			/* S phase */
                                                          lsnc++;
                                                          break;
                                                          case 3:			/* G2 phase */
                                                          lg2nc++;
                                                          break;
                                                          case 4:			/* M phase */
                                                          lmnc++;
                                                          break;
                                                          case 5:			/* Necrotic cells */
                                                          lnnc++;
                                                          break;
                                                        }

                                                        if (cells[i].tumor == 1) {
                                                          lcnc++;
                                                          lcancer = 1;
                                                        }

                                                        if (cells[i].ctype == 1)
                                                        lvc++;

                                                      }

                                                      for(i=0; i<lnc; i++)
                                                      nscinst[cells[i].scstage]+=1;

                                                      localID = lnc;

                                                      randomStreamInit();

                                                      MPI_Allreduce(localCellCount, totalCellCount, numberOfCounts,
                                                        MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
                                                        MPI_Allreduce(nscinst,gnscinst,nscstages,
                                                          MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
                                                          MPI_Allreduce(&localbc,&globalbc,1,MPI_INT64_T,MPI_SUM,MPI_COMM_WORLD);
                                                          MPI_Allreduce(&lcancer, &cancer, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

                                                          nsteps = 512;			/* set by default in the case of the restart simulation */

                                                          /* close file */
                                                          MPI_File_close(&fh);

                                                          MPI_Barrier(MPI_COMM_WORLD);
                                                          //decompositionInit(argc, argv, MPI_COMM_WORLD);

                                                        }

                                                        /*!
                                                        * This function defines output for global fields.
                                                        */
                                                        void ioDefineOutputGlobalFields()
                                                        {
                                                          int f;
                                                          /* output fields */
                                                          for (f = 0; f < NFIELDS; f++) {
                                                            if(f==BVES && !bvsim) continue;
                                                            if(f==TEMP && !temperature) continue;
                                                            if(f==OXYG && !oxygen) continue;
                                                            if(f==GLUC && !glucose) continue;
                                                            if(f==HYDR && !hydrogenIon) continue;
                                                            strcpy(nameOut[f], fieldName[f]);
                                                            dimOut[f] = SCALAR;
                                                            typeOut[f] = REAL;
                                                            addrOut[f] = &fieldAddr[f];
                                                            jumpOut[f] = sizeof(double);
                                                          }
                                                          /* output gradient */
                                                          for(f=NFIELDS; f<NFIELDS+NCHEM; f++) {
                                                            if(f-NFIELDS==OXYG && !oxygen) continue;
                                                            if(f-NFIELDS==GLUC && !glucose) continue;
                                                            if(f-NFIELDS==HYDR && !hydrogenIon) continue;
                                                            strcpy(nameOut[f], fieldName[NGLOB+f-NFIELDS]);
                                                            sprintf(nameOut[f]+strlen(nameOut[f]),"Gradient");
                                                            dimOut[f] = VECTOR;
                                                            typeOut[f] = REAL;
                                                            addrOut[f] = &gradAddr[f];
                                                            jumpOut[f] = sizeof(double);
                                                          }
                                                        }

                                                        /*!
                                                        * This function prints global fields data in VisNow data format
                                                        * http://visnow.icm.edu.pl
                                                        */
                                                        void ioWriteFields(int step)
                                                        {
                                                          int j,f;
                                                          MPI_File fh1, fh2;
                                                          FILE *fh3;
                                                          struct floatVector3d *floatVectorField;
                                                          float *floatScalarField;
                                                          int *integerScalarField;
                                                          int64_t size;
                                                          int bdim;
                                                          int gsize[3];
                                                          int bsize[3];			/* box size */
                                                          int bstart[3];		/* box start */
                                                          MPI_Datatype subarray1_t, subarray2_t, float3_t;

                                                          if (!gfields)
                                                          return;

                                                          ioDefineOutputGlobalFields();

                                                          bdim = 3;
                                                          gsize[0] = gridI;
                                                          gsize[1] = gridJ;
                                                          gsize[2] = gridK;
                                                          bsize[0] = gridSize.x;
                                                          bsize[1] = gridSize.y;
                                                          bsize[2] = gridSize.z;
                                                          bstart[0] = gridStartIdx[MPIrank].x;
                                                          bstart[1] = gridStartIdx[MPIrank].y;
                                                          bstart[2] = gridStartIdx[MPIrank].z;

                                                          MPI_Type_vector(1, 3, 0, MPI_FLOAT, &float3_t);
                                                          MPI_Type_commit(&float3_t);

                                                          MPI_Type_create_subarray(bdim, gsize, bsize, bstart, MPI_ORDER_C,
                                                            float3_t, &subarray1_t);
                                                            MPI_Type_commit(&subarray1_t);

                                                            gsize[0] = gridI;
                                                            gsize[1] = gridJ;
                                                            gsize[2] = gridK;
                                                            bsize[0] = gridSize.x;
                                                            bsize[1] = gridSize.y;
                                                            bsize[2] = gridSize.z;
                                                            bstart[0] = gridStartIdx[MPIrank].x;
                                                            bstart[1] = gridStartIdx[MPIrank].y;
                                                            bstart[2] = gridStartIdx[MPIrank].z;

                                                            MPI_Type_create_subarray(bdim, gsize, bsize, bstart, MPI_ORDER_C,
                                                              MPI_FLOAT, &subarray2_t);
                                                              MPI_Type_commit(&subarray2_t);

                                                              for (f = 0; f < NFIELDS+NCHEM; f++) {

                                                                if(f==BVES && !bvsim) continue;
                                                                if(f==TEMP && !temperature) continue;
                                                                if((f==OXYG || f-NFIELDS+NGLOB==OXYG) && !oxygen) continue;
                                                                if((f==GLUC || f-NFIELDS+NGLOB==GLUC) && !glucose) continue;
                                                                if((f==HYDR || f-NFIELDS+NGLOB==HYDR) && !hydrogenIon) continue;

                                                                char fstname1[256];
                                                                char fstname2[256];
                                                                char fstname3[256];

                                                                size = gridSize.x * gridSize.y * gridSize.z;

                                                                floatVectorField =
                                                                (struct floatVector3d *) malloc(size *
                                                                  sizeof(struct floatVector3d));
                                                                  floatScalarField = (float *) malloc(size * sizeof(float));
                                                                  integerScalarField = (int *) malloc(size * sizeof(int));

                                                                  sprintf(fstname1, "%s/%s%08dcoords.bin", outdir, nameOut[f], step);
                                                                  sprintf(fstname2, "%s/%s%08dvalues.bin", outdir, nameOut[f], step);
                                                                  sprintf(fstname3, "%s/%s%08d.vnf", outdir, nameOut[f], step);

                                                                  if(MPIrank==0) {
                                                                    fh3=fopen(fstname3,"w");
                                                                    fprintf(fh3,"#VisNow regular field\n");
                                                                    fprintf(fh3,"field %s%08d, dim %ld %ld %ld, coords\n",nameOut[f],step,gridI,gridJ,gridK);
                                                                    if(dimOut[f]==SCALAR) fprintf(fh3,"component %s float\n",nameOut[f]);
                                                                    else fprintf(fh3,"component %s float, vector 3\n",nameOut[f]);
                                                                    fprintf(fh3,"file=%s%08dcoords.bin binary ",nameOut[f],step);
                                                                    if(!endian) fprintf(fh3,"big\n");
                                                                    else fprintf(fh3,"little\n");
                                                                    fprintf(fh3,"coords\n");
                                                                    fprintf(fh3,"file=%s%08dvalues.bin binary ",nameOut[f],step);
                                                                    if(!endian) fprintf(fh3,"big\n");
                                                                    else fprintf(fh3,"little\n");
                                                                    fprintf(fh3,"%s\n",nameOut[f]);
                                                                    fclose(fh3);
                                                                  }

                                                                  MPI_File_open(MPI_COMM_WORLD, fstname1,
                                                                    MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1);
                                                                    MPI_File_set_view(fh1, 0, MPI_FLOAT, subarray1_t, "native",
                                                                    MPI_INFO_NULL);
                                                                    /* truncate the first file */
                                                                    MPI_File_set_size(fh1, 0);

                                                                    for (j = 0; j < size; j++) {
                                                                      floatVectorField[j].x = (float) (gridBuffer[j].x);
                                                                      floatVectorField[j].y = (float) (gridBuffer[j].y);
                                                                      floatVectorField[j].z = (float) (gridBuffer[j].z);
                                                                    }
                                                                    if (!endian)
                                                                    swap_Nbyte((char *) floatVectorField, size * 3, sizeof(float));
                                                                    MPI_File_write(fh1, floatVectorField, 3 * size, MPI_FLOAT,
                                                                      MPI_STATUS_IGNORE);
                                                                      MPI_File_close(&fh1);

                                                                      MPI_File_open(MPI_COMM_WORLD, fstname2,
                                                                        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2);
                                                                        if(dimOut[f]==SCALAR)
                                                                        MPI_File_set_view(fh2, 0, MPI_FLOAT, subarray2_t, "native",
                                                                        MPI_INFO_NULL);
                                                                        else
                                                                        MPI_File_set_view(fh2, 0, MPI_FLOAT, subarray1_t, "native",
                                                                        MPI_INFO_NULL);

                                                                        /* truncate the second file */
                                                                        MPI_File_set_size(fh2, 0);

                                                                        if(dimOut[f]==SCALAR) {
                                                                          for (j = 0; j < size; j++)
                                                                          floatScalarField[j] = fieldAddr[f][j];
                                                                          if (!endian)
                                                                          swap_Nbyte((char *) floatScalarField, size, sizeof(float));
                                                                          MPI_File_write(fh2, floatScalarField, size, MPI_FLOAT,
                                                                            MPI_STATUS_IGNORE);
                                                                          } else {
                                                                            for (j = 0; j < size; j++) {
                                                                              floatVectorField[j].x = (float) (gradAddr[f-NFIELDS][j*3]);
                                                                              floatVectorField[j].y = (float) (gradAddr[f-NFIELDS][j*3+1]);
                                                                              floatVectorField[j].z = (float) (gradAddr[f-NFIELDS][j*3+2]);
                                                                              //if(floatVectorField[j].x!=0.0) printf("J:%d\n",j);
                                                                              //if(floatVectorField[j].y!=0.0) printf("J:%d\n",j);
                                                                              //if(floatVectorField[j].z!=0.0) printf("J:%d\n",j);
                                                                            }
                                                                            if (!endian)
                                                                            swap_Nbyte((char *) floatVectorField, size*3, sizeof(float));
                                                                            MPI_File_write(fh2, floatVectorField, size*3, MPI_FLOAT,
                                                                              MPI_STATUS_IGNORE);
                                                                            }

                                                                            MPI_File_close(&fh2);

                                                                            free(floatVectorField);
                                                                            free(floatScalarField);
                                                                            free(integerScalarField);

                                                                          }
                                                                          MPI_Type_free(&subarray1_t);
                                                                          MPI_Type_free(&subarray2_t);
                                                                        }

                                                                        /*!
                                                                        * This function redirects stdout to a given file.
                                                                        */
                                                                        void switchStdOut(const char *newStream)
                                                                        {
                                                                          fflush(stdout);
                                                                          dup2(fileno(stdout), fdSave);
                                                                          fflush(stdout);
                                                                          fdNew = open(newStream, O_WRONLY | O_CREAT | O_TRUNC, 0644);
                                                                          dup2(fdNew, fileno(stdout));
                                                                          close(fdNew);
                                                                        }

                                                                        /*!
                                                                        * This function brings back the stdout.
                                                                        */
                                                                        void revertStdOut()
                                                                        {
                                                                          fflush(stdout);
                                                                          dup2(fdSave, fileno(stdout));
                                                                          close(fdSave);
                                                                        }

                                                                        /*!
                                                                        * This function defines colormaps for PovRay output.
                                                                        */
                                                                        void defineColormaps()
                                                                        {
                                                                          int numberOfColormaps;

                                                                          /* on the ocassion - initialize the rotation angle */
                                                                          beta = 0.0;

                                                                          numberOfColormaps = 6;
                                                                          cmaps = (colormap *) malloc(numberOfColormaps * sizeof(colormap));

                                                                          /* medical colormap i=0 */
                                                                          strcpy(cmaps[0].name, "medical");
                                                                          cmaps[0].ncp = 4;
                                                                          cmaps[0].cp =
                                                                          (colormapPoint *) malloc(cmaps[0].ncp * sizeof(colormapPoint));
                                                                          cmaps[0].cp[0].position = 0.0;
                                                                          cmaps[0].cp[0].r = 0.0;
                                                                          cmaps[0].cp[0].g = 0.0;
                                                                          cmaps[0].cp[0].b = 0.0;
                                                                          cmaps[0].cp[1].position = 0.33;
                                                                          cmaps[0].cp[1].r = 122.0;
                                                                          cmaps[0].cp[1].g = 32.0;
                                                                          cmaps[0].cp[1].b = 32.0;
                                                                          cmaps[0].cp[2].position = 0.66;
                                                                          cmaps[0].cp[2].r = 255.0;
                                                                          cmaps[0].cp[2].g = 179.0;
                                                                          cmaps[0].cp[2].b = 77.0;
                                                                          cmaps[0].cp[3].position = 1.0;
                                                                          cmaps[0].cp[3].r = 255.0;
                                                                          cmaps[0].cp[3].g = 255.0;
                                                                          cmaps[0].cp[3].b = 255.0;

                                                                          /* rainbow colormap i=1 */
                                                                          strcpy(cmaps[1].name, "rainbow");
                                                                          cmaps[1].ncp = 5;
                                                                          cmaps[1].cp =
                                                                          (colormapPoint *) malloc(cmaps[1].ncp * sizeof(colormapPoint));
                                                                          cmaps[1].cp[0].position = 0.0;
                                                                          cmaps[1].cp[0].r = 0.0;
                                                                          cmaps[1].cp[0].g = 0.0;
                                                                          cmaps[1].cp[0].b = 255.0;
                                                                          cmaps[1].cp[1].position = 0.25;
                                                                          cmaps[1].cp[1].r = 0.0;
                                                                          cmaps[1].cp[1].g = 255.0;
                                                                          cmaps[1].cp[1].b = 255.0;
                                                                          cmaps[1].cp[2].position = 0.5;
                                                                          cmaps[1].cp[2].r = 0.0;
                                                                          cmaps[1].cp[2].g = 255.0;
                                                                          cmaps[1].cp[2].b = 0.0;
                                                                          cmaps[1].cp[3].position = 0.75;
                                                                          cmaps[1].cp[3].r = 255.0;
                                                                          cmaps[1].cp[3].g = 255.0;
                                                                          cmaps[1].cp[3].b = 0.0;
                                                                          cmaps[1].cp[4].position = 1.0;
                                                                          cmaps[1].cp[4].r = 255.0;
                                                                          cmaps[1].cp[4].g = 0.0;
                                                                          cmaps[1].cp[4].b = 0.0;

                                                                          /* blue red yellow */
                                                                          strcpy(cmaps[2].name, "bry");
                                                                          cmaps[2].ncp = 4;
                                                                          cmaps[2].cp =
                                                                          (colormapPoint *) malloc(cmaps[2].ncp * sizeof(colormapPoint));
                                                                          cmaps[2].cp[0].position = 0.0;
                                                                          cmaps[2].cp[0].r = 0.0;
                                                                          cmaps[2].cp[0].g = 0.0;
                                                                          cmaps[2].cp[0].b = 255.0;
                                                                          cmaps[2].cp[1].position = 0.33;
                                                                          cmaps[2].cp[1].r = 255.0;
                                                                          cmaps[2].cp[1].g = 0.0;
                                                                          cmaps[2].cp[1].b = 255.0;
                                                                          cmaps[2].cp[2].position = 0.67;
                                                                          cmaps[2].cp[2].r = 255.0;
                                                                          cmaps[2].cp[2].g = 0.0;
                                                                          cmaps[2].cp[2].b = 0.0;
                                                                          cmaps[2].cp[3].position = 1.0;
                                                                          cmaps[2].cp[3].r = 255.0;
                                                                          cmaps[2].cp[3].g = 255.0;
                                                                          cmaps[2].cp[3].b = 0.0;

                                                                          /* hot */
                                                                          strcpy(cmaps[3].name, "hot");
                                                                          cmaps[3].ncp = 5;
                                                                          cmaps[3].cp =
                                                                          (colormapPoint *) malloc(cmaps[3].ncp * sizeof(colormapPoint));
                                                                          cmaps[3].cp[0].position = 0.0;
                                                                          cmaps[3].cp[0].r = 107.0;
                                                                          cmaps[3].cp[0].g = 0.0;
                                                                          cmaps[3].cp[0].b = 0.0;
                                                                          cmaps[3].cp[1].position = 0.35;
                                                                          cmaps[3].cp[1].r = 255.0;
                                                                          cmaps[3].cp[1].g = 102.0;
                                                                          cmaps[3].cp[1].b = 28.0;
                                                                          cmaps[3].cp[2].position = 0.57;
                                                                          cmaps[3].cp[2].r = 250.0;
                                                                          cmaps[3].cp[2].g = 235.0;
                                                                          cmaps[3].cp[2].b = 128.0;
                                                                          cmaps[3].cp[3].position = 0.76;
                                                                          cmaps[3].cp[3].r = 232.0;
                                                                          cmaps[3].cp[3].g = 230.0;
                                                                          cmaps[3].cp[3].b = 230.0;
                                                                          cmaps[3].cp[4].position = 1.0;
                                                                          cmaps[3].cp[4].r = 156.0;
                                                                          cmaps[3].cp[4].g = 161.0;
                                                                          cmaps[3].cp[4].b = 255.0;

                                                                          /* hot1 */
                                                                          strcpy(cmaps[4].name, "hot1");
                                                                          cmaps[4].ncp = 5;
                                                                          cmaps[4].cp =
                                                                          (colormapPoint *) malloc(cmaps[4].ncp * sizeof(colormapPoint));
                                                                          cmaps[4].cp[0].position = 0.0;
                                                                          cmaps[4].cp[0].r = 128.0;
                                                                          cmaps[4].cp[0].g = 0.0;
                                                                          cmaps[4].cp[0].b = 0.0;
                                                                          cmaps[4].cp[1].position = 0.2;
                                                                          cmaps[4].cp[1].r = 255.0;
                                                                          cmaps[4].cp[1].g = 0.0;
                                                                          cmaps[4].cp[1].b = 0.0;
                                                                          cmaps[4].cp[2].position = 0.4;
                                                                          cmaps[4].cp[2].r = 255.0;
                                                                          cmaps[4].cp[2].g = 255.0;
                                                                          cmaps[4].cp[2].b = 0.0;
                                                                          cmaps[4].cp[3].position = 0.7;
                                                                          cmaps[4].cp[3].r = 255.0;
                                                                          cmaps[4].cp[3].g = 255.0;
                                                                          cmaps[4].cp[3].b = 255.0;
                                                                          cmaps[4].cp[4].position = 1.0;
                                                                          cmaps[4].cp[4].r = 128.0;
                                                                          cmaps[4].cp[4].g = 128.0;
                                                                          cmaps[4].cp[4].b = 255.0;

                                                                          /* my */
                                                                          strcpy(cmaps[5].name, "my");
                                                                          cmaps[5].ncp = 5;
                                                                          cmaps[5].cp =
                                                                          (colormapPoint *) malloc(cmaps[5].ncp * sizeof(colormapPoint));
                                                                          cmaps[5].cp[0].position = 0.0;
                                                                          cmaps[5].cp[0].r = 107.0;
                                                                          cmaps[5].cp[0].g = 0.0;
                                                                          cmaps[5].cp[0].b = 0.0;
                                                                          cmaps[5].cp[1].position = 0.35;
                                                                          cmaps[5].cp[1].r = 0.0;
                                                                          cmaps[5].cp[1].g = 100.0;
                                                                          cmaps[5].cp[1].b = 255.0;
                                                                          cmaps[5].cp[2].position = 0.57;
                                                                          cmaps[5].cp[2].r = 250.0;
                                                                          cmaps[5].cp[2].g = 235.0;
                                                                          cmaps[5].cp[2].b = 128.0;
                                                                          cmaps[5].cp[3].position = 0.76;
                                                                          cmaps[5].cp[3].r = 232.0;
                                                                          cmaps[5].cp[3].g = 230.0;
                                                                          cmaps[5].cp[3].b = 230.0;
                                                                          cmaps[5].cp[4].position = 1.0;
                                                                          cmaps[5].cp[4].r = 156.0;
                                                                          cmaps[5].cp[4].g = 161.0;
                                                                          cmaps[5].cp[4].b = 255.0;

                                                                        }


                                                                        /*!
                                                                        * This function prints PovRay output with cellular data.
                                                                        */
                                                                        void ioWriteStepPovRay(int step, int type)
                                                                        {
                                                                          int c;
                                                                          int i;
                                                                          char fstname[256];
                                                                          MPI_File fh;
                                                                          MPI_Status status;
                                                                          MPI_Datatype subarray_t;
                                                                          float cr, cg, cb;
                                                                          char *const fmt =
                                                                          "sphere{ <%10.4lf,%10.4lf,%10.4lf>,%10.4lf texture { pigment { color rgb <%6.4f,%6.4f,%6.4f> } finish { phong 0.2 ambient .1 }}}         \n";
                                                                          char testBuffer[512];
                                                                          char *txtData;
                                                                          char *txtData_p;
                                                                          char *txtHeaderData;
                                                                          int numCharsPerCell;
                                                                          int headerLen;
                                                                          int headerSize;
                                                                          int error;
                                                                          int gdims;
                                                                          int gsize[1];
                                                                          int istart[1];
                                                                          int isize[1];
                                                                          int64_t printed = 0;
                                                                          int64_t tPrinted[MPIsize];
                                                                          double fmin, fmax, fepsilon;
                                                                          int cm;
                                                                          int cmReverse = 0;
                                                                          int cmShift = 0;
                                                                          double cmPad;
                                                                          double minCorner[3], maxCorner[3];
                                                                          double gMinCorner[3], gMaxCorner[3];

                                                                          double middlePointLocal[3];
                                                                          double middlePointGlobal[3];
                                                                          double lmass, gmass;

                                                                          /* type: 0 - denisty, 1 - oxygen, 2 - phases, 3 - slice & phases */

                                                                          minCorner[0] = DBL_MAX;
                                                                          minCorner[1] = DBL_MAX;
                                                                          minCorner[2] = DBL_MAX;
                                                                          maxCorner[0] = -DBL_MAX;
                                                                          maxCorner[1] = -DBL_MAX;
                                                                          maxCorner[2] = -DBL_MAX;

                                                                          for (i = 0; i < lnc; i++) {
                                                                            minCorner[0] =
                                                                            (cells[i].x - cells[i].size <
                                                                              minCorner[0] ? cells[i].x - cells[i].size : minCorner[0]);
                                                                              maxCorner[0] =
                                                                              (cells[i].x + cells[i].size >
                                                                                maxCorner[0] ? cells[i].x + cells[i].size : maxCorner[0]);
                                                                                minCorner[1] =
                                                                                (cells[i].y - cells[i].size <
                                                                                  minCorner[1] ? cells[i].y - cells[i].size : minCorner[1]);
                                                                                  maxCorner[1] =
                                                                                  (cells[i].y + cells[i].size >
                                                                                    maxCorner[1] ? cells[i].y + cells[i].size : maxCorner[1]);
                                                                                    minCorner[2] =
                                                                                    (cells[i].z - cells[i].size <
                                                                                      minCorner[2] ? cells[i].z - cells[i].size : minCorner[2]);
                                                                                      maxCorner[2] =
                                                                                      (cells[i].z + cells[i].size >
                                                                                        maxCorner[2] ? cells[i].z + cells[i].size : maxCorner[2]);
                                                                                      }
                                                                                      MPI_Allreduce(minCorner, gMinCorner, 3, MPI_DOUBLE, MPI_MIN,
                                                                                        MPI_COMM_WORLD);
                                                                                        MPI_Allreduce(maxCorner, gMaxCorner, 3, MPI_DOUBLE, MPI_MAX,
                                                                                          MPI_COMM_WORLD);

                                                                                          middlePointLocal[0] = 0.0;
                                                                                          middlePointLocal[1] = 0.0;
                                                                                          middlePointLocal[2] = 0.0;
                                                                                          lmass = 0.0;
                                                                                          gmass = 0.0;

                                                                                          for (c = 0; c < lnc; c++) {
                                                                                            middlePointLocal[0] += cells[c].size * cells[c].x;
                                                                                            middlePointLocal[1] += cells[c].size * cells[c].y;
                                                                                            middlePointLocal[2] += cells[c].size * cells[c].z;
                                                                                            lmass += cells[c].size;
                                                                                          }

                                                                                          MPI_Allreduce(middlePointLocal, middlePointGlobal, 3, MPI_DOUBLE,
                                                                                            MPI_SUM, MPI_COMM_WORLD);
                                                                                            MPI_Allreduce(&lmass, &gmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                                                                                            middlePointGlobal[0] /= gmass;
                                                                                            middlePointGlobal[1] /= gmass;
                                                                                            middlePointGlobal[2] /= gmass;

                                                                                            middlePointGlobal[0] =
                                                                                            gMinCorner[0] + (gMaxCorner[0] - gMinCorner[0]) / 2;
                                                                                            middlePointGlobal[1] =
                                                                                            gMinCorner[1] + (gMaxCorner[1] - gMinCorner[1]) / 2;
                                                                                            middlePointGlobal[2] =
                                                                                            gMinCorner[2] + (gMaxCorner[2] - gMinCorner[2]) / 2;

                                                                                            /* write data to test buffer */
                                                                                            cr = 0.0;
                                                                                            cg = 0.0;
                                                                                            cb = 0.0;
                                                                                            cm = 0;
                                                                                            numCharsPerCell =
                                                                                            sprintf(testBuffer, fmt, cells[1].x, cells[1].y, cells[1].z,
                                                                                              cells[1].size, cr, cg, cb);
                                                                                              if (!
                                                                                                (txtData =
                                                                                                  (char *) malloc(numCharsPerCell * (lnc + 16) * sizeof(char))))
                                                                                                  stopRun(106, "txtData", __FILE__, __LINE__);
                                                                                                  txtData_p = txtData;

                                                                                                  for (i = 0; i < MPIsize; i++)
                                                                                                  tPrinted[i] = 0;

                                                                                                  switch (type) {
                                                                                                    case 0:
                                                                                                    sprintf(fstname, "%s/step%08ddensity.pov", outdir, step);
                                                                                                    cm = 5;
                                                                                                    cmReverse = 1;
                                                                                                    break;
                                                                                                    case 1:
                                                                                                    sprintf(fstname, "%s/step%08doxygen.pov", outdir, step);
                                                                                                    MPI_Allreduce(&fieldMin[OXYG], &fmin, 1, MPI_DOUBLE, MPI_MIN,
                                                                                                      MPI_COMM_WORLD);
                                                                                                      MPI_Allreduce(&fieldMax[OXYG], &fmax, 1, MPI_DOUBLE, MPI_MAX,
                                                                                                        MPI_COMM_WORLD);
                                                                                                        fepsilon = (fmax - fmin) * 0.1;
                                                                                                        fmin -= fepsilon;
                                                                                                        fmax += fepsilon;
                                                                                                        cm = 0;
                                                                                                        break;
                                                                                                        case 2:
                                                                                                        sprintf(fstname, "%s/step%08dphases.pov", outdir, step);
                                                                                                        cm = 1;
                                                                                                        break;
                                                                                                        case 3:
                                                                                                        sprintf(fstname, "%s/step%08dslice.pov", outdir, step);
                                                                                                        cm = 1;
                                                                                                        break;
                                                                                                        case 4:
                                                                                                        sprintf(fstname, "%s/step%08doutside.pov", outdir, step);
                                                                                                        cm = 0;
                                                                                                        cmReverse = 1;
                                                                                                        cmShift = 1;
                                                                                                        cmPad = 0.15;
                                                                                                        break;
                                                                                                      }

                                                                                                      /* open file */
                                                                                                      error =
                                                                                                      MPI_File_open(MPI_COMM_WORLD, fstname,
                                                                                                        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
                                                                                                        if (error != MPI_SUCCESS)
                                                                                                        if (MPIrank == 0)
                                                                                                        stopRun(113, NULL, __FILE__, __LINE__);

                                                                                                        MPI_File_set_size(fh, 0);

                                                                                                        /* write header */
                                                                                                        if (MPIrank == 0) {

                                                                                                          float cameraLocation[3];
                                                                                                          float lookAt[3];
                                                                                                          float lightSource1[3];
                                                                                                          float lightSource2[3];
                                                                                                          float a1, a2, a3, b1, b2, b3, dist;
                                                                                                          float ss, cc;
                                                                                                          float corner;

                                                                                                          a1 = (gMaxCorner[0] - gMinCorner[0]) / 2;
                                                                                                          a2 = (gMaxCorner[1] - gMinCorner[1]) / 2;
                                                                                                          a3 = (gMaxCorner[2] - gMinCorner[2]) / 2;

                                                                                                          b1 = a1 / (tan(30.0 * M_PI / 180.0));
                                                                                                          b2 = a2 / (tan(30.0 * M_PI / 180.0));
                                                                                                          b3 = a3 / (tan(30.0 * M_PI / 180.0));

                                                                                                          dist = (b1 >= b2 ? b1 : b2);
                                                                                                          dist = (dist >= b3 ? dist : b3);

                                                                                                          ss = sin(beta * M_PI / 180.0);
                                                                                                          cc = cos(beta * M_PI / 180.0);

                                                                                                          corner = sqrt(pow(a1, 2) + pow(a3, 2));

                                                                                                          cameraLocation[0]=middlePointGlobal[0]-2;
                                                                                                          cameraLocation[1]=middlePointGlobal[1]-0.8*dist;
                                                                                                          cameraLocation[2]=middlePointGlobal[2]-5;

                                                                                                          /*cameraLocation[0] =
                                                                                                          -(middlePointGlobal[2] - corner - dist -
                                                                                                          middlePointGlobal[2]) * ss + middlePointGlobal[0];
                                                                                                          cameraLocation[1] = middlePointGlobal[1];
                                                                                                          cameraLocation[2] =
                                                                                                          (middlePointGlobal[2] - corner - dist -
                                                                                                          middlePointGlobal[2]) * cc + middlePointGlobal[2];*/

                                                                                                          lookAt[0] = middlePointGlobal[0]-2;
                                                                                                          lookAt[1] = middlePointGlobal[1];
                                                                                                          lookAt[2] = middlePointGlobal[2]-5;

                                                                                                          /*    lookAt[0] = middlePointGlobal[0];
                                                                                                          lookAt[1] = middlePointGlobal[1];
                                                                                                          lookAt[2] = middlePointGlobal[2];
                                                                                                          */
                                                                                                          lightSource1[0] = cameraLocation[0];
                                                                                                          lightSource1[1] = cameraLocation[1];
                                                                                                          lightSource1[2] = cameraLocation[2];

                                                                                                          lightSource2[0] = lightSource1[0];
                                                                                                          lightSource2[1] = lightSource1[1];
                                                                                                          lightSource2[2] = lightSource1[2];

                                                                                                          txtHeaderData = (char *) malloc(1024 * sizeof(char));

                                                                                                          headerLen =
                                                                                                          sprintf(txtHeaderData,
                                                                                                            "#include \"colors.inc\"\ncamera { location <%f,%f,%f> look_at <%f,%f,%f> angle 60}\nlight_source { <%f,%f,%f> color White }\nlight_source { <%f,%f,%f> color White }\nbackground { color White }\n",
                                                                                                            cameraLocation[0], cameraLocation[1], cameraLocation[2],
                                                                                                            lookAt[0], lookAt[1], lookAt[2], lightSource1[0],
                                                                                                            lightSource1[1], lightSource1[2], lightSource2[0],
                                                                                                            lightSource2[1], lightSource2[2]);

                                                                                                            headerSize = headerLen * sizeof(char);

                                                                                                            MPI_File_write(fh, txtHeaderData, headerLen, MPI_CHAR, &status);

                                                                                                            free(txtHeaderData);
                                                                                                          }

                                                                                                          MPI_Bcast(&headerSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

                                                                                                          for (c = 0; c < lnc; c++) {
                                                                                                            float color = 0.0;
                                                                                                            int jump;
                                                                                                            if(beta==360.0 && (cells[c].y<1.8 || cells[c].y>4.5)) continue;
                                                                                                            /*    if (beta == 360.0
                                                                                                            && (cells[c].z > middlePointGlobal[2] + 4.0 * csize
                                                                                                            || cells[c].z < middlePointGlobal[2] - 4.0 * csize))
                                                                                                            continue; */
                                                                                                            if (type == 0) {
                                                                                                              color = ((cells[c].density) / (8*densityCriticalLevel2)); /* UWAGA! BUG!!! */
                                                                                                              if (cells[c].tumor)
                                                                                                              color = 0.0;
                                                                                                            }				//2.5);
                                                                                                            if (type == 1)
                                                                                                            color = ((cellFields[OXYG][c] - fmin) / (fmax - fmin));
                                                                                                            if (type == 2 || type == 3 || type == 4)
                                                                                                            color = (((float) cells[c].phase) / 5.0);
                                                                                                            if (type == 2)
                                                                                                            color = (((float) MPIrank) / 512.0);

                                                                                                            if (type == 0)
                                                                                                            color = color / 4 + 0.75;

                                                                                                            if (cmReverse)
                                                                                                            color = 1.0 - color;
                                                                                                            if (cmShift) {
                                                                                                              color = color * (1 - cmPad);
                                                                                                              if (color == 1 - cmPad)
                                                                                                              color = 1.0;
                                                                                                            }

                                                                                                            if(color<0.0) color=0.0;

                                                                                                            if(cells[c].ctype==0) color=0.2;
                                                                                                            else color=0.0;

                                                                                                            for (i = 1; i < cmaps[cm].ncp; i++) {
                                                                                                              float d, dr, dg, db;
                                                                                                              d = cmaps[cm].cp[i].position - cmaps[cm].cp[i - 1].position;
                                                                                                              dr = cmaps[cm].cp[i].r - cmaps[cm].cp[i - 1].r;
                                                                                                              dg = cmaps[cm].cp[i].g - cmaps[cm].cp[i - 1].g;
                                                                                                              db = cmaps[cm].cp[i].b - cmaps[cm].cp[i - 1].b;
                                                                                                              if (color <= cmaps[cm].cp[i].position) {
                                                                                                                cr = cmaps[cm].cp[i - 1].r +
                                                                                                                ((color - cmaps[cm].cp[i - 1].position) / d) * dr;
                                                                                                                cg = cmaps[cm].cp[i - 1].g +
                                                                                                                ((color - cmaps[cm].cp[i - 1].position) / d) * dg;
                                                                                                                cb = cmaps[cm].cp[i - 1].b +
                                                                                                                ((color - cmaps[cm].cp[i - 1].position) / d) * db;
                                                                                                                break;
                                                                                                              }

                                                                                                            }

                                                                                                            cr /= 255.0;
                                                                                                            cg /= 255.0;
                                                                                                            cb /= 255.0;

                                                                                                            jump =
                                                                                                            sprintf(txtData_p, fmt, cells[c].x, cells[c].y, cells[c].z,
                                                                                                              cells[c].size, cr, cg, cb);
                                                                                                              txtData_p += jump;
                                                                                                              printed += 1;

                                                                                                            }

                                                                                                            if (printed == 0) {
                                                                                                              float coords[3];
                                                                                                              float size = 0.0;
                                                                                                              float color[3];
                                                                                                              int jump;
                                                                                                              // artificial cell far away from the scene
                                                                                                              coords[0] = -512.0;
                                                                                                              coords[1] = -512.0;
                                                                                                              coords[2] = -512.0;
                                                                                                              color[0] = 0.0;
                                                                                                              color[1] = 0.0;
                                                                                                              color[2] = 0.0;

                                                                                                              printed = 1;

                                                                                                              jump =
                                                                                                              sprintf(txtData_p, fmt, coords[0], coords[1], coords[2], size,
                                                                                                                color[0], color[1], color[2]);
                                                                                                                txtData_p += jump;
                                                                                                              }

                                                                                                              gdims = 1;
                                                                                                              tPrinted[MPIrank] = printed;
                                                                                                              MPI_Allgather(&printed, 1, MPI_INT64_T, tPrinted, 1, MPI_INT64_T,
                                                                                                                MPI_COMM_WORLD);

                                                                                                                gsize[0] = 0;
                                                                                                                istart[0] = 0;
                                                                                                                isize[0] = printed * numCharsPerCell;
                                                                                                                for (i = 0; i < MPIrank; i++)
                                                                                                                istart[0] += tPrinted[i] * numCharsPerCell;

                                                                                                                for (i = 0; i < MPIsize; i++)
                                                                                                                gsize[0] += tPrinted[i] * numCharsPerCell;

                                                                                                                MPI_Type_create_subarray(gdims, gsize, isize, istart, MPI_ORDER_C,
                                                                                                                  MPI_CHAR, &subarray_t);
                                                                                                                  MPI_Type_commit(&subarray_t);

                                                                                                                  MPI_File_set_view(fh, headerSize, MPI_CHAR, subarray_t, "native",
                                                                                                                  MPI_INFO_NULL);

                                                                                                                  MPI_File_write(fh, txtData, isize[0], MPI_CHAR, &status);

                                                                                                                  free(txtData);

                                                                                                                  MPI_File_close(&fh);
                                                                                                                  MPI_Type_free(&subarray_t);

                                                                                                                  if (beta < 360.0)
                                                                                                                  beta += 1.0;
                                                                                                                  else
                                                                                                                  beta = 360.0;

                                                                                                                }
