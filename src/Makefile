#* **************************************************************************
# * This file is part of Timothy
# *
# * Copyright (c) 2014/15 Maciej Cytowski
# * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# *
# * *************************************************************************/

#SYSTYPE = $(shell hostname)
#SYSTYPE = ibm-p775
#SYSTYPE = ibm-bgq
#SYSTYPE = intel-mic
SYSTYPE = gnu-x86

# My laptop configuration (gnu-x86)
ZOLTANLIB=-L/opt/zoltan/lib -lzoltan
SPRNGLIB=-L/opt/sprng/2.0/lib -lsprng
GMPLIB=-lgmp
HYPRELIB=-L/opt/hypre/lib -lHYPRE -lblas -llapack
ZOLTANINC=-I/opt/zoltan/include
SPRNGINC=-I/opt/sprng/2.0/include
HYPREINC=-I/opt/hypre/include
MPIINC=-I/usr/include/mpi

# Nostromo (ibm-bgq) system configuration
#ZOLTANLIB=-L/opt/zoltan/3.8/lib -lzoltan
#SPRNGLIB=-L/opt/sprng/2.0/lib -lsprng
#GMPLIB=-L/opt/gmp/5.1.1/lib -lgmp
#HYPRELIB=-L/opt/hypre/2.9.0b_bigint_omp/lib -lHYPRE
#ZOLTANINC=-I/opt/zoltan/3.8/include
#SPRNGINC=-I/opt/sprng/2.0/include
#HYPREINC=-I/opt/hypre/2.9.0b_bigint_omp/include
#MPIINC=

# MICLAB system (intel-mic) configuration
#ZOLTANLIB=-L/home/mcytowski/libs/zoltan/lib -lzoltan
#SPRNGLIB=-L/home/mcytowski/libs/sprng/lib -lsprng
#GMPLIB=-L/home/mcytowski/libs/gmp/lib -lgmp
#HYPRELIB=-L/home/mcytowski/libs/hypre/lib -lHYPRE
#ZOLTANINC=-I/home/mcytowski/libs/zoltan/include
#SPRNGINC=-I/home/mcytowski/libs/sprng/include
#HYPREINC=-I/home/mcytowski/libs/hypre/include
#MPIINC=

# Boreasz (ibm-p775) system configuration 
#ZOLTANLIB=-L/p7data/opt/zoltan/3.6/lib -lzoltan
#SPRNGLIB=-L/p7data/opt/sprng/2.0/lib -lsprng
#GMPLIB=-L/p7data/opt/gmp/5.0.5/lib -lgmp
#HYPRELIB=-L/opt/hypre/2.9.0b/bigint_omp/lib -lHYPRE
#ZOLTANINC=-I/p7data/opt/zoltan/3.6/include
#SPRNGINC=-I/p7data/opt/sprng/2.0/include
#HYPREINC=-I/opt/hypre/2.9.0b/bigint_omp/include
#MPIINC=

#############################################################################
# Prepare configuration for your system 
#
#ifeq ($(SYSTYPE),HOSTNAME_OF_YOUR_SYSTEM)
#CC = mpicc
#OPT = -O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline -g
#DEFS =
#INCLUDES = -I/usr/include/mpi -I/usr/include/sprng -I. -I/opt/zoltan/3.81/include -I/opt/sprng/include -I/opt/hypre/2.9.0b/include
#LIBS = -lm -L/opt/zoltan/3.81/lib -lzoltan -L/opt/sprng/lib -lsprng -lgmp -L/opt/hypre/2.9.0b/lib -lHYPRE -lblas -llapack
#endif
#############################################################################

# IBM BG/Q configuration
ifeq ($(SYSTYPE),ibm-bgq)
CC = mpixlc_r
OPT = -O3 -qstrict -qtune=qp -qarch=qp -qsimd=auto -qsmp=omp -g -qinline -qkeyword=inline
DEFS = -DPROF
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB) -L/opt/ibmmath/essl/5.1/lib64 -lesslbg -L/opt/ibmcmp/xlf/bg/14.1/lib64 -lxlf90_r -lxlsmp
endif

# General GNU x86 configuration
ifeq ($(SYSTYPE),gnu-x86)
CC = mpicc -g #-DDEBUG
OPT = -O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline -g
DEFS = #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)
endif

# Intel compiler, Intel MIC architecture
ifeq ($(SYSTYPE),intel-mic)
CC = mpiicc 
OPT = -O3 -openmp 
DEFS = #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)
endif

# IBM Power775 configuration
ifeq ($(SYSTYPE),ibm-p775)
CC = mpcc_r -q64
OPT = -O3 -qinline -qunroll -qhot -qsimd=auto -qtune=pwr7 -qarch=pwr7 -qsmp=omp -qcpluscmt -DPROF -qkeyword=inline #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lmassv -lmass -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB) -lessl -lxlf90_r -lxlsmp
endif

INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)


EXEC = ../timothy

OPTIONS = $(OPT) $(DEFS)

OBJS = main.o io.o utils.o cells.o domdec.o octree.o potential.o comm.o stats.o random.o init.o grid.o interp.o fields.o tempf.o chemf.o vessel.o bone.o compute.o

INCL = global.h io.h inline.h fields.h potential.h octree.h

CFLAGS = $(OPTIONS) $(INCLUDES) 

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS); rm -f $(EXEC)

doc:
	doxygen Doxyfile


