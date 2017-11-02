# ***************************************************************************
# Timothy - Tissue Modelling Framework
# Copyright (C) 2014-2018 Maciej Cytowski
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# ***************************************************************************

CC=mpicc
CFLAGS=-O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline
LDFLAGS = -fopenmp
DEFS=-DDEBUG

ZOLTANLIB=-L/opt/zoltan/lib -lzoltan
ZOLTANINC=-I/opt/zoltan/include

HYPRELIB=-L/opt/hypre/lib -lHYPRE -lblas -llapack
HYPREINC=-I/opt/hypre/include

GMPLIB=-lgmp
SPRNGLIB=-L/opt/sprng/2.0/lib -lsprng
SPRNGINC=-I/opt/sprng/2.0/include

CFLAGS += -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC)
LDFLAGS += -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)

OBJ = main.o io.o utils.o cells.o lb.o octree.o potential.o comm.o stats.o random.o init.o grid.o interp.o fields.o tempf.o chemf.o compute.o
DEPS = global.h io.h inline.h fields.h potential.h octree.h

all: timothy

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(DEFS) -o $@ -c $<

timothy: $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $@

clean:
	rm -f $(OBJ); rm -f timothy

doc:
	doxygen Doxyfile
