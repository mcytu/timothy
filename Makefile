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

.PHONY: all clean doc

all:
	@ if command -v gmake >/dev/null 2>&1; then \
	echo "Using gmake to build Timothy"; \
	cd src/; gmake; \
        else \
        echo "Using make to build Timothy"; \
        cd src/; make; \
        fi
clean:
	@ if command -v gmake >/dev/null 2>&1; then \
        echo "Using gmake to clean Timothy"; \
        cd src/; gmake clean; \
	cd ../doc; rm -fr html/ latex/; \
        else \
        echo "Using make to clean Timothy"; \
        cd src/; make clean; \
	cd ../doc; rm -fr html/ latex/; \
        fi
doc:
	@ if command -v doxygen >/dev/null 2>&1; then \
	    if command -v gmake >/dev/null 2>&1; then \
	      cd src/; gmake doc;\
	    else \
	      cd src/; make doc; \
	    fi \
	else \
	  echo "Doxygen not found"; \
	fi 
