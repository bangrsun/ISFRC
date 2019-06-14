# Makefile for building Green function
# 
# Modified by Shuying SUN for EXL 2019/05/24 
#
MACHINE := $(shell uname)
VPATH =../src
RUN_DIR = ../run

ifeq ($(MACHINE), Darwin)
  FC = gfortran
  FFLAGS = -g -Wall -fcheck=all
  LDFLAGS =
  LIB =
endif
ifeq ($(MACHINE),Linux)
  FC = pgf90
  FFLAGS = -g
  LDFLAGS = -Mlfs
  LIB =
endif

SOURCES = modules.f90 elliptic.f90 utilities.f90 \
					main.f90 setup.f90 load.f90 subroutines.f90

OBJECTS=$(SOURCES:.f90=.o)

gfun:$(OBJECTS) 
	$(FC) $(LDFLAGS) -o $@ $^ $(LIB)
	cp $@ $(RUN_DIR)
	@echo "GFUN BUILT"

$(filter-out modules.o,$(OBJECTS)):modules.o

main.o: main.f90 setup.o load.o subroutines.o
	$(FC) $(FFLAGS) -c $< $(LIB)

utilities.o:utilities.f90 elliptic.o
	$(FC) $(FFLAGS) -c $< $(LIB)

%.o:%.f90
	$(FC) $(FFLAGS) -c $< $(LIB)

.PHONY: clean

clean:
	rm *.o *.mod gfun
	@echo "CLEAN DONE"

