# Makefile for building EFUND
# 
# Modified by Shuying SUN for EXL 2019/05/24 
#
# Revisions:
# $Log: Makefile,v $
# Revision 1.1.2.2  2008/12/10 23:08:33  radhakri
# *** empty log message ***
#
# Revision 1.1.2.1  2008/11/18 22:22:46  radhakri
# *** empty log message ***
#
# Revision 1.1  2007/06/01 05:19:19  renq
#
# Makefile: added by Qilong Ren
#   Makefile for EFUND.
#
# Revision 1.2  2007/05/05 05:28:13  renq
# using make pattern rule
#
# Revision 1.1  2007/04/24 06:03:27  renq
# Initial revision
#
#srini added few lines to include machine dependency
MACHINE := $(shell uname)
VPATH =../src
RUN_DIR = ../run

#default HP-Unix
FFLAGS = -K +U77 +O1 -w -v -c
LDFLAGS= -K +O1 -w -v -o
FC=f90
FCOM=

ifeq ($(MACHINE),Linux)
  FC=pgf90
  FFLAGS = -Mbyteswapio -fast -c
  LDFLAGS = -Mlfs -o
  ifeq ($(FCOM),lf95)
    FC=/usr/local/lf9562/bin/lf95
    FFLAGS = -c
    LDFLAGS = -o
  endif
endif
ifeq ($(MACHINE), Darwin)
  FC=ifort
  FFLAGS= -c -g
  LDFLAGS= -o
  LIB =
	# FC=pgf90
	# FFLAGS = -Mbyteswapio -fast -c
	# LDFLAGS = -Mlfs -o
	# LIB = 
endif

MODULES = exparm.f90 bfgrid.f90 cacoil.f90 cecoil.f90 filech.f90 \
          coilsp.f90 consta.f90 cvesel.f90 fcoil.f90 vas_sub1.f90 \
          fshift.f90 input.f90 mprobe.f90 nio.f90 \
          pmodel.f90 rogowl.f90 siloop.f90 efundud.f90

OBJECTS=$(MODULES:.f90=.o)

efundud:$(OBJECTS) 
	$(FC) $(LDFLAGS) $@ $^ $(LIB)
	cp $@ $(RUN_DIR)
	@echo "EFUND BUILT"

%.o:%.f90
	$(FC) $(FFLAGS) $< $(LIB)

clean:
	rm *.o *.mod efundud
	@echo "CLEAN DONE"
