###############################################################################
# # g95
# FC=g95
# FFLAGS=
# F90=g95
# F90FLAGS=
# LD=$(FC)
# LDFLAGS=$(FFLAGS)
# CPP=cpp
# CPPFLAGS= -P -C -traditional
# LIBS=
###############################################################################
# ifort
# FC=ifort
# FFLAGS=
# F90=ifort
# F90FLAGS=
# LD=$(FC)
# LDFLAGS=$(FFLAGS)
# CPP=/usr/bin/cpp
# CPPFLAGS= -P -C -traditional
# LIBS=
###############################################################################
# ifort + MPI
FC=mpif90
FFLAGS=
F90=mpif90
F90FLAGS=
LD=$(FC)
LDFLAGS=$(FFLAGS)
CPP=/usr/bin/cpp
CPPFLAGS= -P -C -traditional -D__PARA
LIBS=
###############################################################################

.SUFFIXES: .o .f .f90 .F90

OBJS=m_vdw_phi.o vdwphi_drsll.o
PROG=vdwphigen_drsll

all: $(PROG)

$(PROG):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $@ $(LIBS)

.F90.o:
	$(CPP) $(CPPFLAGS) $*.F90 > $*.f90; $(F90) -c $(F90FLAGS) $*.f90; \
        rm -f $*.f90
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(PROG) *.mod *.o
