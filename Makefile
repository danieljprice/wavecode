.KEEP_STATE:

FORTRAN = gfortran-mp-4.7
FC = gfortran-mp-4.7

FFLAGS = -O3 -g -Wall -ffast-math -finline-limit=1500 -ftree-vectorize -funroll-loops

SOURCES= waveutils.f90 binary.f90 blackhole.f90 wavecode.f90

OBJECTS = $(SOURCES:.f90=.o) 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

all: wave

wave: $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)
	dsymutil $@

clean:
	\rm -f *.o *.mod wave
	\rm -rf wave.dSYM