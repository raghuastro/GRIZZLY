SHELL = /bin/sh

FC=mpif90

LDLIBS= -L/disk/dawn-1/gragh/software/fftw3/lib/ -I/disk/dawn-1/gragh/software/fftw3/include/ -lfftw3 -lm 
FFLAGS= -O3 -fPIC -mcmodel=large #-g -traceback -check all -fp-stack-check


PP= param.o adaptint.o sort.o func.o
OBJS= subr_support.o subr_readfile.o read_halofile.o cal_ion.o subr_overlap.o calpha.o subr_grizzly.o grizzly_main.o

run: GR.exe


GR.exe: $(PP) $(OBJS)
	$(FC) $(FFLAGS) -openmp  $^ -o $@ $(LDLIBS)

$(OBJS): $(PP)

%.o: %.f90
	$(FC) $(FFLAGS)  -c $< $(LDLIBS)

clean:
	rm -f *.o *.mod *~ GR.exe *nfs000*






