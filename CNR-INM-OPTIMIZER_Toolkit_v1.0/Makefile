# Intel
FC       = mpiifort 
FCOPT    = -O5
FCOPTOMP = -qopenmp
LIB      = -mkl

# GNU
#FC       = mpif90 
#FCOPT    = -O3
#FCOPTOMP = -fopenmp
#LIB      = -llapack -lblas

export FC
export FCOPT
export FCOPTOMP
export LIB

all:
	cd src/optimizer/        ; make; cd ../../
	cd src/parallel/         ; make; cd ../../
	cd src/sampling/         ; make; cd ../../
	$(FC) $(FCOPT) $(FCOPTOMP) -o bin/POT_v1.0 \
	       		              src/*.f90 \
 			              lib/liboptimizer.a \
                                      lib/libparallel.a \
			              lib/libsampling.a \
			              $(LIB)

clean:
	rm src/*/*.o
	rm lib/*.a
	rm bin/*


