LIB   =optimizer
LNAME =lib$(LIB).a
LOBJS =optimizer.o SODPSO.o DPO.o linesearch.o MODPSO.o closest.o dominance.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

optimizer.o: optimizer.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c optimizer.f90

SODPSO.o: SODPSO.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c SODPSO.f90

DPO.o: DPO.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c DPO.f90

linesearch.o: linesearch.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c linesearch.f90

MODPSO.o: MODPSO.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c MODPSO.f90

closest.o: closest.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c closest.f90

dominance.o: dominance.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c dominance.f90


