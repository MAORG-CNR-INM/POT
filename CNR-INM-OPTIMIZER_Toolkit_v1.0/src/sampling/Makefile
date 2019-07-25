LIB   =sampling
LNAME =lib$(LIB).a
LOBJS =doe.o prime.o D_HSS.o B_HSS.o Ortho_Init.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

doe.o: doe.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c doe.f90

prime.o: prime.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c prime.f90

D_HSS.o: D_HSS.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c D_HSS.f90

B_HSS.o: B_HSS.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c B_HSS.f90

Ortho_Init.o: Ortho_Init.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c Ortho_Init.f90


