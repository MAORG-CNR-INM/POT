LIB   =parallel
LNAME =lib$(LIB).a
LOBJS =parallel.o 
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

parallel.o: parallel.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c parallel.f


