FC = gfortran
#FC = g95

#FFLAGS = -O3 -fopenmp
FFLAGS = -mcmodel=medium -O3
#FFLAGS = -O0 -ggdb

#LFLAGS = -O3 -fopenmp
LFLAGS = -mcmodel=large -O3		
#LFLAGS = -O0 -ggdb

#LIBS = -lgomp /usr/lib64/liblapack.a /usr/lib64/libblas.a
LIBS = /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libblas.a

#DEBUG = -fsanitize=address
DEBUG =

OBJECTS = interface.o eigensg.o susbicub.o bcucof.o bcuint.o

MODULES = interface.mod 

DATA = eigenrsgn02.dat saida.dat

.PHONY: clean

eigenrsgn02.dat: eigensg.exe
	./eigensg.exe > saida.dat
#	gdb eigensg.exe

eigensg.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o eigensg.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) eigensg.exe

help:
	@echo "Valid targets:"
	@echo "  eigensg.exe"
	@echo "  eigensg.o"
	@echo "  clean:  removes .o and .exe files"
