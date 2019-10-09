#  Output file name
# EXC = abox7.out
  EXC = exec.out
#  OUT_FILES = *.txt
# 
# loading routines to compile. 
#  OPTIONS here.   
#

FOBJS = sam.f90 gandalf.f90 eru.f90 istari.f90 peregrin.f90 strider.f90 main.f90

FC = gfortran
FFLAGS = -fdefault-real-8 -fdefault-double-8


#  Compile

$(EXC): 	$(FOBJS)
		$(FC) $(FFLAGS) $(FOBJS) $(USE_INCL) $(LIBS) -o $@

$(USE_INCL):	$(INCLUDES)

.f.o:
		$(FC) $(FFLAGS) -c $*.f

clean:
		rm -f core $(EXC) $(USE_INCL) $(OUT_FILES) *.o *.mod



     

