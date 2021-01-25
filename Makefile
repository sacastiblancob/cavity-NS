#
# NAVIER-STOKES FINITE-DIFFERENCES SOLVER MAKE FILE
#
# ######################################################################
#
#brief      MAKEFILE
#
#history    Sergio Castiblanco
#+          13/01/2021
########################################################################
#
# SET COMPILER OPTIONS
#
FC = gfortran
#
# COMPILER FLAGS
#
FCFLAGS = -m64 -O2
#
# MOD FLAG FOR .mod FILES LOCATION
#
MODFLAG = -J
#
# MAIN DIRECTORY
#
NSDIR = $(CURDIR)
#
# DIRECTORIES FOR SOURCE, OBJECTS, MODULES AND BINARIES
#
SRCDIR = $(NSDIR)/src
OBJDIR = $(NSDIR)/obj
MODDIR = $(NSDIR)/mod
BINDIR = $(NSDIR)/bin
#
# COMPILER AND LINKER FLAGS - GNU FORTRAN
#
CFLAGS = -c -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)
LFLAGS = -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)
#
# FLAGS FOR MODULE DIRECTORY
#
CFLAGS += $(MODFLAG) $(MODDIR)
LFLAGS += $(MODFLAG) $(MODDIR)
#########################################################################
#
# SOURCE FILES
#
SRCMAIN = $(SRCDIR)/home_ns.f
SRCNUME = $(SRCDIR)/declarations_numerical.f
SRCPHYS = $(SRCDIR)/declarations_physic.f
SRCCSCS = $(SRCDIR)/csc_storage.f
SRCDCSC = $(SRCDIR)/declarations_csc.f
SRCWRHE = $(SRCDIR)/write_headers.f
SRCPONS = $(SRCDIR)/point_ns.f
SRCALLC = $(SRCDIR)/all_csc.f
SRCDIAG = $(SRCDIR)/csc_diag.f
SRCKRON = $(SRCDIR)/csc_kron.f
SRCSUMC = $(SRCDIR)/csc_sum.f
SRCTRAN = $(SRCDIR)/csc_trans.f
SRCDIFF = $(SRCDIR)/diffusion_matrix.f
SRCLAPL = $(SRCDIR)/laplacian_matrix.f
SRCPREG = $(SRCDIR)/point_regularization.f
# SRCPRHS = $(SRCDIR)/point_rhs.f
SRCCMMV = $(SRCDIR)/csc_mmatvec.f
SRCCCGS = $(SRCDIR)/csc_cg.f
# SRCUPDA = $(SRCDIR)/update_and_write.f
#
# OBJECT FILES
#
OBJMAIN = $(OBJDIR)/home_ns.o
OBJNUME = $(OBJDIR)/declarations_numerical.o
OBJPHYS = $(OBJDIR)/declarations_physic.o
OBJCSCS = $(OBJDIR)/csc_storage.o
OBJDCSC = $(OBJDIR)/declarations_csc.o
OBJWRHE = $(OBJDIR)/write_headers.o
OBJPONS = $(OBJDIR)/point_ns.o
OBJALLC = $(OBJDIR)/all_csc.o
OBJDIAG = $(OBJDIR)/csc_diag.o
OBJKRON = $(OBJDIR)/csc_kron.o
OBJSUMC = $(OBJDIR)/csc_sum.o
OBJTRAN = $(OBJDIR)/csc_trans.o
OBJDIFF = $(OBJDIR)/diffusion_matrix.o
OBJLAPL = $(OBJDIR)/laplacian_matrix.o
OBJPREG = $(OBJDIR)/point_regularization.o
# OBJPRHS = $(OBJDIR)/point_rhs.o
OBJCMMV = $(OBJDIR)/csc_mmatvec.o
OBJCCGS = $(OBJDIR)/csc_cg.o
# OBJUPDA = $(OBJDIR)/update_and_write.o
#
OBJECTS = $(OBJDCSC) $(OBJCSCS) $(OBJPHYS) $(OBJNUME) $(OBJWRHE) \
			   	$(OBJPONS) $(OBJALLC) $(OBJDIAG) $(OBJKRON) $(OBJSUMC) \
				 	$(OBJTRAN) $(OBJDIFF) $(OBJLAPL) $(OBJPREG) $(OBJCMMV) \
					$(OBJCCGS) \
				 	$(OBJMAIN)
#
#
# MODULE FILES
#
MODNUME = $(MODDIR)/declarations_numerical.mod
MODPHYS = $(MODDIR)/declarations_physic.mod
MODDCSC = $(MODDIR)/declarations_csc.mod
MODCSCS = $(MODDIR)/csc_storage.mod
#
MODULES = $(MODDCSC) $(MODCSCS) $(MODPHYS) $(MODNUME)
#
#
# BINARY FILES
#
BINNSFD = $(BINDIR)/ns_df.out
#########################################################################
#
# !!!COMPILING!!!
#
$(BINNSFD): $(OBJECTS)
	    $(FC) $(LFLAGS) $(OBJECTS) -o $(BINNSFD)

$(OBJDCSC): $(SRCDCSC)
	    $(FC) $(CFLAGS) $(SRCDCSC) -o $(OBJDCSC)

$(OBJCSCS): $(SRCCSCS)
	    $(FC) $(CFLAGS) $(SRCCSCS) -o $(OBJCSCS)

$(OBJPHYS): $(SRCPHYS)
	    $(FC) $(CFLAGS) $(SRCPHYS) -o $(OBJPHYS)

$(OBJNUME): $(SRCNUME)
	    $(FC) $(CFLAGS) $(SRCNUME) -o $(OBJNUME)

$(OBJWRHE): $(SRCWRHE)
	    $(FC) $(CFLAGS) $(SRCWRHE) -o $(OBJWRHE)

$(OBJPONS): $(SRCPONS)
	    $(FC) $(CFLAGS) $(SRCPONS) -o $(OBJPONS)

$(OBJALLC): $(SRCALLC)
	    $(FC) $(CFLAGS) $(SRCALLC) -o $(OBJALLC)

$(OBJDIAG): $(SRCDIAG)
	    $(FC) $(CFLAGS) $(SRCDIAG) -o $(OBJDIAG)

$(OBJKRON): $(SRCKRON)
	    $(FC) $(CFLAGS) $(SRCKRON) -o $(OBJKRON)

$(OBJSUMC): $(SRCSUMC)
	    $(FC) $(CFLAGS) $(SRCSUMC) -o $(OBJSUMC)

$(OBJTRAN): $(SRCTRAN)
	    $(FC) $(CFLAGS) $(SRCTRAN) -o $(OBJTRAN)

$(OBJDIFF): $(SRCDIFF)
	    $(FC) $(CFLAGS) $(SRCDIFF) -o $(OBJDIFF)

$(OBJLAPL): $(SRCLAPL)
	    $(FC) $(CFLAGS) $(SRCLAPL) -o $(OBJLAPL)

$(OBJPREG): $(SRCPREG)
	    $(FC) $(CFLAGS) $(SRCPREG) -o $(OBJPREG)

#$(OBJPRHS): $(SRCPRHS)
#	    $(FC) $(CFLAGS) $(SRCPRHS) -o $(OBJPRHS)

$(OBJCMMV): $(SRCCMMV)
	    $(FC) $(CFLAGS) $(SRCCMMV) -o $(OBJCMMV)

$(OBJCCGS): $(SRCCCGS)
	    $(FC) $(CFLAGS) $(SRCCCGS) -o $(OBJCCGS)

#$(OBJUPDA): $(SRCUPDA)
#	    $(FC) $(CFLAGS) $(SRCUPDA) -o $(OBJUPDA)

$(OBJMAIN): $(MODULES) $(SRCMAIN)
	    $(FC) $(CFLAGS) $(SRCMAIN) -o $(OBJMAIN)

clean:
	    rm -f $(OBJECTS)
			rm -f $(MODULES)
			rm -f $(BINNSFD)
#END



