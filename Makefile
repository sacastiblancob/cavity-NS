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
SRCLECB = $(SRCDIR)/lecbound.f
SRCALLC = $(SRCDIR)/all_csc.f
SRCDIAG = $(SRCDIR)/csc_diag.f
SRCDIAA = $(SRCDIR)/csc_diaga.f
SRCKRON = $(SRCDIR)/csc_kron.f
SRCSUMC = $(SRCDIR)/csc_sum.f
SRCTRAN = $(SRCDIR)/csc_trans.f
SRCPSOR = $(SRCDIR)/csc_presor.f
SRCCMMV = $(SRCDIR)/csc_mmatvec.f
SRCPCSO = $(SRCDIR)/csc_preconSSOR.f
SRCPCIL = $(SRCDIR)/csc_preconilu0.f
SRCSPLU = $(SRCDIR)/csc_solpacklu.f
SRCHOUV = $(SRCDIR)/householderv.f
SRCARHO = $(SRCDIR)/csc_arnoldihouse.f
SRCCCGS = $(SRCDIR)/csc_cg.f
SRCSORS = $(SRCDIR)/csc_sor.f
SRCPCGM = $(SRCDIR)/csc_pcg.f
SRCBICG = $(SRCDIR)/csc_bicgstab.f
SRCGMRE = $(SRCDIR)/csc_gmres.f
SRCDIFF = $(SRCDIR)/diffusion_matrix.f
SRCLAPL = $(SRCDIR)/laplacian_matrix.f
SRCPREG = $(SRCDIR)/point_regularization.f
SRCDIV2 = $(SRCDIR)/diver2d.f
SRCGRA2 = $(SRCDIR)/grad2d.f
SRCPOIS = $(SRCDIR)/poisson.f
SRCDSOL = $(SRCDIR)/diffusion.f
SRCUPDA = $(SRCDIR)/update_and_write.f
SRCKILL = $(SRCDIR)/killemall.f
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
OBJLECB = $(OBJDIR)/lecbound.o
OBJALLC = $(OBJDIR)/all_csc.o
OBJDIAG = $(OBJDIR)/csc_diag.o
OBJDIAA = $(OBJDIR)/csc_diaga.o
OBJKRON = $(OBJDIR)/csc_kron.o
OBJSUMC = $(OBJDIR)/csc_sum.o
OBJTRAN = $(OBJDIR)/csc_trans.o
OBJPSOR = $(OBJDIR)/csc_presor.o
OBJCMMV = $(OBJDIR)/csc_mmatvec.o
OBJPCSO = $(OBJDIR)/csc_preconSSOR.o
OBJPCIL = $(OBJDIR)/csc_preconilu0.o
OBJSPLU = $(OBJDIR)/csc_solpacklu.o
OBJHOUV = $(OBJDIR)/householderv.o
OBJARHO = $(OBJDIR)/csc_arnoldihouse.o
OBJCCGS = $(OBJDIR)/csc_cg.o
OBJSORS = $(OBJDIR)/csc_sor.o
OBJPCGM = $(OBJDIR)/csc_pcg.o
OBJBICG = $(OBJDIR)/csc_bicgstab.o
OBJGMRE = $(OBJDIR)/csc_gmres.o
OBJDIFF = $(OBJDIR)/diffusion_matrix.o
OBJLAPL = $(OBJDIR)/laplacian_matrix.o
OBJPREG = $(OBJDIR)/point_regularization.o
OBJDIV2 = $(OBJDIR)/diverd2.o
OBJGRA2 = $(OBJDIR)/grad2d.o
OBJPOIS = $(OBJDIR)/poisson.o
OBJDSOL = $(OBJDIR)/diffusion.o
OBJUPDA = $(OBJDIR)/update_and_write.o
OBJKILL = $(OBJDIR)/killemall.o
#
OBJECTS = $(OBJDCSC) $(OBJCSCS) $(OBJPHYS) $(OBJNUME) $(OBJWRHE) \
			   	$(OBJPONS) $(OBJLECB) $(OBJALLC) $(OBJDIAG) $(OBJDIAA) \
				 	$(OBJKRON) $(OBJSUMC) $(OBJTRAN) $(OBJPSOR) $(OBJCMMV) \
					$(OBJPCSO) $(OBJPCIL) $(OBJSPLU) $(OBJHOUV) $(OBJARHO) \
				 	$(OBJCCGS) $(OBJSORS) $(OBJPCGM) $(OBJBICG) $(OBJGMRE) \
					$(OBJDIFF) $(OBJLAPL) $(OBJPREG) $(OBJDIV2) \
					$(OBJGRA2) $(OBJPOIS) $(OBJDSOL) $(OBJUPDA) $(OBJKILL) \
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

$(OBJLECB): $(SRCLECB)
	    $(FC) $(CFLAGS) $(SRCLECB) -o $(OBJLECB)

$(OBJALLC): $(SRCALLC)
	    $(FC) $(CFLAGS) $(SRCALLC) -o $(OBJALLC)

$(OBJDIAG): $(SRCDIAG)
	    $(FC) $(CFLAGS) $(SRCDIAG) -o $(OBJDIAG)

$(OBJDIAA): $(SRCDIAA)
	    $(FC) $(CFLAGS) $(SRCDIAA) -o $(OBJDIAA)

$(OBJKRON): $(SRCKRON)
	    $(FC) $(CFLAGS) $(SRCKRON) -o $(OBJKRON)

$(OBJSUMC): $(SRCSUMC)
	    $(FC) $(CFLAGS) $(SRCSUMC) -o $(OBJSUMC)

$(OBJTRAN): $(SRCTRAN)
	    $(FC) $(CFLAGS) $(SRCTRAN) -o $(OBJTRAN)

$(OBJPSOR): $(SRCPSOR)
	    $(FC) $(CFLAGS) $(SRCPSOR) -o $(OBJPSOR)

$(OBJCMMV): $(SRCCMMV)
	    $(FC) $(CFLAGS) $(SRCCMMV) -o $(OBJCMMV)

$(OBJPCSO): $(SRCPCSO)
	    $(FC) $(CFLAGS) $(SRCPCSO) -o $(OBJPCSO)

$(OBJPCIL): $(SRCPCIL)
	    $(FC) $(CFLAGS) $(SRCPCIL) -o $(OBJPCIL)

$(OBJSPLU): $(SRCSPLU)
	    $(FC) $(CFLAGS) $(SRCSPLU) -o $(OBJSPLU)

$(OBJHOUV): $(SRCHOUV)
	    $(FC) $(CFLAGS) $(SRCHOUV) -o $(OBJHOUV)

$(OBJARHO): $(SRCARHO)
	    $(FC) $(CFLAGS) $(SRCARHO) -o $(OBJARHO)

$(OBJCCGS): $(SRCCCGS)
	    $(FC) $(CFLAGS) $(SRCCCGS) -o $(OBJCCGS)

$(OBJSORS): $(SRCSORS)
	    $(FC) $(CFLAGS) $(SRCSORS) -o $(OBJSORS)

$(OBJPCGM): $(SRCPCGM)
	    $(FC) $(CFLAGS) $(SRCPCGM) -o $(OBJPCGM)

$(OBJBICG): $(SRCBICG)
	    $(FC) $(CFLAGS) $(SRCBICG) -o $(OBJBICG)

$(OBJGMRE): $(SRCGMRE)
	    $(FC) $(CFLAGS) $(SRCGMRE) -o $(OBJGMRE)

$(OBJDIFF): $(SRCDIFF)
	    $(FC) $(CFLAGS) $(SRCDIFF) -o $(OBJDIFF)

$(OBJLAPL): $(SRCLAPL)
	    $(FC) $(CFLAGS) $(SRCLAPL) -o $(OBJLAPL)

$(OBJPREG): $(SRCPREG)
	    $(FC) $(CFLAGS) $(SRCPREG) -o $(OBJPREG)

$(OBJDIV2): $(SRCDIV2)
	    $(FC) $(CFLAGS) $(SRCDIV2) -o $(OBJDIV2)

$(OBJGRA2): $(SRCGRA2)
	    $(FC) $(CFLAGS) $(SRCGRA2) -o $(OBJGRA2)

$(OBJPOIS): $(SRCPOIS)
	    $(FC) $(CFLAGS) $(SRCPOIS) -o $(OBJPOIS)

$(OBJDSOL): $(SRCDSOL)
	    $(FC) $(CFLAGS) $(SRCDSOL) -o $(OBJDSOL)

$(OBJUPDA): $(SRCUPDA)
	    $(FC) $(CFLAGS) $(SRCUPDA) -o $(OBJUPDA)

$(OBJKILL): $(SRCKILL)
	    $(FC) $(CFLAGS) $(SRCKILL) -o $(OBJKILL)

$(OBJMAIN): $(MODULES) $(SRCMAIN)
	    $(FC) $(CFLAGS) $(SRCMAIN) -o $(OBJMAIN)

clean:
	    rm -f $(OBJECTS)
			rm -f $(MODULES)
			rm -f $(BINNSFD)
#END



