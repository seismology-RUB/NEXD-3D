#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#
#  General definitions
# Check Operating system:
OS := $(shell uname)
#
#  set the SHELL
#   MSB MSB: Actually nothing should be specified here! (or SHELL = /bin/sh if really necessary for any reason...)
SHELL = /bin/tcsh
#
#  Paths
#
srcdir = ./src
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
AR = ar
F95 = mpif90

ifeq ($(notdir $(F95)),g95)
        FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
#-march=opteron
#-march=core2

        LINKFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow -Wunused
        FFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow  -Wunused -O3

# profiling and debugging
#        LINKFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow -Wall -Wextra -pg -fprofile-arcs -ftest-coverage
#        FFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow -Wall -Wextra -pg	-fprofile-arcs -ftest-coverage
# for debugging
#        LINKFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow -g
#        FFLAGS =   -J$(moduledir)  -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunderflow \
	 -TENV:simd_zmask=OFF -TENV:simd_imask=OFF -TENV:simd_omask=OFF -g	

#ifort
#	LINKFLAGS = -module $(moduledir) -heap-arrays 
#	FFLAGS = -module $(moduledir) -heap-arrays 
#	LINKFLAGS = -O3 -module $(moduledir) -heap-arrays -diag-enable sc3
#	FFLAGS = -O3 -module $(moduledir) -heap-arrays -diag-enable sc3
#	FFLAGS = -O3 -module $(moduledir) -heap-arrays -funroll-loops 
#	LINKFLAGS= -O3 -xP -vec-report0 -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe0 -ftz -ftrapuv -check all -traceback -module $(moduledir)
#	FFLAGS= -O3 -xP -vec-report0 -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe0 -ftz -ftrapuv -check all -traceback -module $(moduledir)

#	LINKFLAGS= -O3 -xP -vec-report0 -e95 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe0 -ftz -ftrapuv -check all -traceback -module $(moduledir)
#	FFLAGS= -O3 -xP -vec-report0 -e95 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe0 -ftz -ftrapuv -check all -traceback -module $(moduledir)




#	FFLAGS= /F64000000 -O3 -xP -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe0 -ftz -ftrapuv -check all -traceback -module $(moduledir)

endif

obj_prog = 	program_dg3d.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	timeloopMod.o \
	matrixMod.o \
	waveMod.o \
	dtMod.o \
	stfMod.o \
	imageMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	fileParameterMod.o \
        anelasticMod.o\
        linearSystemMod.o\
	mpiMod.o\
        pmlMod.o\
        rungekuttaMod.o\
        allocateMod.o\
        typeMod.o

#obj_mesher = program_mesher.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	matrixMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o

obj_mesher = program_mesher.o \
	constantsMod.o \
	parameterMod.o \
	tetraNeighborMod.o \
	meshMod.o \
	vtkMod.o \
	nodesMod.o \
	gllMod.o \
	warpfactorMod.o \
	matrixMod.o \
	tetTrafoMod.o \
	rosettaGammaMod.o \
	jacobiMod.o \
	vandermondeMod.o\
	simplexMod.o \
	dmatricesMod.o \
	liftMod.o \
	geometricFactorsMod.o \
	normalsMod.o \
	sourceReceiverMod.o \
	fileParameterMod.o \
        anelasticMod.o\
        linearSystemMod.o\
	logoMod.o\
        allocateMod.o\
        typeMod.o\
        mpiMod.o\
        stfMod.o

obj_solver = program_solver.o \
	mpiincludeMod.o \
	timeloopMod.o \
	constantsMod.o \
	parameterMod.o \
	tetraNeighborMod.o \
	meshMod.o \
	vtkMod.o \
	nodesMod.o \
	gllMod.o \
	warpfactorMod.o \
	matrixMod.o \
	tetTrafoMod.o \
	rosettaGammaMod.o \
	jacobiMod.o \
	vandermondeMod.o\
	simplexMod.o \
	dmatricesMod.o \
	liftMod.o \
	geometricFactorsMod.o \
	normalsMod.o \
	sourceReceiverMod.o \
	mpiMod.o \
	stfMod.o \
	vtkMod.o \
	plotMod.o \
	waveMod.o \
	fileParameterMod.o \
        anelasticMod.o\
        linearSystemMod.o\
	logoMod.o\
        pmlMod.o\
        rungekuttaMod.o\
        allocateMod.o\
        typeMod.o

obj_movie = program_movie.o \
	constantsMod.o \
	parameterMod.o \
	collectMovieMod.o \
	tetraNeighborMod.o \
	meshMod.o \
	vtkMod.o \
	nodesMod.o \
	gllMod.o \
	warpfactorMod.o \
	matrixMod.o \
	tetTrafoMod.o \
	rosettaGammaMod.o \
	jacobiMod.o \
	vandermondeMod.o\
	simplexMod.o \
	dmatricesMod.o \
	liftMod.o \
	geometricFactorsMod.o \
	normalsMod.o \
	sourceReceiverMod.o \
	logoMod.o \
	vtkMod.o \
	fileParameterMod.o \
        anelasticMod.o\
        linearSystemMod.o\
	plotMod.o \
        allocateMod.o\
        typeMod.o\
        mpiMod.o\
        stfMod.o



#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
#vpath %.f90 $(srcdir) ./src/include
#vpath %.f $(srcdir) ./src/include
#vpath %.c $(srcdir) ./src/include

#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = 
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: $(srcdir)/%.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
obstringtest = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Library paths
#
#ifeq ($(OS), Linux)
#        pgplot = -lpgplot -L/usr/lib -lpng -lz -lX11
#endif
#ifeq ($(OS), Darwin)
#        pgplot = -lpgplot -L/usr/local/lib -lpng -lz -lX11
#endif
la = -llapack -lblas
#lib = /usr/lib/libmetis.a /usr/lib/libscotch.a /usr/lib/libscotcherr.a
lib = metis-4.0.3/libmetis.a
#-L/opt/ompi/lib -lmpi
#---------------------------------------------------------
.PHONY: all

#
#  create dependencies on modules 
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	./tools/scripts/makeDepFromUseInclude.py $(srcdir) $(DEPDIRS) > $@
-include make.incdep
#
#       Targets
#
all: required mesher solver movie

clean :
	rm -f $(bindir)/* $(moduledir)/*.mod $(obsdir)/*.o out/* out2/* *.mod tmp.log make.incdep

dg3d: all

mesher: $(obj_mesher)
	$(F95) $(LINKFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

solver: $(obj_solver)
	$(F95) $(LINKFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

movie: $(obj_movie)
	$(F95) $(LINKFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

bin:
	mkdir -p $(bindir)

mod:
	mkdir -p $(moduledir)

obj:
	mkdir -p $(obsdir)

required: mod obj bin
