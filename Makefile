#
FC= gfortran
EXECUTEABLE=GeoHBond.out

CFLAGS=-cpp -c -CC -ffree-form -DISPNTLIST=1 
DFLAGS=-ggdb3 -pedantic -Wall -ffpe-trap=overflow -fcheck=all -fbacktrace
OFLAGS=-funroll-loops  -O2  #Turn on Compiler optimization
LDFLAGS=

WD=$(shell pwd)
SRC=$(addprefix $(WD)/,src)
BUILD=$(addprefix $(WD)/,build)
#a set of variables and names required for Program version documentation
COMPILE_TIME=$(shell date +"%a, %Y.%m.%d %H:%M:%S %Z")
MACHINE=$(shell echo $$HOSTNAME)
GIT_COMMIT=$(shell git rev-parse HEAD)
CFLAGS+=-D__PROGRAM_COMPILE_TIME="'$(COMPILE_TIME)'" \
	    -D__MACHINE_NAME="'$(MACHINE)'" -D__GIT_VERSION="'$(GIT_COMMIT)'"

#for GenStorMod define the CompileFlags to build AtomStorMod
CFLAGS+=-D__TYPE="type(Atom)" -D__USE_DECLARATION="use AtomMod, only : Atom" \
	   	-D__MNAME="'AtomStorMod'" -D__MODULE_NAME="AtomStorMod"

#define empty files to indicate whether previous run was a debug or compile run
DBGFLAG=$(addprefix $(BUILD)/,.dbg)
CPLFLAG=$(addprefix $(BUILD)/,.cpl)

INCDIR=$(addprefix $(WD)/,../repo/include)
LIBDIR=$(addprefix $(WD)/,../repo/lib)
LIBS=-lfftw3 -lbase

T2AMODS= 
T2BMODS= AtomListMod.F Molecule.F ReadPSFMod.F DetMolsMod.F 

STRCMODS= GenGrids.F
LOCMODS= ProvideGeoMod.F DistancesMod.F Main.F
		 
		 
		 
		 
		 
		 
INCLUDE= $(addprefix $(SRC)/, Main_subroutines.F)
SOURCES=$(addprefix $(SRC)/, $(T2AMODS) $(T2BMODS) $(STRCMODS) $(LOCMODS))
OBJECTS=$(addprefix $(BUILD)/,$(T2AMODS:.f95=.o) $(T2BMODS:.F=.o) \
	    $(STRCMODS:.F=.o) $(LOCMODS:.F=.o))


.PHONY: directories clean all
all: CFLAGS:=$(CFLAGS) $(OFLAGS)
all: LDFLAGS:=$(OFLAGS)
all: directories $(CPLFLAG) $(SOURCES) $(EXECUTEABLE) 
	@echo ----------------------
	@echo  ***Build complete***	
	@echo ----------------------

$(EXECUTEABLE): $(OBJECTS)
	$(FC) $(LDFLAGS) $^ -L$(LIBDIR) $(LIBS) -o $(BUILD)/$@
	ln --force $(BUILD)/$(EXECUTEABLE) $(WD)

#bit of a hack. could be done better.
$(BUILD)/Main.o: $(SRC)/Main.F $(INCLUDE)
	$(FC) $(CFLAGS)  -I$(INCDIR) $< -J $(BUILD) -o $@
$(BUILD)/%.o: $(SRC)/%.f95
	$(FC) $(CFLAGS)  -I$(INCDIR) $< -J $(BUILD) -o $@ 
$(BUILD)/%.o: $(SRC)/%.F
	$(FC) $(CFLAGS)  -I$(INCDIR) $< -J $(BUILD) -o $@

directories:
	@if [ ! \( -e $(BUILD) \) ];then mkdir $(BUILD);fi;
clean:
	@if [ \( -e $(BUILD) \) ];then rm -r build;fi;
	@echo ***Cleanup complete***
debug: CFLAGS:=$(CFLAGS) $(DFLAGS)
debug: LDFLAGS:=$(DFLAGS)
debug: directories $(DBGFLAG) $(SOURCES) $(EXECUTEABLE)
	@echo '--------------------------'
	@echo '  ***Build complete***	 '
	@echo '--------------------------'

$(DBGFLAG):
	rm -f $(CPLFLAG) $(OBJECTS)
	touch $(DBGFLAG)
$(CPLFLAG):
	rm -f $(DBGFLAG) $(OBJECTS)
	touch $(CPLFLAG)
