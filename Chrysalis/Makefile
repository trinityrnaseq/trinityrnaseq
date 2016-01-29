###################################################################################
#
# OVERVIEW OF THE SPINES MAKEFILE
#
# I. General comments
#
# This Makefile works by running the Spines program MakeDepend to determine
# dependencies, which are placed in a file Makefile_auto, every time you do a 
# make.	 Consequently, the first time you run this Makefile, it will appear to be 
# slow, because it compiles MakeDepend and its prerequisites, silently.
# 
# To utilize MakeDepend fully, you should adhere to the convention that for every
# .cc file which is not a main program, there is a corresponding .h file, which
# is included by any .cc file which accesses the corresponding routines.
#
# II. Major targets
# 
# make: makes all executables
# make d (where d is a subdirectory of Spines): make all executables in d
# make assembler: make all executables needed to run Spines
#
# III. Supported platforms
#
# compiler: gcc 4.1.1
#
# hardware: ia64 (Itanium)   - we use it all the time
#           i686 (Pentium)   - we use it a lot of the time
#           x86_64 (Opteron) - other people tell us that it works
#           alpha            - we used to use it, probably unusable now
#           macosx           - not fully functional
#
# Note. These are all compiled with -O3, but with for ia64, we also use
# -fno-rename-registers, to compensate for a gcc bug 
# (see testing/IsWeirdness.cc, testing/IsWeirdness.run, gcc bug report #14815).
# We also compile with -fno-strict-aliasing, without which errors can result,
# and these errors are extremely difficult to trace.
#
# IV. Common options
#
# FORCE_DEBUG=yes: compile in debugging mode
# OPTIM=none: compile without optimization
# UNSUPPORTED=yes: override check for supported platform
# ALL=yes: compile files, even if they have the line // MakeDepend: archived
# -k: keep going in the event of error (a make option)
# -jn: use up to n parallel processors to do the make (a make option)
#
###################################################################################


##############################################################################
# DEFINE GENERAL VARIABLES
##############################################################################

DATE = $(shell date)
OS_NAME = $(shell uname -s)
NODE_NAME = $(shell uname -n)
OS_RELEASE = $(shell uname -r)
OS_MACHINE = $(shell uname -m)
OS_VENDOR = $(VENDOR)
USER = $(shell whoami)
RELEASE = 3.0

ifeq ($(OS_MACHINE),Power Macintosh)
  OS_MACHINE = macosx
endif

# this ought to be done more elegantly, but this ugly hack suffices for now
GPP_VERSION := $(shell g++ -dumpversion)
GPP_MAJOR_VERSION := $(shell g++ -dumpversion | tr '.' '\n' | head -1)

##############################################################################
# CHECK FOR SUPPORTED PLATFORM
##############################################################################

ifneq ($(UNSUPPORTED),yes)
     ifneq ($(SUBRUN),1)

	  ifneq ($(GPP_VERSION),4.1.1)
               $(warning ----------------------------------------------------------)
               $(warning using g++ version $(GPP_VERSION) )
               $(warning --------------------------------------------------------)
	  endif

          ACCEPT_HARDWARE := no
	  ifeq ($(OS_MACHINE),x86_64)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(OS_MACHINE),ia64)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(OS_MACHINE),i686)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(OS_MACHINE),macosx)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(OS_MACHINE),alpha)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(OS_MACHINE),i86pc)
	       ACCEPT_HARDWARE := yes
          endif
	  ifeq ($(ACCEPT_HARDWARE),no)
               $(warning ----------------------------------------------------------)
               $(warning Platform $(OS_MACHINE) is UNSUPPORTED)
               $(error --------------------------------------------------------)
	  endif

     endif
else
#     $(warning ------------------------------------------------------------)
#     $(warning WARNING --- Unsupported platform --- proceed with caution)
#     $(warning ------------------------------------------------------------)
endif


##############################################################################
## DIRECTORY SETUP:
##############################################################################

SRC = .
OBJ = obj

ifdef SPINES_BIN_DIR
     BIN = $(SPINES_BIN_DIR)
else
     BIN = .
endif

LIB = .

# Make cross-compiling in one directory easier.
ifdef CROSS
  CROSS_SUFFIX = $(OS_MACHINE)_$(OS_VENDOR)
  OBJ = obj_$(CROSS_SUFFIX)
  BIN = bin_$(CROSS_SUFFIX)
  # If DEBUGGING==yes, put the object and binary files in a separate dir.
  # The DEBUGGING option also has other effects later in the Makefile.
  ifeq ($(DEBUGGING),yes)
    OBJ = obj_$(CROSS_SUFFIX)_debug
    BIN = bin_$(CROSS_SUFFIX)_debug
  else
    # If PROFILE==yes, put the object and binary files in a separate dir.
    # The PROFILE option also has other effects later in the Makefile.
    ifeq ($(PROFILE),yes)
      OBJ = obj_$(CROSS_SUFFIX)_prof
      BIN = bin_$(CROSS_SUFFIX)_prof
    endif
  endif
endif

vpath %.cc $(SRC)
vpath %.c $(SRC)
vpath %.h $(SRC)
vpath %.hh $(SRC)
vpath %.o $(OBJ) .
vpath % $(BIN)

## This part should be removed before distribution vvv ##

##############################################################################
## RECURSIVE MAKE DEPEND
##############################################################################

ifndef SUBRUN
  ifdef MAKEDEPEND_STATUS
    $(warning Finding dependencies...)
  endif
  # result := $(shell bash -c "if [ ! { $(MAKE) CXX=$(CXX) CROSS=$(CROSS) BIN=$(BIN) OBJ=$(OBJ) ALL=$(ALL) FULL=$(FULL) OPTIM=$(OPTIM) SUBRUN=1 depend } ]; then echo MakeDepend FAILED fi")
  result := $(shell bash -c " $(MAKE) CXX=$(CXX) CROSS=$(CROSS) BIN=$(BIN) OBJ=$(OBJ) ALL=$(ALL) FULL=$(FULL) OPTIM=$(OPTIM) SUBRUN=1 depend ")
  ifneq ($(findstring FAILED,$(result)),)
    $(error MakeDepend failed)
  endif
  # Note this really checks for any WARNING output, but currently 
  # the only way that text can be output is by checkLock --checkSetLock.
  ifneq ($(findstring WARNING,$(result)),)
    $(warning checkSetLock failed. Running processes are not protected against overwrite from other computers.)
  endif
  ifdef MAKEDEPEND_STATUS
    $(warning Done.)
  endif
endif

## This part should be removed before distribution ^^^ ##

##############################################################################
## COMPILE-TIME CONFIGURATION:
##############################################################################

## These can be overridden from the `make' command line.  For example,
##
##   % make FORCE_DEBUG=yes
##
## builds debugging executables.
##
##  In most cases, one must run `make clean' before changing the defaults
## below to ensure that executables are rebuilt correctly.  Generally, the
## compiler-specific makefiles sort out these options and turn on the right
## compiler flags to get the desired behavior.

# Which compiler to use; can be `g++', `cxx', `icpc', or `sunCC'.
COMPILER	= g++
#COMPILER	= cxx
#COMPILER	= sunCC
#COMPILER	= icpc

# QUIET=yes causes the compiler to be less verbose about warnings.
QUIET		= no

# DEBUG=yes forces out-of-date executables to be built with symbolic
# debugging information.  Any other value is ignored.
DEBUG		= yes

# OPTIM=<compiler optimization> forces user-defined optimizations to be used,
# rather than the compiler-specific default.  
# If OPTIM=fast, then less optimization is done to speed compilation. 
# If OPTIM=none, then no optimization is used.	
# If unset, defaults are used.
OPTIM		=

# PROFILE=yes builds sampled-profiling (i.e. gprof) executables. Any other
# value is ignored.
PROFILE		= no

# FORCE_DEBUG=yes clears NDEBUG (i.e. turns on macros from <assert.h>).	 Any
# other value is ignored.
FORCE_DEBUG	= no

# ASSERT=FAST turns on fast(er) assertion checking.  Any other value is
# ignored.
ASSERT		= SLOW

# OPEN_MP=yes enables OpenMP directives.  Any other value is ignored.
OPEN_MP		= yes

# Enables Posix multithreading.	 This is currently only supported for CXX
# compilations, and partially for GCC v3.  This is enabled by default, since
# it has no effect in Gcc compilations.	 Any other value is ignored.
PTHREADS	= yes

# GCOV=yes builds test-coverage analyzing (gcov) executables. Suggested use
# is with OPTIM=none
GCOV		= no

# TRACK_MEMORY=yes turns on memory tracking (see MemTracker.h) in all modules
# which include that header.
TRACK_MEMORY    = no

# DEBUGGING=yes turns off optimization and turns on FORCE_DEBUG.
# DEBUGGING also has an effect earlier in the Makefile.
ifeq ($(DEBUGGING),yes)
  OPTIM = none
  FORCE_DEBUG = yes
endif

# flags to pass to $(AR). POSIX default is "-rv".
ARFLAGS = -qc 

#Add here additional options to be turned on by MakeDepend's
#directive mechanism. For example:
TEST_DIRECTIVE 	   =  -O0 -ffast-math -funroll-loops
OPT1               =  -O1


##############################################################################
## BASE SYSTEM SETUP:
##############################################################################

## These macros are used in Makefile and Makefile_auto to compile SPINES.
## The compiler-specific makefile augments these macros before they are used
## used to build the system.  These are here to ensure that the macros are
## instantiated, and as a reminder that the macros exist.

# Location of Xerces XML parser library:
XERCES_INC	=
XERCES_LIB	=

# Base definitions:
SYS_DEFS	= \
		  -DMAKE_DATE='"$(DATE)"' \
		  -DMAKE_OS_RELEASE='"$(OS_RELEASE)"' \
		  -DMAKE_RELEASE='"$(RELEASE)"' \
		  -DNEW_MAKEFILE

# Warning control:
SYS_WARN	=

# Optimization control:
SYS_OPT		=

# Language control (e.g. template instantiations):
SYS_LANG	=

# Where to find system-specific includes:
ifeq ($(COMPILER),g++)
  ifeq ($(GPP_MAJOR_VERSION),2)
       SYS_INCS	= -I. -I- $(XERCES_INC)
  endif
  ifeq ($(GPP_MAJOR_VERSION),3)
       SYS_INCS	= -I. -I- $(XERCES_INC)
  endif
  ifeq ($(GPP_MAJOR_VERSION),4)
       SYS_INCS	= -iquote . $(XERCES_INC)
  endif
  ifeq ($(GPP_MAJOR_VERSION),5)
       SYS_INCS        = -iquote . $(XERCES_INC)
  endif
endif

# Linking control (e.g. to link templates):
SYS_LINK	=

# Required libraries:
SYS_LIBS	= -lm -pthread

# OpenMP macros:
ifeq ($(OPEN_MP),yes)
  OMP_FLAGS		=
  OMP_INCS		=
  OMP_DEFS		=

  OMP_OPTIONS		= \
			$(OMP_FLAGS) \
			$(OMP_INCS) \
			$(OMP_DEFS)

  OMP_LINK		=
  OMP_LIBS		=
endif

# Posix threads macros:
ifeq ($(PTHREADS),yes)
  PTHREAD_FLAGS		=
  PTHREAD_INCS		=
  PTHREAD_DEFS		= -DHAVE_PTHREADS

  PTHREAD_OPTIONS	= \
			$(PTHREAD_FLAGS) \
			$(PTHREAD_INCS) \
			$(PTHREAD_DEFS)

  PTHREAD_LINK	=
  PTHREAD_LIBS	=
endif

ifeq ($(OPTIM),fast)
  FAST		= yes
endif


##############################################################################
## COMPILATION:
##############################################################################

# Bring in the compiler-specific makefile:

ifeq ($(COMPILER),g++)
include Makefile_g++
endif

ifeq ($(COMPILER),cxx)
include Makefile_cxx
endif

ifeq ($(COMPILER),sunCC)
include Makefile_sunCC
endif

ifeq ($(COMPILER),icpc)
include Makefile_icpc
endif

# Setup final compilation options:

CPP_OPTIONS	= \
		$(SYS_WARN) \
		$(SYS_OPT) \
		$(SYS_DEBUG) \
		$(SYS_DEFS) \
		$(SYS_LANG) \
		$(SYS_INCS) \
		$(OMP_OPTIONS) \
		$(PTHREAD_OPTIONS)

LINK_OPTIONS	= \
		$(SYS_DEBUG) \
		$(SYS_LINK) \
		$(OMP_LINK) \
		$(PTHREAD_LINK)

LINK_LIBS	= \
		$(SYS_LIBS) \
		$(OMP_LIBS) \
		$(PTHREAD_LIBS) \

CPPO		= $(LINK_OPTIONS)
CPPC		= $(CPP_OPTIONS) $(SYS_MEMTRACK)

##############################################################################
## COMPILATION TARGET SETUP:
##############################################################################

# Standard targets:

include Makefile_targets

# Makefile_auto: depend

ALL_ARG =
ifeq ($(ALL),yes)
     ALL_ARG = -A
endif

FULL_RECURSION_ARG =
ifeq ($(FULL),yes)
     FULL_RECURSION_ARG = -f
endif

depend: MakeDepend checkLock
	$(BIN)/checkLock --checkSetLock
	@$(BIN)/MakeDepend $(ALL_ARG) $(FULL_RECURSION_ARG) $(SRC)

# Import automatic dependencies here, but not if doing a make depend; the
# `-' before the include prevents a warning from being issued.
# N.B.: `-' also prevents an _error_ from being issued if Makefile_auto
#	couldn't be found or built!
ifeq ($(findstring depend,$(MAKECMDGOALS)),)
-include Makefile_auto
endif

clean_execs:
	for file in $(EXECUTABLES); do rm -f $(BIN)/$$file; done

clean:
	for file in $(EXECUTABLES); do rm -f $(BIN)/$$file; done
	rm -f MakeDepend $(BIN)/MakeDepend contigs.out my.permanent.log.file \
             core a.out Makefile.bak bsubin BasevectorTables.h $(BIN)/checkLock 
	find $(OBJ) -name '*.o' -exec rm {} \;
	rm -rf cxx_repository
	rm -f lib_*_temp.a

execs_size:
	du -sk $(ASSEMBLY_EXECS)
	du -sk $(ASSEMBLY_EXECS) | awk '{ sum += $$1 } END { print sum }'

# Non-standard targets:
MakeDepend: MakeDepend.cc
	@ mkdir -p $(BIN)
	$(CPLUSPLUS) $(CPPO) $(CPPC) $(SYS_LANG) $(MAKEDEPEND_OPTS) -o $(BIN)/MakeDepend $(SRC)/MakeDepend.cc

checkLock: $(SRC)/util/checkLock.c
	$(CC) $(SRC)/util/checkLock.c -o $(BIN)/checkLock 

Spines/random/Random.o: Spines/random/Random.cc
	@ mkdir -p $(OBJ)/${@D}
	cp $(SRC)/Spines/random/Random.cc $(SRC)/Spines/random/Random.c
	$(CC) -c $(SRC)/Spines/random/Random.c -w -Drandom=randomx -Dsrandom=srandomx -o $(OBJ)/Spines/random/Random.o
	rm $(SRC)/Spines/random/Random.c

MemTracker.o: MemTracker.cc
	@ mkdir -p $(OBJ)
	$(CPLUSPLUS) $(CPP_OPTIONS) -o $(OBJ)/MemTracker.o -c $(SRC)/MemTracker.cc

procbuf_local.o: procbuf_local.cc procbuf_local.h

# How to compile C++ files:
.cc.o:
	@ mkdir -p $(OBJ)/${@D}
	$(CPLUSPLUS) $(CPPC) -c $(SRC)/$*.cc -o $(OBJ)/$*.o

##############################################################################
## ADDITIONAL TARGETS:
##############################################################################

## Prepping SPINES for download:
## For example:
##
##    % make download DOWNLOAD_DIR=<destination>
##
## This will create three tarred, gzipped directories: the binaries, the
## source, and a valid data directory with sample project.

ifndef $(DOWNLOAD_DIR)
  DOWNLOAD_DIR	= for_download
endif

BIN_DL_DIR	= $(DOWNLOAD_DIR)/Spines_bin
SRC_DL_DIR	= $(DOWNLOAD_DIR)/Spines_src
DATA_DL_DIR	= $(DOWNLOAD_DIR)/Spines_data

download: bin_download src_download data_download

bin_download: $(EXECUTABLES) ViewAlignments.cc Time doc/Manual.html LICENSE.txt lib$(XERCES_LIB).so LICENSE.xerces.txt
	rm -rf $(BIN_DL_DIR)
	mkdir -p $(BIN_DL_DIR)
	cp $(EXECUTABLES) Time $(BIN_DL_DIR)
	cp doc/Manual.html doc/ReleaseNotes.html $(BIN_DL_DIR)
	cp LICENSE.txt $(BIN_DL_DIR)
	cp ViewAlignments.cc $(BIN_DL_DIR)
	cp lib$(XERCES_LIB).so LICENSE.xerces.txt $(BIN_DL_DIR)
	( cd $(DOWNLOAD_DIR); \
	  tar -cvf Spines_bin.tar Spines_bin; \
	  rm -rf Spines_bin; \
	  gzip -vf Spines_bin.tar )

src_download:
	rm -rf $(SRC_DL_DIR)
	mkdir -p $(SRC_DL_DIR)
	find . \( -name "*.c" -o \
		  -name "*.cc" -o \
		  -name "*.cpp" -o \
		  -name "*.h" -o \
		  -name "*.hh" -o \
		  -name "*.hpp" \) -print | cpio -pd $(SRC_DL_DIR)
	cp Makefile $(SRC_DL_DIR)
	cp LICENSE.txt $(SRC_DL_DIR)
	cp LICENSE.xerces.txt $(SRC_DL_DIR)/xerces_include
	( cd $(DOWNLOAD_DIR); \
	  tar -cvf Spines_src.tar --exclude=CVS Spines_src; \
	  rm -rf Spines_src; \
	  gzip -vf Spines_src.tar )

data_download:
	( cd $(DOWNLOAD_DIR); \
	  rm -rf Spines_data; \
	  cvs checkout Spines_data; \
	  tar -cvf Spines_data.tar --exclude=CVS Spines_data; \
	  rm -rf Spines_data; \
	  gzip -vf Spines_data.tar )

kecho:
	@echo $(command) 

gcov_clean:
	@touch tmp.bb	
	rm -f *.bb *.bbg *.da *.gcov

diff: CvsListFiles
	@./kcomp `CvsListFiles LEVELS=1`
