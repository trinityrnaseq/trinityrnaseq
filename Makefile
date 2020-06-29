###################################################################
#
# The default compiler is GNU gcc/g++.
# Run
#  make TRINITY_COMPILER=intel
# to build Inchworm and Chrysalis with the Intel compiler.
#


export TRINITY_HOME = $(PWD)

ifeq ($(TRINITY_COMPILER),intel)
 INCHWORM_CONFIGURE_FLAGS = CXX=icpc CXXFLAGS="-fast"
 CHRYSALIS_MAKE_FLAGS = COMPILER=icpc
else
 override TRINITY_COMPILER=gnu
endif


all: inchworm_target chrysalis_target trinity_essentials
	sh ./util/support_scripts/trinity_install_tests.sh


no_bamsifter: inchworm_target chrysalis_target
	cd trinity-plugins && $(MAKE) no_bamsifter

install:
	util/support_scripts/trinity_installer.py

inchworm_target:
	@echo Using $(TRINITY_COMPILER) compiler for Inchworm and Chrysalis
	cd Inchworm && $(MAKE)

chrysalis_target:
	cd Chrysalis && $(MAKE)


trinity_essentials:
	cd trinity-plugins && $(MAKE) trinity_essentials


plugins:
	cd trinity-plugins && $(MAKE) plugins
	sh ./util/support_scripts/plugin_install_tests.sh

test:
	@echo
	@echo "Checking for Trinity essentials (built from 'make all'):"
	sh ./util/support_scripts/trinity_install_tests.sh
	@echo
	@echo "Checking for plugins (built from 'make plugins'):"
	sh ./util/support_scripts/plugin_install_tests.sh
	@echo "Run 'make test_trinity' if you want to test Trinity execution on a small data set"

clean:
	cd Inchworm && $(MAKE) clean
	cd Chrysalis && $(MAKE) clean 
	cd trinity-plugins && $(MAKE) clean 
	cd sample_data/ && $(MAKE) clean


test_trinity:
	cd sample_data/test_Trinity_Assembly && make test


# note 'test_all': ** this is for a more advanced installation including devel features **

test_all:
	cd sample_data/ && make test_all
	./__pull_trinity_ext_sample_data.sh
	cd trinity_ext_sample_data/ && make test

test_clean:
	cd sample_data/ && make clean
	cd trinity_ext_sample_data/ && make clean

###################################################################


