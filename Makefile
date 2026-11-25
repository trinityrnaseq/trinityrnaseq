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


# Optional: rebuild Trinity's Java programs (Butterfly.jar, ExitTester.jar) from source.
# The release already includes these jars, so 'make' only needs a Java runtime; 'make jars'
# needs a JDK (javac and jar, version 8 or later). Override JAVAC=/path/to/javac JAR=/path/to/jar if needed.
JAVAC ?= javac
JAR   ?= jar
JAVAC_RELEASE := $(shell $(JAVAC) --release 8 -version >/dev/null 2>&1 && echo "--release 8" || echo "-source 8 -target 8")

jars: butterfly_jar exittester_jar

butterfly_jar:
	cd Butterfly/Butterfly && $(MAKE) JAVAC="$(JAVAC)" JAR="$(JAR)"
	cp Butterfly/Butterfly/build/Butterfly.jar Butterfly/Butterfly.jar

exittester_jar:
	rm -rf util/support_scripts/build_ExitTester
	mkdir -p util/support_scripts/build_ExitTester
	$(JAVAC) -nowarn $(JAVAC_RELEASE) -d util/support_scripts/build_ExitTester util/support_scripts/ExitTester.java
	$(JAR) cfe util/support_scripts/ExitTester.jar ExitTester -C util/support_scripts/build_ExitTester ExitTester.class
	rm -rf util/support_scripts/build_ExitTester

.PHONY: jars butterfly_jar exittester_jar

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
	cd sample_data/test_Trinity_Assembly && $(MAKE) test


# note 'test_all': ** this is for a more advanced installation including devel features **

test_all:
	cd sample_data/ && $(MAKE) test_all
	./__pull_trinity_ext_sample_data.sh
	cd trinity_ext_sample_data/ && $(MAKE) test

test_clean:
	cd sample_data/ && $(MAKE) clean
	cd trinity_ext_sample_data/ && $(MAKE) clean

###################################################################


