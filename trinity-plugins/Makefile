
JELLYFISH_CODE="jellyfish-2.1.4"
HTSLIB_CODE="htslib.github.Aug182015.tar.gz"
COLLECTL_CODE="collectl-4.1.0"
FASTOOL_CODE=fstrozzi-Fastool-7c3e034f05
PARAFLY_CODE=parafly-code
TRIMMOMATIC_CODE=Trimmomatic-0.32

# set ltinfo use for samtools if on linux and it exists on the system.
ifeq ($(shell uname -s),"Linux")
	LTINFO_SEARCH=$(shell ldconfig -p |grep tinfo)
endif

ifneq (,$(findstring tinfo,${LTINFO_SEARCH}))
	LTINFO="LIBPATH=-ltinfo"
endif



trinity_essentials: jellyfish scaffold_iworm_contigs_target seqtk_target parafly_target trimmomatic_target samtools

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


samtools:
	tar xvf samtools-0.1.19.tar.bz2
	cd samtools-0.1.19 && $(MAKE) $(LTINFO)
	mv samtools-0.1.19/samtools ./BIN/.

seqtk_target:
	tar xvf seqtk-trinity.v0.0.2.tar.gz
	cd seqtk-trinity-0.0.2 && $(MAKE)
	mv seqtk-trinity-0.0.2/seqtk-trinity ./BIN/.

jellyfish:
	tar -zxvf ${JELLYFISH_CODE}.tar.gz && ln -sf ${JELLYFISH_CODE} tmp.jellyfish
	cd ./tmp.jellyfish/ && ./configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` && $(MAKE) LDFLAGS="-lpthread -all-static" AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`"/include"
	mv tmp.jellyfish jellyfish

htslib_target:
#	tar xjvf htslib-1.2.1.tar.bz2
#	cd htslib-1.2.1 && ./configure && $(MAKE)
	tar xvf ${HTSLIB_CODE} && cd htslib && $(MAKE)

scaffold_iworm_contigs_target: htslib_target
	cd scaffold_iworm_contigs && $(MAKE)


parafly_target:
	cd ${PARAFLY_CODE} && sh ./configure --prefix=`pwd` && $(MAKE) install
	ln -sf ${PARAFLY_CODE} parafly


##################################
### Needed for downstream analyses
##################################

plugins: slclust_target collectl_target
	@echo "\n\n** Done building plugins **\n\n"


slclust_target:
	cd slclust && $(MAKE) install

collectl_target:
	cd COLLECTL && tar xvf ${COLLECTL_CODE}.src.tar.gz && ln -sf ${COLLECTL_CODE} collectl

clean:
	rm -rf ${JELLYFISH_CODE} ./jellyfish # has tar.gz
	rm -rf ./samtools-0.1.19
	rm -rf ./seqtk-trinity-0.0.2
	rm -rf ./htslib
	cd scaffold_iworm_contigs && $(MAKE) clean
	cd parafly-code && $(MAKE) clean
	rm -f ./parafly # rm symlink
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd COLLECTL && rm -rf ${COLLECTL_CODE} && rm -f collectl
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


