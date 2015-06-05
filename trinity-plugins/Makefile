

JELLYFISH_CODE="jellyfish-2.1.4"
RSEM_CODE="rsem-1.2.19"

RSEM_LIBCURSES=-lcurses

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	RSEM_LIBCURSES += -ltinfo 
endif


#######################
### Required by Trinity
#######################

FASTOOL_CODE=fstrozzi-Fastool-7c3e034f05
PARAFLY_CODE=parafly-code
TRIMMOMATIC_CODE=Trimmomatic-0.32

trinity_essentials: jellyfish scaffold_iworm_contigs fastool_target parafly_target trimmomatic_target

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


jellyfish:
	tar -zxvf ${JELLYFISH_CODE}.tar.gz && ln -sf ${JELLYFISH_CODE} tmp.jellyfish
	cd ./tmp.jellyfish/ && ./configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` && $(MAKE) LDFLAGS="-lpthread -all-static" AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`"/include"
	mv tmp.jellyfish jellyfish

htslib:
	tar xjvf htslib-1.2.1.tar.bz2
	cd htslib-1.2.1 && ./configure && $(MAKE)

scaffold_iworm_contigs: htslib
	cd scaffold_iworm_contigs && $(MAKE)

fastool_target:
	cd ${FASTOOL_CODE} && $(MAKE)
	ln -sf ${FASTOOL_CODE} fastool


parafly_target:
	cd ${PARAFLY_CODE} && sh ./configure --prefix=`pwd` && $(MAKE) install
	ln -sf ${PARAFLY_CODE} parafly


##################################
### Needed for downstream analyses
##################################

plugins: slclust_target collectl_target rsem
	@echo "\n\n** Done building plugins **\n\n"


rsem:
	tar -zxvf ${RSEM_CODE}.tar.gz
	cd ${RSEM_CODE} && $(MAKE) LIBCURSES="${RSEM_LIBCURSES}"
	ln -sf ${RSEM_CODE} rsem


# no longer bundling transdecoder
#transdecoder_target:
#	tar -zxvf ${TRANSDECODER_CODE}.tar.gz 
#	cd ${TRANSDECODER_CODE} && ${MAKE} simple
#	ln -sf ${TRANSDECODER_CODE} transdecoder

slclust_target:
	cd slclust && $(MAKE) install

collectl_target:
	cd collectl && ./build_collectl.sh

clean:
	rm -rf ./${JELLYFISH_CODE} ./jellyfish # has tar.gz
	rm -rf ./${RSEM_CODE} ./rsem # has tar.gz
	rm -rf ./${TRANSDECODER_CODE} ./transdecoder  # has tar.gz
	cd ${FASTOOL_CODE} && $(MAKE) clean
	rm -f ./fastool # rm symlink
	cd parafly && $(MAKE) clean
	rm -f ./parafly # rm symlink
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd collectl && rm -rf bin doc man
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


