

JELLYFISH_CODE="jellyfish-2.1.4"
HTSLIB_CODE="htslib.github.Aug182015.tar.gz"

#######################
### Required by Trinity
#######################

FASTOOL_CODE=fstrozzi-Fastool-7c3e034f05
PARAFLY_CODE=parafly-code
TRIMMOMATIC_CODE=Trimmomatic-0.32

trinity_essentials: jellyfish scaffold_iworm_contigs_target fastool_target parafly_target trimmomatic_target samtools

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


samtools:
	tar xvf samtools-0.1.19.tar.bz2
	cd samtools-0.1.19 && $(MAKE)
	mv samtools-0.1.19/samtools ./BIN/.

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

fastool_target:
	cd ${FASTOOL_CODE} && $(MAKE)
	ln -sf ${FASTOOL_CODE} fastool


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
	cd collectl && ./build_collectl.sh

clean:
	rm -rf ${JELLYFISH_CODE} ./jellyfish # has tar.gz
	cd ${FASTOOL_CODE} && $(MAKE) clean
	rm -rf ./samtools-0.1.19
	rm -f ./fastool # rm symlink
	rm -rf ./htslib
	cd scaffold_iworm_contigs && $(MAKE) clean
	cd parafly && $(MAKE) clean
	rm -f ./parafly # rm symlink
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd collectl && rm -rf bin doc man
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


