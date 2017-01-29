
OS := $(shell uname)

JELLYFISH_CODE="jellyfish-2.1.4"
COLLECTL_CODE="collectl-4.1.0"
PARAFLY_CODE=parafly-code
TRIMMOMATIC_CODE=Trimmomatic-0.32



trinity_essentials: jellyfish seqtk_target parafly_target trimmomatic_target samtools_target

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


samtools_target:
	tar xvf samtools-1.3.1.tar.bz2
	cd samtools-1.3.1 && ./configure --prefix=`pwd`
ifeq ($(OS),Darwin)
	sed -i.bak s/-rdynamic//g samtools-1.3.1/Makefile
endif
	cd samtools-1.3.1 && $(MAKE)
	mv samtools-1.3.1/samtools ./BIN/.

seqtk_target:
	tar xvf seqtk-trinity.v0.0.2.tar.gz
	cd seqtk-trinity-0.0.2 && $(MAKE)
	mv seqtk-trinity-0.0.2/seqtk-trinity ./BIN/.

jellyfish:
	tar -zxvf ${JELLYFISH_CODE}.tar.gz && ln -sf ${JELLYFISH_CODE} tmp.jellyfish
	cd ./tmp.jellyfish/ && ./configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` && $(MAKE) LDFLAGS="-lpthread -all-static" AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`"/include"
	mv tmp.jellyfish jellyfish


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
	rm -rf ./seqtk-trinity-0.0.2
	rm -rf ./samtools-1.3.1/
	cd scaffold_iworm_contigs && $(MAKE) clean
	cd parafly-code && $(MAKE) clean
	rm -f ./parafly # rm symlink
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd COLLECTL && rm -rf ${COLLECTL_CODE} && rm -f collectl
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


