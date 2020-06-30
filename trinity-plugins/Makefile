
OS := $(shell uname)

CXX = g++
CC = gcc
COLLECTL_CODE="collectl-4.1.0"
PARAFLY_CODE="ParaFly-0.1.0"
TRIMMOMATIC_CODE=Trimmomatic-0.36



trinity_essentials: seqtk_target parafly_target trimmomatic_target bamsifter_target

no_bamsifter: seqtk_target parafly_target trimmomatic_target


trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


seqtk_target:
	cd seqtk-trinity && $(MAKE) CXX=$(CXX) CC=$(CC)
	mv seqtk-trinity/seqtk-trinity ./BIN/.

parafly_target:
	tar -zxvf ${PARAFLY_CODE}.tar.gz && \
	cd ${PARAFLY_CODE} && sh ./configure --prefix=`pwd` CXX=$(CXX) CC=$(CC) CFLAGS="-fopenmp" CXXFLAGS="-fopenmp" && $(MAKE) install && \
	cp bin/ParaFly ../BIN/


bamsifter_target:
	cd bamsifter && make

##################################
### Needed for downstream analyses
##################################

plugins: slclust_target collectl_target
	@echo "\n\n** Done building plugins **\n\n"


slclust_target:
	cd slclust && $(MAKE) CXX=$(CXX) CC=$(CC) install

collectl_target:
	cd COLLECTL && tar xvf ${COLLECTL_CODE}.src.tar.gz && ln -sf ${COLLECTL_CODE} collectl

clean:
	cd seqtk-trinity && $(MAKE) clean
	cd scaffold_iworm_contigs && $(MAKE) clean
	rm -rf ${PARAFLY_CODE}
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd COLLECTL && rm -rf ${COLLECTL_CODE} && rm -f collectl
	cd htslib && make clean
	cd bamsifter && make clean
	rm -f ./htslib.ok
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


