
OS := $(shell uname)

COLLECTL_CODE="collectl-4.1.0"
PARAFLY_CODE="ParaFly-0.1.0"
TRIMMOMATIC_CODE=Trimmomatic-0.36



trinity_essentials: seqtk_target parafly_target trimmomatic_target

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


seqtk_target:
	cd seqtk-trinity-0.0.2 && $(MAKE) CXX=g++ CC=gcc
	mv seqtk-trinity-0.0.2/seqtk-trinity ./BIN/.

parafly_target:
	tar -zxvf ${PARAFLY_CODE}.tar.gz && \
	cd ${PARAFLY_CODE} && sh ./configure --prefix=`pwd` CXX=g++ CC=gcc CFLAGS="-fopenmp" CXXFLAGS="-fopenmp" && $(MAKE) install && \
	cp bin/ParaFly ../BIN/


##################################
### Needed for downstream analyses
##################################

plugins: slclust_target collectl_target
	@echo "\n\n** Done building plugins **\n\n"


slclust_target:
	cd slclust && $(MAKE) CXX=g++ CC=gcc install

collectl_target:
	cd COLLECTL && tar xvf ${COLLECTL_CODE}.src.tar.gz && ln -sf ${COLLECTL_CODE} collectl

clean:
	cd seqtk-trinity-0.0.2 && $(MAKE) clean
	cd scaffold_iworm_contigs && $(MAKE) clean
	rm -rf ${PARAFLY_CODE}
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd COLLECTL && rm -rf ${COLLECTL_CODE} && rm -f collectl
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


