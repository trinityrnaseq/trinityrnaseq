
OS := $(shell uname)

COLLECTL_CODE="collectl-4.1.0"
PARAFLY_CODE="ParaFly-0.1.0"
TRIMMOMATIC_CODE=Trimmomatic-0.36



trinity_essentials: seqtk_target parafly_target trimmomatic_target

trimmomatic_target:
	ln -sf ${TRIMMOMATIC_CODE} Trimmomatic


seqtk_target:
	tar xvf seqtk-trinity.v0.0.2.tar.gz
	cd seqtk-trinity-0.0.2 && $(MAKE)
	mv seqtk-trinity-0.0.2/seqtk-trinity ./BIN/.

parafly_target:
	tar -zxvf ${PARAFLY_CODE}.tar.gz && \
	cd ${PARAFLY_CODE} && sh ./configure --prefix=`pwd` && $(MAKE) install && \
	cp bin/ParaFly ../BIN/


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
	rm -rf ./seqtk-trinity-0.0.2
	cd scaffold_iworm_contigs && $(MAKE) clean
	cd ${PARAFLY_CODE} && $(MAKE) clean
	rm -f ./parafly # rm symlink
	rm -f ./Trimmomatic # rm symlink
	cd slclust && $(MAKE) clean
	cd COLLECTL && rm -rf ${COLLECTL_CODE} && rm -f collectl
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


