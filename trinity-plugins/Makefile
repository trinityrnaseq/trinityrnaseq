
TARGETS=jellyfish rsem transdecoder others


JELLYFISH_CODE="jellyfish-2.1.4"
RSEM_CODE="rsem-1.2.19"
TRANSDECODER_CODE="TransDecoder_r20140704"

RSEM_LIBCURSES=-lcurses

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	RSEM_LIBCURSES += -ltinfo 
endif

all: ${TARGETS}
	@echo "\n\n** Done building plugins **\n\n"

jellyfish:
	tar -zxvf ${JELLYFISH_CODE}.tar.gz && ln -sf ${JELLYFISH_CODE} tmp.jellyfish
	cd ./tmp.jellyfish/ && ./configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` && $(MAKE) LDFLAGS="-lpthread -all-static" AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`"/include"
	mv tmp.jellyfish jellyfish



rsem:
	tar -zxvf ${RSEM_CODE}.tar.gz && ln -sf ${RSEM_CODE} tmp.rsem
	cd ./tmp.rsem && $(MAKE) LIBCURSES="${RSEM_LIBCURSES}"
	mv tmp.rsem rsem

transdecoder:
	tar -zxvf ${TRANSDECODER_CODE}.tar.gz && ln -sf ${TRANSDECODER_CODE} tmp.transdecoder
	cd ./tmp.transdecoder && ${MAKE} simple
	mv ./tmp.transdecoder transdecoder

others:
	cd fastool && $(MAKE)
	cd parafly && ./configure --prefix=`pwd` && $(MAKE) install
	cd slclust && $(MAKE) install
	cd collectl && ./build_collectl.sh

clean:
	rm -rf ./${JELLYFISH_CODE} ./jellyfish
	rm -rf ./${RSEM_CODE} ./rsem
	rm -rf ./${TRANSDECODER_CODE} ./transdecoder
	cd fastool && $(MAKE) clean
	cd parafly && $(MAKE) clean
	cd slclust && $(MAKE) clean
	cd collectl && rm -rf bin doc man
	@echo "\n\n** Done cleaning plugins area **"


###################################################################


