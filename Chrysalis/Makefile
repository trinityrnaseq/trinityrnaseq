all:
	mkdir -p build
	cd build && cmake -DCMAKE_INSTALL_PREFIX="" ../ && make DESTDIR=../ install

clean:
	@echo cleaning
	(cd build && make clean) || :
	rm -rf ./build ./bin
