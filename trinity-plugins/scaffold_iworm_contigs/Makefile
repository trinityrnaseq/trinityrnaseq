CXX    = g++
prefix = ../htslib

ScaffoldIwormContigs:
	$(CXX) -I$(prefix) -L$(prefix) ScaffoldIwormContigs.cpp error_checker.cpp -lhts -lz -o scaffold_iworm_contigs

clean:
	rm scaffold_iworm_contigs


