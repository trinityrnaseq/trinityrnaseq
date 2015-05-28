CXX    = g++
prefix = ../htslib-1.2.1

ScaffoldIwormContigs:
	$(CXX) -I$(prefix) -L$(prefix) ScaffoldIwormContigs.cpp error_checker.cpp -lhts -lz -o scaffold_iworm_contigs

clean:
	rm scaffold_iworm_contigs


