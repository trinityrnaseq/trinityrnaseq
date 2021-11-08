CXX    ?= g++
prefix = ../htslib

ScaffoldIwormContigs:
	$(CXX) $(LDFLAGS) -I$(prefix) -L$(prefix) $(CPPFLAGS) $(CXXFLAGS) ScaffoldIwormContigs.cpp error_checker.cpp -lhts -o scaffold_iworm_contigs

clean:
	rm -f scaffold_iworm_contigs


