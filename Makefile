CXX=g++
CXXFLAGS=-g -c -I /gpfs/home/cxs1031/backup/bin/include/
CXXLINK=-static -L/gpfs/home/cxs1031/backup/bin/lib

SRCFILES=data.cpp  functions.cpp  roc.cpp  securegenome.cpp
OBJFILES=$(SRCFILES:%.cpp=%.o)

SRCSIMPLE=data.cpp functions.cpp roc.cpp simplefunction.cpp
OBJSIMPLE=$(SRCSIMPLE:%.cpp=%.o)

all:securegenome

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $< 

securegenome:$(OBJFILES)
	$(CXX) $(CXXLINK) -o $@ $^ -lgsl -lgslcblas
	
simplefunction:$(OBJSIMPLE)
	$(CXX) $(CXXLINK) -o $@ $^ -lgsl -lgslcblas

clean:
	rm  -f *.o securegenome simplefunction
