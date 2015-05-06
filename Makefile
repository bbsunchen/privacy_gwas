CXX=g++
CXXFLAGS=-g -c -I /gpfs/home/cxs1031/backup/bin/include/
CXXLINK=-static -L/gpfs/home/cxs1031/backup/bin/lib

SRCSIMPLE=beta.cpp functions.cpp roc.cpp privacy.cpp gaussian.cpp
OBJSIMPLE=$(SRCSIMPLE:%.cpp=%.o)

all:privacy

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $< 
	
privacy:$(OBJSIMPLE)
	$(CXX) $(CXXLINK) -o $@ $^ -lgsl -lgslcblas

clean:
	rm  -f *.o privacy
