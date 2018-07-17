OBJS = myTree.o analysis.o
EXE = analysis

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -m32  -Wall -g -std=c++11 $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -m32 -o $(EXE) $(OBJS)  $(ROOTFLAGS) $(ROOTLIBS)
%.o: %.cxx
	$(COMPILE)  $< 

%.o: %.C
	$(COMPILE)  $< 


