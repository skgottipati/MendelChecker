# GNU GCC compiler is required
CC=g++
# FLAGS for optimization
CFLAGS=-std=c++11 -O3 -Wall
SRCDIR=src
SOURCES= $(SRCDIR)/mcc_series_optimization_sexbased_clean.cpp $(SRCDIR)/fileread3.cpp $(SRCDIR)/computeLikelihood.cpp $(SRCDIR)/cart.cpp $(SRCDIR)/genotypeLikelihood.cpp $(SRCDIR)/GenotypeInfo.cpp optparse/OptionParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MendelChecker
all: 
	$(CC) $(CFLAGS) $(SOURCES) -o $@
