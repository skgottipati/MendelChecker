# GNU GCC compiler is required
CC=g++
# FLAGS for optimization
CFLAGS=-std=c++11 -O3 -Wall
SOURCES= mcc_series_optimization_sexbased_clean.cpp fileread3.cpp computeLikelihood.cpp cart.cpp genotypeLikelihood.cpp GenotypeInfo.cpp optparse/OptionParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MendelChecker
all: 
	$(CC) $(CFLAGS) $(SOURCES) -o $@
