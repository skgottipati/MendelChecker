# GNU GCC compiler is required
CC=g++
# FLAGS for optimization
CFLAGS=-std=c++11 -O3 -Wall -Wno-unused-private-field
SRCDIR=src
SOURCES= $(SRCDIR)/mendel_checker.cpp $(SRCDIR)/split.cpp $(SRCDIR)/genoped.cpp $(SRCDIR)/verify_snps.cpp $(SRCDIR)/fileread.cpp $(SRCDIR)/computeLikelihood.cpp $(SRCDIR)/cart.cpp $(SRCDIR)/genotypeLikelihood.cpp $(SRCDIR)/GenotypeInfo.cpp optparse/OptionParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MendelChecker
$(EXECUTABLE): 
	$(CC) $(CFLAGS) $(SOURCES) -o $@
