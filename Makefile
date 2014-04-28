# GNU GCC compiler is required
CC=g++
# FLAGS for optimization
UNAME := $(shell uname)
CFLAGS=-std=c++11 -O3 -Wall
ifeq ($UNAME, Linux)
    CFLAGS=-std=c++0x -O3 -Wall
endif
ifeq ($UNAME, Darwin)
    CFLAGS=-std=c++11 -O3 -Wall -Wno-unused-private-field
endif
ifneq (,$(findstring MINGW,$(UNAME)))
    CFLAGS=-std=c++0x -O3 -Wall
endif
ifneq (,$(findstring CYGWIN,$(UNAME)))
    CFLAGS=-std=c++11 -O3 -Wall
endif
#-Wno-unused-variable
#-Wno-unused-private-field
SRCDIR=src
SOURCES= $(SRCDIR)/mendel_checker.cpp $(SRCDIR)/split.cpp $(SRCDIR)/genoped.cpp $(SRCDIR)/verify_snps.cpp $(SRCDIR)/fileread.cpp $(SRCDIR)/computeLikelihood.cpp $(SRCDIR)/cart.cpp $(SRCDIR)/genotypeLikelihood.cpp $(SRCDIR)/GenotypeInfo.cpp optparse/OptionParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MendelChecker
$(EXECUTABLE): 
	$(CC) $(CFLAGS) $(SOURCES) -o $@
