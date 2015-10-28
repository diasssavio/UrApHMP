# System architecture
SYSTEM     = x86-64_linux

# Static library format for Cplex
LIBFORMAT  = static_pic

# Source code folder
SRC	= src
INCLUDE = include

# Machine hostname
MACHINE = $(shell hostname)

# Library type(STATIC or DYNAMIC)
MERGE = DYNAMIC

##### Folders
# Temp folders
TMP_ILS = ./tmp/ILS
TMP_STATIC = ./tmp/lib/static
# Perm folders
DAT_DOXYFILE = ./dat/doxyfile
DAT_INSTANCES = ./dat/instances
DAT_LP_MODELS = ./dat/lp_models
DAT_RESULTS = ./dat/results


# Cplex directory
CPLEXDIR	  = /opt/ibm/ILOG/CPLEX_Studio1261/cplex

# Concert directory
CONCERTDIR	  = /opt/ibm/ILOG/CPLEX_Studio1261/concert

# Compiler
CCC = g++-4.8

# Compilation parameters (Add afterward: --coverage -pg -ftree-vectorize -mfpmath=sse -march=native)
CCOPT = -std=gnu++0x -O3 -ftree-vectorize -mfpmath=sse -march=native -march=native -flto -g -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD

# Cplex static libraries directory
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# Concert static libraries directory
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# Include libraries identifiers
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread

# Cplex header's directory
CPLEXINCDIR   = $(CPLEXDIR)/include

# Concert header's directory
CONCERTINCDIR = $(CONCERTDIR)/include

# Header's include path
# CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
CCFLAGS = $(CCOPT)

# Executable name
CPP_EX = ILS-UrApHMP

# Compiling
all:
	mkdir -p $(TMP_ILS)
	mkdir -p $(TMP_STATIC)
	mkdir -p $(DAT_DOXYFILE)
	mkdir -p $(DAT_INSTANCES)
	mkdir -p $(DAT_LP_MODELS)
	mkdir -p $(DAT_RESULTS)
	make -j8 $(CPP_EX);

# Executing
execute: $(CPP_EX)
	./$(CPP_EX)

# Cleaning
clean:
	# /bin/rm -rf $(CPP_EX)
	/bin/rm -rf ./tmp
	/bin/rm -rf ./dat


########################## GENERATING OBJECT's ######################################################

# CONFIGURATION - INSTANCES
$(TMP_ILS)/UrApHMP.o: $(SRC)/UrApHMP.cpp $(INCLUDE)/UrApHMP.h
	$(CCC) -c $(CCFLAGS) $(SRC)/UrApHMP.cpp -o $(TMP_ILS)/UrApHMP.o

# STRUCTURE - SOLUTION
$(TMP_ILS)/solution.o: $(SRC)/solution.cpp $(INCLUDE)/solution.h
	$(CCC) -c $(CCFLAGS) $(SRC)/solution.cpp -o $(TMP_ILS)/solution.o

# STRUCTURE - TIMER
$(TMP_ILS)/FWChrono.o: $(SRC)/FWChrono.cpp $(INCLUDE)/FWChrono.h
	$(CCC) -c $(CCFLAGS) $(SRC)/FWChrono.cpp -o $(TMP_ILS)/FWChrono.o

# STRUCTURE - RANDOM GEN
$(TMP_ILS)/mt19937ar.o: $(SRC)/mt19937ar.c $(INCLUDE)/mt19937ar.h
	$(CCC) -c $(CCFLAGS) $(SRC)/mt19937ar.c -o $(TMP_ILS)/mt19937ar.o

# ILS
$(TMP_ILS)/ils.o: $(SRC)/ils.cpp $(INCLUDE)/ils.h
	$(CCC) -c $(CCFLAGS) $(SRC)/ils.cpp -o $(TMP_ILS)/ils.o

# ILS
$(TMP_ILS)/grasp.o: $(SRC)/grasp.cpp $(INCLUDE)/grasp.h
	$(CCC) -c $(CCFLAGS) $(SRC)/grasp.cpp -o $(TMP_ILS)/grasp.o

# MAIN
$(TMP_ILS)/main.o: $(SRC)/main.cpp
	$(CCC) -c $(CCFLAGS) $(SRC)/main.cpp -o $(TMP_ILS)/main.o

########################## OBJECT's LIBRARIES #######################################################
# CONFIGURATION
$(TMP_ILS)/Configuration.o:  $(TMP_ILS)/UrApHMP.o
	gcc -Wl,-r  $(TMP_ILS)/UrApHMP.o -o $(TMP_ILS)/Configuration.o -nostdlib

# STRUCTURE & TIMER
$(TMP_ILS)/Structure.o: $(TMP_ILS)/solution.o $(TMP_ILS)/FWChrono.o $(TMP_ILS)/mt19937ar.o
	gcc -Wl,-r $(TMP_ILS)/solution.o $(TMP_ILS)/FWChrono.o $(TMP_ILS)/mt19937ar.o -o $(TMP_ILS)/Structure.o -nostdlib

# ILS
$(TMP_ILS)/ILS.o: $(TMP_ILS)/ils.o
	gcc -Wl,-r $(TMP_ILS)/ils.o -o $(TMP_ILS)/ILS.o -nostdlib

# GRASP
$(TMP_ILS)/GRASP.o: $(TMP_ILS)/grasp.o
	gcc -Wl,-r $(TMP_ILS)/grasp.o -o $(TMP_ILS)/GRASP.o -nostdlib

########################## LINKANDO TUDO ########################################################

$(CPP_EX): $(TMP_ILS)/Configuration.o $(TMP_ILS)/Structure.o $(TMP_ILS)/ILS.o $(TMP_ILS)/GRASP.o $(TMP_ILS)/main.o
	$(CCC)  $(CCFLAGS) $(TMP_ILS)/Configuration.o $(TMP_ILS)/Structure.o $(TMP_ILS)/ILS.o $(TMP_ILS)/GRASP.o $(TMP_ILS)/main.o -L$(TMP_STATIC) -o $(CPP_EX)
#endif
