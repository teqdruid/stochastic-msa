##
## Makefile
##  
## Made by John Demme
## Login   <teqdruidteqBook>
##
## 
##############################
# Complete this to make it ! #
##############################
SRC	= storage.cpp pairwisealign.cpp msa.cpp
TESTS   = tests/Set1.cpp tests/Alignment.cpp
PROGS   = test.spp speedTest.cpp msa_main.cpp pairalign.cpp
INCL  	= storage.h algo.h msa.h util.h
################
# Optional add #
################
CXXFLAGS = -I. -Wall -g -fopenmp -O3 -msse3 -finline --save-temps
LDFLAGS  = -fopenmp -g -O3

#####################
# Macro Definitions #
#####################
OBJS 	= $(SRC:.cpp=.o)
PROGOBJS = $(PROGS:.cpp=.o)
TESTOBJS =$(TESTS:.cpp=.o)
RM 	= /bin/rm -f 	
STRIP	= strip

all:	test speedTest msa pairalign bali_score sstest

bali_score:
	gcc -o bali_score bali_score.c -O3

sstest: sstest.cpp
	g++ -o sstest sstest.cpp -I. -O3

test: test.o $(OBJS) $(TESTOBJS) $(INCL)
	g++ -o $@ test.cpp $(OBJS) $(TESTOBJS) $(CXXFLAGS) $(LDFLAGS) -lcppunit

speedTest: speedTest.o $(OBJS) $(INCL)
	g++ -o $@ speedTest.cpp $(OBJS) $(LDFLAGS)  $(CXXFLAGS)

msa: msa_main.o $(OBJS) $(INCL)
	g++ -o $@ msa_main.cpp $(OBJS) $(LDFLAGS) $(CXXFLAGS)

pairalign: pairalign.o $(OBJS) $(INCL)
	g++ -o $@ pairalign.cpp $(OBJS) $(LDFLAGS) $(CXXFLAGS)

clean:
	-$(RM) test *.o $(OBJS) $(TESTOBJS) bali_score *~

.cpp.o:
	g++ $(CXXFLAGS) -c $< -o $@

jup:
	rsync --exclude=Sequences -av . jdd2127@compute:stochastic-msa/
	rsync --exclude=Sequences -av jdd2127@compute:stochastic-msa/ComputeOutputs/* ComputeOutputs/


################
# Dependencies #
################

$(TESTOBJS): $(INCL)
$(OBJS): $(INCL)
