##
## Makefile
##  
## Made by John Demme
## Login   <teqdruidteqBook>
##
## Started on  Sat Mar 28 14:45:01 2009 John Demme
## Last update Sat Mar 28 14:45:01 2009 John Demme
## 
##############################
# Complete this to make it ! #
##############################
SRC	= storage.cpp pairwisealign.cpp 
TESTS   = tests/Set1.cpp tests/Alignment.cpp
INCL  	= storage.h algo.h msa.h
################
# Optional add #
################
CXXFLAGS = -I. -Wall -g -O3 -fopenmp
LDFLAGS  = -g -fopenmp

#####################
# Macro Definitions #
#####################
OBJS 	= $(SRC:.cpp=.o)
TESTOBJS =$(TESTS:.cpp=.o)
RM 	= /bin/rm -f 	
STRIP	= strip

all:	test speedTest msa

test: test.o $(OBJS) $(TESTOBJS) $(INCL)
	g++ -o $@ test.o $(OBJS) $(TESTOBJS) $(LDFLAGS) -lcppunit

speedTest: speedTest.o $(OBJS) $(INCL)
	g++ -o $@ speedTest.o $(OBJS) $(LDFLAGS)

msa: msa.o $(OBJS) $(INCL)
	g++ -o $@ msa.o $(OBJS) $(LDFLAGS)

clean:
	-$(RM) test *.o $(OBJS) $(TESTOBJS) *~

.cpp.o:
	g++ $(CXXFLAGS) -c $< -o $@

jup:
	rsync --exclude=Sequences -av . jdd2127@compute:stochastic-msa/


################
# Dependencies #
################

$(TESTOBJS): $(INCL)
$(OBJS): $(INCL)