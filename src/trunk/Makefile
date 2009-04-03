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
INCL  	= storage.h algo.h
################
# Optional add #
################
CXXFLAGS = -I. -Wall -g # -O2
LDFLAGS  = -g

#####################
# Macro Definitions #
#####################
OBJS 	= $(SRC:.cpp=.o)
TESTOBJS =$(TESTS:.cpp=.o)
RM 	= /bin/rm -f 	
STRIP	= strip

all:	test

test: test.o $(OBJS) $(TESTOBJS) $(INCL)
	g++ -o $@ test.o $(OBJS) $(TESTOBJS) -lcppunit

clean:
	-$(RM) test $(OBJS) $(TESTOBJS) *~

.cpp.o:
	g++ $(CXXFLAGS) -c $< -o $@


################
# Dependencies #
################

$(TESTOBJS): $(INCL)
$(OBJS): $(INCL)