/*
** algo.h
** 
*/

#ifndef   	ALGO_H_
# define   	ALGO_H_

#include "storage.h"

template <class M, class A, class B> //A and B must be a Sequences
    double alignmentScore(M& m, A& a, B& b);

#endif 	    /* !ALGO_H_ */
