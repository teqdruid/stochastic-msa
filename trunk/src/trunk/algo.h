/*
** algo.h
** 
*/

#ifndef   	ALGO_H_
# define   	ALGO_H_

#include "storage.h"
#include <vector>

template <class M, class A, class B> //A and B must be a Sequences
    double alignmentScore(M& m, A& a, B& b, SiteInformation* si = NULL);

template<class M, class T>
vector<T>* nwAlignment(M& S,
		       ImmutableSequence<T>& a,
		       ImmutableSequence<T>& b);

#endif 	    /* !ALGO_H_ */
