/*
** algo.h
** 
*/

#ifndef   	ALGO_H_
# define   	ALGO_H_

#include "storage.h"
#include <vector>

template <class M, class A, class B> //A and B must be a Sequences
    double alignmentScore(M& m, A& a, B& b);

template<class M, class A, class B>
const char*** getAlignment(M& m, A& a, B& b);

template<class A>
void freeAlignment(const char***, A& a);

template<class A, class B>
void reconstructAlignment(ostream& os, A& a, B& b, const char*** directions);

template<class M, class T>
vector<T>* nwAlignment(M& S,
		       ImmutableSequence<T>& a,
		       ImmutableSequence<T>& b);

#endif 	    /* !ALGO_H_ */
