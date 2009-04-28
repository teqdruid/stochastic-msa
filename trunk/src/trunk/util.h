/*
** util.h
** 
** Made by John Demme
** Login   <teqdruid@teqMini>
** 
** Started on  Tue Apr 14 14:09:09 2009 John Demme
** Last update Tue Apr 14 14:13:18 2009 John Demme
*/

#ifndef   	UTIL_H_
# define   	UTIL_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <msa.h>
#include <algo.h>

__inline__ ImmutableSequence<GeneticSymbols>* longRndSeq(size_t size) {
    GeneticSymbols* buffer = (GeneticSymbols*)malloc(size*sizeof(GeneticSymbols));
    
    for (size_t i=0; i<size; i++) {
	int r = rand() % 4;
	switch (r) {
	case 0:
	    buffer[i] = A;
	    break;
	case 1:
	    buffer[i] = C;
	    break;
	case 2:
	    buffer[i] = G;
	    break;
	case 3:
	    buffer[i] = T;
	    break;
	}
    }

    return new ImmutableSequence<GeneticSymbols>(buffer, size);
}

__inline__ double get_runtime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return ((double)t.tv_sec) + (((double)t.tv_usec) / 1e6);
}

class Timer {
public:
    double start, end;
    string name;
    Timer(string name) {
	this->name = name;
	start = get_runtime();
    }

    Timer() {
	name = "";
	start = get_runtime();
    }

    ~Timer() {
	end = get_runtime();

	printf("----%15s Timer (s):\t%lf\n", name.c_str(), (end-start));
    }
};

template<class A>
size_t stochasticSelect(A* arr, size_t len) {
    A total = 0;
    A cumArr[len];
    for (size_t i=0; i<len; ++i) {
	total += arr[i];
	cumArr[i] = total;
    }

    A sel = random() % ((long long)total);
    for (size_t i=0; i<len; ++i) {
	if (cumArr[i] > sel)
	    return i;
    }

    return len - 1;
}

#endif 	    /* !UTIL_H_ */
