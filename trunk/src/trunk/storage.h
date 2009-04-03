/*
** Storage.h
**
*/

#ifndef   	STORAGE_H_
#define   	STORAGE_H_

#include <malloc.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>
#include <exception>

using namespace std;

class MsaException : public std::exception {
	public:
		string msg;
		MsaException(const char* msg)
		{
			this->msg = msg;
		}

		virtual ~MsaException() throw() {}
		virtual const char* what()
		{
			return msg.c_str();
		}
};

typedef enum
{
    A = 0,
    C,
    G,
    T
} GeneticSymbols;

static char gs2cLookup[] = {'A', 'C', 'G', 'T'};

__inline__ char toChar(GeneticSymbols gs)
{
    return gs2cLookup[(int)gs];
}

__inline__ GeneticSymbols fromChar(char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
    	return A;
    case 'C':
    case 'c':
    	return C;
    case 'G':
    case 'g':
    	return G;
    case 'T':
    case 't':
    	return T;
    default:
    	throw new MsaException("Invalid GS character");
    }
}


template <class T>
void parse(T** data, size_t* size, istream&);

template <class T>
class ImmutableSequence
{
    T* seq;
    size_t seqSize;

public:
    ImmutableSequence(istream s)
    {
    	parse(&seq, &seqSize, s);
    }

    ImmutableSequence(string seqstr)
    {
    	seqSize = seqstr.length();
    	stringstream ss(stringstream::in | stringstream::out);
    	ss << seqSize << seqstr;
    	parse(&seq, &seqSize, ss);
    }

    ~ImmutableSequence()
    {
    	free(this->seq);
    }

    T at(size_t idx)
    {
    	assert(idx <= seqSize);
    	return seq[idx];
    }

    size_t length()
    {
    	return seqSize;
    }

    string toString()
    {
    	char buf[seqSize+1];
    	for (size_t i=0; i<seqSize; i++)
    	{
    		buf[i] = toChar(seq[i]);
    	}
    	buf[seqSize] = 0;
    	return string(buf);
    }

    bool operator==(ImmutableSequence& o)
    {
    	if (o.seqSize != seqSize)
    		return false;
    	return (memcmp(o.seq, seq, seqSize) == 0);
    }
};

template <class T>
class MutableSequence
{
public:

    MutableSequence(ImmutableSequence<T>& seq);
    T at(size_t idx);
};

template <class T, int S>
class ScoringMatrix {
public:
    double matrix[S*S];
    double gap;
    float alpha;
    float beta;

    ScoringMatrix() {}

    /*ScoringMatrix(double alpha, double gap)
    {
    	double diag = 1.0 - alpha * S;
    	for (size_t i = 0; i<(S*S); i++)
    	{
    		matrix[i] = alpha;
    	}

    	for (size_t i=0; i<S; i++)
    	{
    		matrix[i*S + i] = diag;
    	}
    }*/

    ScoringMatrix(float alpha, float beta)
    {
    	this->alpha = alpha;
    	this->beta = beta;
    }

    double score(T a, T b)
    {
    	if(a == b)
    	{
    		return 1;
    	}
    	else
    	{
    		return 0;
    	}
    	//return matrix[a*S + b];
    }

    double affineScore(int gapSize)
    {
    	if(gapSize < 1)
    	{
    		return 0;
    	}
    	return -1 * (alpha + beta*(gapSize-1));
    }
};

typedef ScoringMatrix<GeneticSymbols, 4> GenScores;
typedef ImmutableSequence<GeneticSymbols> GISeq;

#endif 	    /* !STORAGE_H_ */
