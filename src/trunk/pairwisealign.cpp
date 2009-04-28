#include "algo.h"
#include "storage.h"
#include <iostream>
#include <algorithm>
#include <cmath>

template<class T>
class Matrix {
    T* mat;
    size_t xs, ys;

public:
    Matrix(size_t xs, size_t ys, T init): xs(xs), ys(ys) {
	mat = new T[xs*ys];
	for (size_t i=0; i<xs*ys; ++i) {
	    mat[i] = init;
	}
    }

    ~Matrix() {
	delete[] mat;
    }

    __inline__ T& operator()(size_t x, size_t y) {
	return mat[x + y*xs];
    }
};

/* 
 * Much of this code copied and adapted from the Wikipedia entry:
 *  http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
 *
 * b is the sequence, a is the profile
 */
template<class M, class A, class B>
double alignmentScore(M& S, A& a, B& b, SiteInformation* si = NULL)
{
    vector<double> Fc(b.length()+1), Fb(b.length()+1);
    vector<size_t> hGapsB(b.length()+1), hGapsC(b.length()+1);

    //Backtracking directions
    Matrix<char>   dirct(a.length()+1, b.length()+1, 'S');

    size_t vGaps = 0;
    hGapsB[0] = 0;

    for (size_t j=0; j<=b.length(); ++j) {
	Fb[j] = S.vAffineScore(j);
    }

    for (size_t i = 1; i <= a.length(); ++i) {
	Fc[0] = S.hAffineScore(i);
	hGapsC[0] = i; 

	for (size_t j = 1; j <= b.length(); ++j) {
	    double diag = Fb[j-1] + S(a[i-1], b[j-1]);
	    double horz = Fb[j]   + S.hAffineScore(hGapsB[j] + 1);
	    double vert = Fc[j-1] + S.vAffineScore(vGaps + 1);

	    if (diag > horz && diag > vert) {
		Fc[j] = diag;
		dirct(i, j) = 'D';
		hGapsC[j] = 0;
		vGaps = 0;
	    } else if (horz > vert) {
		Fc[j] = horz;
		hGapsC[j] = hGapsB[j] + 1;
		vGaps = 0;
		dirct(i, j) = 'H';
	    } else {
		Fc[j] = vert;
		vGaps++;
		hGapsC[j] = 0;
		dirct(i, j) = 'V';
	    }
	}

	Fc.swap(Fb);
	hGapsC.swap(hGapsB);
    }

    if (si != NULL) {

	size_t delMul = round(S.vAlpha/S.hAlpha);
	if (delMul == 0)
	    delMul = 1;

	size_t i = a.length();
	size_t j = b.length();
	
	while (i > 0 && j > 0)
	{
	    switch (dirct(i, j)) {
	    case 'D':
		//aln->push_back(b[j-1]);
		if (a[i-1] != b[j-1])
		    si->subst[i-1]++;
		--i;
		--j;
		break;
	    case 'H':
		si->dels[i-1] += delMul;
		--i;
		break;
	    case 'V':
		si->ins[i-1]++;
		--j;
		break;
	    default:
		assert(false);
	    }
	}
    
	while (i > 0)
	{
	    si->dels[i-1] += delMul;
	    --i;
	}

	while (j > 0)
	{
	    si->ins[0]++;
	    --j;
	}
    }

    return Fb[b.length()];
}


/* 
 * Much of this code copied and adapted from the Wikipedia entry:
 *  http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
 *
 * b is the sequence, a is the profile
 */
template<class M, class T>
vector<T>* nwAlignment(M& S,
		       ImmutableSequence<T>& a,
		       ImmutableSequence<T>& b) {

    vector<T>* aln = new vector<T>();
    aln->reserve(a.length());

    Matrix<double> F(a.length()+1, b.length()+1, 0.0);
    Matrix<size_t> vGaps(a.length()+1, b.length()+1, 0);
    Matrix<size_t> hGaps(a.length()+1, b.length()+1, 0);
    Matrix<char>   dirct(a.length()+1, b.length()+1, 'S');

    for (size_t i=0; i<=a.length(); ++i) {
	F(i, 0) = S.hAffineScore(i);
	hGaps(i, 0) = i;
    }

    for (size_t j=0; j<=b.length(); ++j) {
	F(0, j) = S.vAffineScore(j);
	vGaps(0, j) = j;
    }

    for (size_t i = 1; i <= a.length(); ++i) {
	for (size_t j = 1; j <= b.length(); ++j) {
	    double diag = F(i-1, j-1) + S(a[i-1], b[j-1]);
	    double horz = F(i-1, j)   + S.hAffineScore(hGaps(i-1, j) + 1);
	    double vert = F(i, j-1)   + S.vAffineScore(vGaps(i, j-1) + 1);

	    if (diag > horz && diag > vert) {
		F(i, j) = diag;
		dirct(i, j) = 'D';
	    } else if (horz > vert) {
		F(i, j) = horz;
		hGaps(i, j) = hGaps(i-1, j) + 1;
		dirct(i, j) = 'H';
	    } else {
		F(i, j) = vert;
		vGaps(i, j) = vGaps(i, j-1) + 1;
		dirct(i, j) = 'V';
	    }
	}
    }

    size_t i = a.length();
    size_t j = b.length();

    while (i > 0 && j > 0)
    {
	switch (dirct(i, j)) {
	case 'D':
	    aln->push_back(b[j-1]);
	    --i;
	    --j;
	    break;
	case 'H':
	    aln->push_back(Dash);
	    --i;
	    break;
	case 'V':
	    //aln->push_back(Dash);
	    --j;
	    break;
	default:
	    assert(false);
	}
    }
    
    while (i > 0)
    {
	aln->push_back(Dash);
	--i;
    }

    /*while (j > 0)
      {
      //aln->push_back(b[j-1]);
      --j;
      }*/

    reverse(aln->begin(), aln->end());

    return aln;
}

template double 
alignmentScore<GenScores, GISeq, GISeq> 
(GenScores&, GISeq&, GISeq&, SiteInformation* si = NULL);

template vector<GeneticSymbols>*
nwAlignment<GenScores, GeneticSymbols> (GenScores& S, GISeq& a, GISeq& b);
