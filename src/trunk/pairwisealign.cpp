#include "algo.h"
#include "storage.h"
#include <iostream>
#include <algorithm>

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

    T& operator()(size_t x, size_t y) {
	return mat[x + y*xs];
    }
};


template<class M, class A, class B>
double alignmentScore(M& S, A& a, B& b)
{
    Matrix<double> F(a.length()+1, b.length()+1, 0.0);
    Matrix<size_t> vGaps(a.length()+1, b.length()+1, 0);
    Matrix<size_t> hGaps(a.length()+1, b.length()+1, 0);
    Matrix<char>   dirct(a.length()+1, b.length()+1, 'S');

    for (size_t i=0; i<=a.length(); ++i) {
	F(i, 0) = S.affineScore(i);
	hGaps(i, 0) = i;
    }

    for (size_t j=0; j<=b.length(); ++j) {
	F(0, j) = S.affineScore(j);
	vGaps(0, j) = j;
    }

    for (size_t i = 1; i <= a.length(); ++i) {
	for (size_t j = 1; j <= b.length(); ++j) {
	    double diag = F(i-1, j-1) + S(a[i-1], b[j-1]);
	    double horz = F(i-1, j)   + S.affineScore(hGaps(i-1, j) + 1);
	    double vert = F(i, j-1)   + S.affineScore(vGaps(i, j-1) + 1);

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

    return F(a.length(), b.length());
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
	F(i, 0) = S.affineScore(i);
	hGaps(i, 0) = i;
    }

    for (size_t j=0; j<=b.length(); ++j) {
	F(0, j) = S.affineScore(j);
	vGaps(0, j) = j;
    }

    for (size_t i = 1; i <= a.length(); ++i) {
	for (size_t j = 1; j <= b.length(); ++j) {
	    double diag = F(i-1, j-1) + S(a[i-1], b[j-1]);
	    double horz = F(i-1, j)   + S.affineScore(hGaps(i-1, j) + 1);
	    double vert = F(i, j-1)   + S.affineScore(vGaps(i, j-1) + 1);

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

template<class M, class A, class B>
const char*** getAlignment(M& m, A& a, B& b)
{
	//Visualize 'b' along the row.
	double scores[b.length()+1];
	int verticalGaps[b.length()+1];
	int horizontalGaps[b.length()+1];
	const char*** directions;

	directions = new const char**[a.length()+1];
	for(size_t i = 0; i < a.length() + 1; i++)
	{
		directions[i] = new const char*[b.length()+1];
	}
	double tempScore;

	//Initialization;
	scores[0] = 0;
	directions[0][0] = "S";
	verticalGaps[0] = 0;
	for(size_t i = 1; i < b.length()+1; i++)
	{
		scores[i] = m.affineScore(i);
		directions[0][i] = "L";
		verticalGaps[i] = 0;
	}

	//Recursion;
	for(size_t i = 1; i < a.length()+1; i++)
	{
		tempScore = m.affineScore(i);
		directions[i][0] = "U";
		horizontalGaps[0] = 0;
		for(size_t j = 1; j < b.length()+1; j++)
		{
			horizontalGaps[j] = horizontalGaps[j-1] + 1;
			verticalGaps[j] = 0;
			double diagonalScore = scores[j-1] + m.score(b.at(j-1), a.at(i-1));
			double verticalScore = scores[j] + m.affineScore(verticalGaps[j]+1);
			double horizontalScore = tempScore + m.affineScore(horizontalGaps[j]);
			double maxScore = horizontalScore;
			const char* direction = "L";
			if(diagonalScore > maxScore)
			{
				maxScore = diagonalScore;
				horizontalGaps[j] = 0;
				direction = "D";
			}
			if(verticalScore > maxScore)
			{
				maxScore = verticalScore;
				horizontalGaps[j] = horizontalGaps[j-1] - 1;
				verticalGaps[j] = verticalGaps[j]+1;
				direction = "U";
			}
			scores[j-1] = tempScore;
			tempScore = maxScore;
			directions[i][j] = direction;
		}
		scores[b.length()] = tempScore;
	}
	return directions;
}

template<class A>
void freeAlignment(const char*** directions, A& a) {
    for(size_t i = 0; i < a.length() + 1; i++)
    {
	delete directions[i];
    }
    delete directions;
}

template<class A, class B>
void reconstructAlignment(ostream& os, A& a, B& b, const char*** directions)
{
	const char* C = "A";
	string aAlignment = "";
	string bAlignment = "";
	int aIndex = a.length();
	int bIndex = b.length();

	const char* S = "S";
	const char* U = "U";
	const char* D = "D";
	const char* L = "L";

	while(strcmp(C, S))
	{
		if(!strcmp(directions[aIndex][bIndex], U))
		{
			aIndex--;
			aAlignment = toChar(a.at(aIndex)) + aAlignment;
			bAlignment = "-" + bAlignment;
		}
		else if(!strcmp(directions[aIndex][bIndex], D))
		{
			aIndex--;
			bIndex--;
			aAlignment = toChar(a.at(aIndex)) + aAlignment;
			bAlignment = toChar(b.at(bIndex)) + bAlignment;
		}
		else if(!strcmp(directions[aIndex][bIndex], L))
		{
			bIndex--;
			bAlignment = "-" + bAlignment;
			bAlignment = toChar(b.at(bIndex)) + bAlignment;
		}
		else
		{
			C = "S";
		}
	}
	//os << aAlignment << "\n";

	for (size_t i = 0; i < bAlignment.length(); ++i) {
	    if ((i % 70) == 0)
		os << endl;
	    os << bAlignment[i];
	}
	os << endl;
}


template double alignmentScore<GenScores, GISeq, GISeq>(GenScores&, GISeq&, GISeq&);

template const char*** getAlignment<GenScores, GISeq, GISeq>(GenScores&, GISeq&, GISeq&);

template void reconstructAlignment<GISeq, GISeq>(ostream&, GISeq&, GISeq&, const char***);

template void freeAlignment<GISeq>(const char*** directions, GISeq& a);

template vector<GeneticSymbols>*
nwAlignment<GenScores, GeneticSymbols> (GenScores& S, GISeq& a, GISeq& b);
