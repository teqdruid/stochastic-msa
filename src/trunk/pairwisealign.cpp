#include "algo.h"
#include "storage.h"
#include <iostream>

// The flag parameter tells us whether we need the alignment to be computed.
// The direction of the alignment is stored as follows:
// S: Start of alignment
// D: Diagonal
// L: Left
// U: Up
template<class M, class A, class B>
double alignmentScore(M& m, A& a, B& b)
{
	//Visualize 'b' along the row.
	double scores[b.length()+1];
	int verticalGaps[b.length()+1];
	int horizontalGaps[b.length()+1];
	double tempScore;

	//Initialization;
	scores[0] = 0;
	verticalGaps[0] = 0;
	for(size_t i = 1; i < b.length()+1; i++)
	{
		scores[i] = m.affineScore(i);
		verticalGaps[i] = 0;
	}

	cout << "For testing purposes\n";
	for(size_t i = 0; i < b.length()+1; i++)
	{
		cout << scores[i] << "\t";
	}
	cout << "\n";

	//Recursion;
	for(size_t i = 1; i < a.length()+1; i++)
	{
		tempScore = m.affineScore(i);
		horizontalGaps[0] = 0;
		for(size_t j = 1; j < b.length()+1; j++)
		{
			horizontalGaps[j] = horizontalGaps[j-1] + 1;
			verticalGaps[j] = 0;
			double diagonalScore = scores[j-1] + m.score(b.at(j-1), a.at(i-1));
			double verticalScore = scores[j] + m.affineScore(verticalGaps[j]+1);
			double horizontalScore = tempScore + m.affineScore(horizontalGaps[j]);
			double maxScore = horizontalScore;
			if(diagonalScore > maxScore)
			{
				maxScore = diagonalScore;
				horizontalGaps[j] = 0;
			}
			if(verticalScore > maxScore)
			{
				maxScore = verticalScore;
				horizontalGaps[j] = horizontalGaps[j-1] - 1;
				verticalGaps[j] = verticalGaps[j]+1;
			}
			scores[j-1] = tempScore;
			tempScore = maxScore;
		}
		scores[b.length()] = tempScore;

		//For testing purposes.
		for(size_t i = 0; i < b.length()+1; i++)
		{
			cout << scores[i] << "\t";
		}
		cout << "\n";
	}

	//Find the maximum score.
	size_t i = 1;
	double maxScore = scores[0];
	while(i < b.length()+1)
	{
		if(scores[i] > maxScore)
		{
			maxScore = scores[i];
		}
		i++;
	}
	cout << "Returning max score.\n";
	return maxScore;
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

template<class A, class B>
void reconstructAlignment(A& a, B& b, const char*** directions)
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
	cout << aAlignment << "\n";
	cout << bAlignment << "\n";
}


//template double alignmentScore<GenScores, GISeq, GISeq>(GenScores&, GISeq&, GISeq&);

