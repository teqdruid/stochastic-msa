#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "msa.h"
#include "algo.h"
#include "util.h"

int main(int argv, char** argc) {
    srand(time(NULL));

    MSA<GeneticSymbols> msa;

    string filename = "";
    long seqA = -1, seqB = -1;
    for (int i=1; i<argv; ++i) {
	string opt = argc[i];
	if (opt[0] != '-') {
	    if (filename == "")
		filename = opt;
	    else if (seqA == -1)
		seqA = atol(opt.c_str());
	    else if (seqB == -1)
		seqB = atol(opt.c_str());
	    continue;
	}
    }

    if (seqB == -1)
	seqB = seqA;

    if (filename == "" || seqA < 0 || seqB < 0) {
	cout << "Usage: " << argc[0] << " [inputFile] [seq num] [seq num]" << endl;
	return 1;
    }

    {Timer a("Read data");
	ifstream inp(filename.c_str(), ios::in);
	if (inp.is_open())
	    msa.read(inp);
	
	if (filename == "genRandom" && msa.sequences.size() == 0) {
	    ImmutableSequence<GeneticSymbols>* a = longRndSeq(seqA);
	    ImmutableSequence<GeneticSymbols>* b = longRndSeq(seqB);

	    assert(a->length() == (size_t)seqA);
	    assert(b->length() == (size_t)seqB);
		
	    msa.sequences.push_back(a);
	    msa.sequences.push_back(b);

	    seqA = 0;
	    seqB = 1;
	}
    }
    
    if ((size_t)seqA >= msa.sequences.size() || (size_t)seqB >= msa.sequences.size()) {
	cerr << "Error: only found " <<
	    msa.sequences.size() << " sequences in file." << endl;
	return 1;
    }

    msa.scores = new GenScores(1, -10, 15.0, 6.66);

    double score;
    {Timer a("Pairwise compute");
	score = alignmentScore(*msa.scores, *msa.sequences[seqA], *msa.sequences[seqB]);
    }

    cout << "Pairwise alignment score: " << score << endl;
    
    return 0;
}
