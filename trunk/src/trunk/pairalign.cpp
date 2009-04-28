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

    string filename = "", outFilename = "";
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

	if (opt == "-out") {
	    ++i;
	    if (i >= argv) {
		cerr << "Need argument after -out" << endl;
		return 1;
	    }
	    outFilename = argc[i];
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

    msa.scores = new GenScores(1, -5, 5.0, 1, 25, 5);

    double score;
    {Timer a("Pairwise score compute");
	score = alignmentScore(*msa.scores, *msa.sequences[seqA], *msa.sequences[seqB]);
    }

    if (outFilename != "") {
	ofstream outp(outFilename.c_str(), ios::out);
	msa.sequences[seqA]->write(outp);
	outp << endl;
	msa.sequences[seqB]->write(outp);
	outp << endl;

	Timer a("Pairwise compute");
//	const char*** aln = 
//    getAlignment(*msa.scores, *msa.sequences[seqA], *msa.sequences[seqB]);
	
	outp << ">Aligned";

	vector<GeneticSymbols>* align =
	    nwAlignment(*msa.scores, *msa.sequences[seqA], *msa.sequences[seqB]);
	for (size_t i=0; i<align->size(); ++i) {
	    if ((i % 70) == 0)
		outp << endl;
	    outp << toChar(align->at(i));
	}
	outp << endl;

	delete align;
//	reconstructAlignment(outp, *msa.sequences[seqA], *msa.sequences[seqB], aln);
    }
    

    cout << "Pairwise alignment score: " << score << endl;
    
    return 0;
}
