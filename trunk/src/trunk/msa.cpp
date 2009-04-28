#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <msa.h>
#include <algo.h>
#include "util.h"

template<class T>
void MSA<T>::read(istream& is) {
    double count = 0;
    while (!is.eof()) {
	ImmutableSequence<T>* seq = new ImmutableSequence<T>(is);
	if (seq->length() > 0) {
	    this->sequences.push_back(seq);
	    count += seq->length();

	    printf("Read %10u residues from %s\n", seq->length(),
		   seq->identifier.c_str());
	} else {
	    delete seq;
	}
    }

    printf("-----------\nRead %u sequences of average length %lf\n",
	   this->sequences.size(), count / ((double)this->sequences.size()));
}

template<class T>
void MSA<T>::execute() {
    {
	cout << "Generating initial profile set..." << endl;
	vector<ImmutableSequence<T>*> iSet = generator->generateSet(*this);
	cout << "Scoring initial profiles..." << endl;
	this->scoreAddProfiles(iSet);
    }

    while (!terminator->terminate(*this)) {
	Timer a("Iteration");
	if (this->scorer->rescore()) {
	    cout << "\tRescoring profiles..." << endl;
	    int size = this->profiles.size();

#pragma omp parallel for
	    for (int i=0; i<size; i++) {
		pair<double, ImmutableSequence<T>*>& p = this->profiles[i];
		p.first = this->scorer->score(*p.second, *this);
	    }
	}
 
	cout << "\tRunning successor function..." << endl;
	vector<ImmutableSequence<T>*> newSet =
	    this->mutator->mutate(*this);

	cout << "\tScoring new profiles..." << endl;
	this->scoreAddProfiles(newSet);

	cout << "\tSelecting new profile set..." << endl;
	vector< pair<double, ImmutableSequence<T>* > > selset = 
	    this->selector->select(*this, this->K);

	this->profiles.swap(selset);

	for (size_t i=0; i<selset.size(); i++) {
	    //This profile was not selected
	    delete selset[i].second;
	}
    }
}

template<class T>
void MSA<T>::index() {
    
}

template<class T>
void MSA<T>::scoreAddProfiles(vector<ImmutableSequence<T>*>& profs) {
    int size = profs.size();

    cout << "\t\tScored: ";
    cout.flush();
#pragma omp parallel for
    for (int i=0; i<size; i++) {
	pair<double, ImmutableSequence<T>*> p;
	p.second = profs[i];
	p.first = this->scorer->score(*p.second, *this);
	cout << i << ":" << p.first << " ";
	cout.flush();
	
#pragma omp critical
	this->profiles.push_back(p);
    }

    cout << endl;
}

//Used to sort profiles by their score
template<class T>
bool pairCmp( pair<double, ImmutableSequence<T>* > a,
	  pair<double, ImmutableSequence<T>* > b ) {
    return a.first > b.first;
}

template<class T>
pair<double, ImmutableSequence<T>* > MSA<T>::best() {
    sort(profiles.begin(), profiles.end(), pairCmp<T>);
    return profiles[0];
}

template<class T>
pair<double, ImmutableSequence<T>* > MSA<T>::worst() {

    sort(profiles.begin(), profiles.end(), pairCmp<T>);

    return *(profiles.end() - 1);
}

template<class T>
MSA<T>::~MSA() {
    for (typeof(sequences.begin()) iter = sequences.begin();
	 iter != sequences.end(); iter++) {
	delete *iter;
    }

    for (typeof(profiles.begin()) iter = profiles.begin();
	 iter != profiles.end(); iter++) {
	delete iter->second;
    }

    delete scores;
    delete scorer;
    delete selector;
    delete mutator;
    delete terminator;
    delete generator;
}

template<class T>
void MSA<T>::output(ostream& os, ImmutableSequence<T>* profile) {
    for (size_t i=0; i<sequences.size(); i++) {
	ImmutableSequence<T>* seq = sequences[i];
	cout << "Writing results for " << seq->identifier << endl;
	os << ">" << seq->identifier << endl;
	vector<T>* align = nwAlignment(*this->scores, *profile, *seq);
	for (size_t i=0; i<align->size(); ++i) {
	    if ((i % 70) == 0)
		os << endl;
	    os << toChar(align->at(i));
	}
	os << endl << endl;
	delete align;
    }
}

/*************
 *   Plugin implementations
 ****************/

template<class T>
class StarScore: public Scorer<T> {
public:
    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa) {
	double total = 0.0;

	long seqSize = msa.sequences.size();
	for (long i = 0; i<seqSize; ++i) {
	    total += alignmentScore(*msa.scores, profile, *msa.sequences[i]);
	}

	return total / ((double)seqSize);
    }

    virtual bool rescore() { return false; }

    virtual void print() { cout << "Using star scorer" << endl; }
};

template<class T>
class SampledStarScore: public Scorer<T> {
public:
    long samples;

    SampledStarScore(long samples, MSA<T>& msa) : samples(samples) {
	srand(time(NULL));
    }

    SampledStarScore(double frac, MSA<T>& msa) {
	srand(time(NULL));
	samples = round(msa.sequences.size() * frac);
    }

    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa) {
	double total = 0.0;

	long seqSize = msa.sequences.size();
	for (long i = 0; i<samples; ++i) {
	    long rIndx = rand() % seqSize;
	    total += alignmentScore(*msa.scores, profile, *msa.sequences[rIndx]);
	}

	return total / ((double)samples);
    }

    virtual bool rescore() { return true; }

    virtual void print() { cout << "Using sampled star scorer, " 
				<< samples << " samples" << endl; }
};


//This guys just selects the K best profiles
template<class T>
class HighSelector: public Selector<T> {
public:
    vector< pair<double, ImmutableSequence<T>* > >
    select(MSA<T>& msa, size_t k) {
	sort(msa.profiles.begin(), msa.profiles.end(), pairCmp<T>);

	cout << "\t\tSelection set best: " << msa.profiles.begin()->first
	     << ", worst: " << (msa.profiles.end()-1)->first << endl;
	vector< pair<double, ImmutableSequence<T>* > > 
	    ret(msa.profiles.begin(), msa.profiles.begin() + k);

	msa.profiles.erase(msa.profiles.begin(), msa.profiles.begin() + k);

	assert(ret.size() == k);

	return ret;
    }

    virtual void print() { cout << "Using high-k selector" << endl; }
};

//This guys just selects profiles stochastically
template<class T>
class StochasticSelector: public Selector<T> {
public:
    StochasticSelector() {
	srandom(time(NULL));
    }

    //TODO: This now takes O(p*k) , p = #of profiles.  Do it faster
    vector< pair<double, ImmutableSequence<T>* > >
    select(MSA<T>& msa, size_t k) {
	vector< pair<double, ImmutableSequence<T>* > > ret;
	sort(msa.profiles.begin(), msa.profiles.end(), pairCmp<T>);

	double worst = (msa.profiles.end()-1)->first;

	double totalScore = 0;
	for (size_t i=0; i<msa.profiles.size(); i++)
	    totalScore += msa.profiles[i].first - worst;
	
	for (size_t j=0; j<k; ++j) {
	    double sel = random() % ((long int)round(totalScore));
	    double count = 0.0;
	    for (size_t i=0; i<msa.profiles.size(); i++) {
		count += msa.profiles[i].first - worst;
		if (count >= sel) {
		    totalScore -= msa.profiles[i].first - worst;
		    ret.push_back(msa.profiles[i]);
		    msa.profiles.erase(msa.profiles.begin() + i);
		    break;
		}
	    }
	}

	assert(ret.size() == k);

	return ret;
    }

    virtual void print() { cout << "Using stochastic selector" << endl; }
};

template<class T>
class RandomPicker: public Generator<T> {
public:
    RandomPicker() { srand(time(NULL)); }
    virtual vector<ImmutableSequence<T>*> generateSet(MSA<T>& msa) {
	vector<ImmutableSequence<T>*> ret;

	for (size_t i=0; i<msa.K; ++i) {
	    size_t s = rand() % msa.sequences.size();
	    ret.push_back(msa.sequences[s]->copy());
	}

	return ret;
    }

    virtual void print() { cout << "Using random picking generator" << endl; }
};

template<class T>
class RandomGenerator: public Generator<T> {
public:
    RandomGenerator() { srand(time(NULL)); }
    virtual vector<ImmutableSequence<T>*> generateSet(MSA<T>& msa) {
	vector<ImmutableSequence<T>*> ret;

	double total = 0.0;
	for (size_t i=0; i<msa.sequences.size(); i++) {
	    total += msa.sequences[i]->length();
	}
	double avg = total / msa.sequences.size();

	for (size_t i=0; i<msa.K; ++i) {
	    ret.push_back(longRndSeq(avg));
	}

	return ret;
    }

    virtual void print() { cout << "Using random generator" << endl; }
};

template<class T>
class IterationsTerminator: public Terminator<T> {
public:
    size_t count;
    size_t max;

    IterationsTerminator(size_t max): count(0), max(max) {}
    virtual bool terminate(MSA<T>& msa) {
	double bestScore = msa.best().first;
	double worstScore = msa.worst().first;
	cout << "Starting iteration: " << count 
	     << " best score: " << msa.best().first
	     << " range: " << bestScore - worstScore
	     << endl;
	cout << "\tScores: ";
	for (size_t i=0; i<msa.profiles.size(); i++) {
	    cout << i << ":" << msa.profiles[i].first << " ";
	}
	cout << endl;
	return count++ >= max;
    };

    virtual void print() { cout << "Using iterations (" << max << ") terminator" << endl; }
};

template<class T>
class TotallyRandomMutator : public Mutator<T> {
public:
    long outputs;
    double mutations;//Expected number of mutations
    TotallyRandomMutator(size_t outputs, size_t mutations) :
	outputs(outputs), mutations(mutations) {
	srand(time(NULL));
    }

    virtual vector<ImmutableSequence<T>*> mutate(MSA<T>& msa) {
	vector<ImmutableSequence<T>*> ret;

	size_t chanceEnd = (1.0 / mutations) * 1.6777216e7;

#pragma omp parallel for shared(ret)
	for (long i=0; i<outputs; ++i) {
	    ImmutableSequence<T>* is =
		msa.sequences[rand() % msa.sequences.size()];
	    MutableSequence<T> m(*is);

	    //Allow no more than 10*expected mutations
	    for (long j=0; j<(mutations*10); j++) {
		size_t loc = rand() % is->length();
		size_t insDel = rand() % 3;
		GeneticSymbols t = (GeneticSymbols) (rand() & 3);

		switch(insDel){
		case 0:
		    m.set(loc, t);
		    break;
		case 1:
		    m.del(loc);
		    break;
		case 2:
		    m.insert(loc, t);
		    break;
		} 

		size_t brk = random() & 0xFFFFFF;
		if (brk <= chanceEnd) {
		    //cout << "Rnd end!! " << j << endl;
		    break;
		}
	    }

	    ImmutableSequence<T>* cm = m.commit();

#pragma omp critical
	    ret.push_back(cm);
	}
	return ret;
    }

    virtual void print() { cout << "Using totally random mutator, outputs = "
				<< outputs << ", mutations ~= " << mutations << endl; }
};



int msa_main(int argv, char** argc) {
    srand(time(NULL));

    MSA<GeneticSymbols> msa;

    string filename = "", outputFile = "", profOutputFile = "";
    for (int i=1; i<argv; ++i) {
	string opt = argc[i];
	if (opt[0] != '-') {
	    filename = opt;
	    break;
	}
    }

    if (filename == "") {
	cout << "Usage: " << argc[0] << " [inputFile] (options)" << endl;
	return 1;
    }

    {Timer a("Read data");
	ifstream inp(filename.c_str(), ios::in);
	if (inp.is_open())
	    msa.read(inp);
	
	if (filename == "genRandom" && msa.sequences.size() == 0) {
#pragma omp parallel for
	    for (int i=0; i < 50; i++) {
		ImmutableSequence<GeneticSymbols>* n = longRndSeq(500);
		
#pragma omp critical
		msa.sequences.push_back(n);
	    }
	}
    }
    
    if (msa.sequences.size() < 3) {
	cerr << "Error: only found " <<
	    msa.sequences.size() << " sequences in file." << endl;
	return 1;
    }

    {Timer a("Indexing data");
	msa.index();
    }

    double a = 2.0, m = 300.0;

    for (int i=1; i<argv; ++i) {
	string opt = argc[i];

	if (opt == "-k") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -k" << endl;
		return 1;
	    }
	    opt = argc[i];
	    msa.K = atol(opt.c_str());
	}

	if (opt == "-a") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -a" << endl;
		return 1;
	    }
	    opt = argc[i];
	    a = atof(opt.c_str());
	}

	if (opt == "-m") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -m" << endl;
		return 1;
	    }
	    opt = argc[i];
	    m = atof(opt.c_str());
	}

	if (opt == "-starscore") {
	    msa.scorer = new StarScore<GeneticSymbols>();
	}

	if (opt == "-samplescore") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -samplescore" << endl;
		return 1;
	    }
	    opt = argc[i];
	    double r = atof(opt.c_str());
	    msa.scorer = new SampledStarScore<GeneticSymbols>(r, msa);
	}

	if (opt == "-iter") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -iter" << endl;
		return 1;
	    }
	    opt = argc[i];
	    long i = atol(opt.c_str());
	    msa.terminator = new IterationsTerminator<GeneticSymbols>(i);
	}

	if (opt == "-highsel") {
	    msa.selector = new HighSelector<GeneticSymbols>();
	}

	if (opt == "-randgen") {
	    msa.generator = new RandomGenerator<GeneticSymbols>();
	}

	if (opt == "-out") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -out" << endl;
		return 1;
	    }
	    opt = argc[i];
	    outputFile = opt;
	}

	if (opt == "-pout") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -pout" << endl;
		return 1;
	    }
	    opt = argc[i];
	    profOutputFile = opt;
	}	
    }

    //Fill in defaults
    if (msa.K == 0)
	msa.K = 10;

    if (msa.scores == NULL)
	msa.scores = new GenScores(1, -.5, 5.0, 1.6);

    if (msa.scorer == NULL)
	msa.scorer = new StarScore<GeneticSymbols>();
    if (msa.selector == NULL)
	msa.selector = new HighSelector<GeneticSymbols>();
    if (msa.generator == NULL)
	msa.generator = new RandomPicker<GeneticSymbols>();
    if (msa.terminator == NULL)
	msa.terminator = new IterationsTerminator<GeneticSymbols>(25);
    if (msa.mutator == NULL)
	msa.mutator = new TotallyRandomMutator<GeneticSymbols>(msa.K*a, m);

    cout << "K = " << msa.K << endl;
    msa.scorer->print();
    msa.selector->print();
    msa.generator->print();
    msa.mutator->print();
    msa.terminator->print();

    {Timer a("Execution");
	msa.execute();
    }

    StarScore<GeneticSymbols> ss;

    {Timer a("Score compute");
	cout << "Star score of best alignment found: " << 
	    ss.score(*msa.best().second, msa) << endl;
    }

    if (profOutputFile != "") {
	msa.best(); //Sort the entries
	ofstream outp(profOutputFile.c_str(), ios::out);

	for (size_t i=0; i<msa.profiles.size(); ++i) {
	    char buffer[128];
	    GISeq* p = msa.profiles[i].second;

	    snprintf(buffer, 128, "Profile %d scoring %lf",
		     i, msa.profiles[i].first);

	    p->identifier = buffer;
	    p->write(outp);
	}
    }

    if (outputFile != "") {
	Timer a("Final alignment");
	cout << "Doing final alignment and writing data to "
	     << outputFile << endl;

	ofstream outp(outputFile.c_str(), ios::out);
	msa.output(outp, msa.best().second);
    }
    
    return 0;
}

MSA<GeneticSymbols> __cpp_sux1;

