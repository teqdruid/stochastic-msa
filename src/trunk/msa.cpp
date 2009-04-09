#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <msa.h>
#include <algo.h>

template<class T>
void MSA<T>::read(istream& is) {
    double count = 0;
    while (!is.eof()) {
	ImmutableSequence<T>* seq = new ImmutableSequence<T>(is);
	if (seq->length() > 0) {
	    this->sequences.push_back(seq);
	    count += seq->length();

	    printf("Read %10u residues from %s\n", seq->length(),
		   seq->identifier.substr(0, 50).c_str());
	}
    }

    printf("-----------\nRead %u sequences of average length %lf\n",
	   this->sequences.size(), count / ((double)this->sequences.size()));
}

template<class T>
void MSA<T>::execute() {
    {
	vector<ImmutableSequence<T>*> iSet = generator->generateSet(*this);
	this->scoreAddProfiles(iSet);
    }

    while (!terminator->terminate(*this)) {
	if (this->scorer->rescore()) {
	    int size = this->profiles.size();

#pragma omp parallel for
	    for (int i=0; i<size; i++) {
		pair<double, ImmutableSequence<T>*>& p = this->profiles[i];
		p.first = this->scorer->score(*p.second, *this);
	    }
	}
 
	vector<ImmutableSequence<T>*> newSet =
	    this->mutator->mutate(*this);

	this->scoreAddProfiles(newSet);

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
void MSA<T>::scoreAddProfiles(vector<ImmutableSequence<T>*>& profs) {
    int size = profs.size();

#pragma omp parallel for
    for (int i=0; i<size; i++) {
	pair<double, ImmutableSequence<T>*> p;
	p.second = profs[i];
	p.first = this->scorer->score(*p.second, *this);
	
#pragma omp critical
	this->profiles.push_back(p);
    }
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

    return *profiles.end();
}

template<class T>
MSA<T>::MSA():
    scores(NULL), scorer(NULL), selector(NULL), mutator(NULL),
    terminator(NULL), generator(NULL), K(0)
{}

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

	return total / ((double)seqSize);;
    }

    virtual bool rescore() { return false; }
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
};


//This guys just selects the K best profiles
template<class T>
class HighSelector: public Selector<T> {
public:
    vector< pair<double, ImmutableSequence<T>* > >
    select(MSA<T>& msa, size_t k) {
	sort(msa.profiles.begin(), msa.profiles.end(), pairCmp<T>);

	vector< pair<double, ImmutableSequence<T>* > > 
	    ret(msa.profiles.begin(), msa.profiles.begin() + k);

	msa.profiles.erase(msa.profiles.begin(), msa.profiles.begin() + k);

	assert(ret.size() == k);

	return ret;
    }
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

	double totalScore = 0;
	for (size_t i=0; i<msa.profiles.size(); i++)
	    totalScore += msa.profiles[i].first;
	
	for (size_t j=0; j<k; ++j) {
	    double sel = random() % ((long int)round(totalScore));
	    double count = 0.0;
	    for (size_t i=0; i<msa.profiles.size(); i++) {
		count += msa.profiles[i].first;
		if (count >= sel) {
		    totalScore -= msa.profiles[i].first;
		    ret.push_back(msa.profiles[i]);
		    msa.profiles.erase(msa.profiles.begin() + i);
		    break;
		}
	    }
	}

	assert(ret.size() == k);

	return ret;
    }
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
};

template<class T>
class IterationsTerminator: public Terminator<T> {
public:
    size_t count;
    size_t max;

    IterationsTerminator(size_t max): count(0), max(max) {}
    virtual bool terminate(MSA<T>& msa) {
	cout << "Starting iteration: " << count 
	     << " best score: " << msa.best().first
	     << " range: " << (double)(msa.best().first - msa.worst().first)
	     << endl;
	return count++ >= max;
    };
};

template<class T>
class TotallyRandomMutator : public Mutator<T> {
public:
    long outputs;
    long mutations;
    TotallyRandomMutator(size_t outputs, size_t mutations) :
	outputs(outputs), mutations(mutations) {
	srand(time(NULL));
    }

    virtual vector<ImmutableSequence<T>*> mutate(MSA<T>& msa) {
	vector<ImmutableSequence<T>*> ret;

#pragma omp parallel for shared(ret)
	for (long i=0; i<outputs; ++i) {
	    ImmutableSequence<T>* is =
		msa.sequences[rand() % msa.sequences.size()];
	    MutableSequence<T> m(*is);

	    for (long j=0; j<mutations; j++) {
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
	    }

#pragma omp critical
	    ret.push_back(m.commit());
	}
	return ret;
    }
};

/********
 *  Primary execution area
 ******/

double get_runtime(void)
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

static ImmutableSequence<GeneticSymbols>* longRndSeq(size_t size) {
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

int main(int argv, char** argc) {
    srand(time(NULL));

    MSA<GeneticSymbols> msa;

    string filename;
    for (int i=1; i<argv; ++i) {
	string opt = argc[1];
	if (opt[0] != '-') {
	    filename = opt;
	    continue;
	}
    }

    if (filename == "") {
	cout << "Usage: " << argc[0] << " [inputFile] (options)" << endl;
	return 1;
    }

    {Timer a("Read data");
	ifstream inp(filename.c_str(), ios::in);
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

    //Fill in defaults
    if (msa.K == 0)
	msa.K = 2;

    if (msa.scores == NULL)
	msa.scores = new GenScores(.8, .3, 1, .1);
    if (msa.scorer == NULL)
	msa.scorer = new StarScore<GeneticSymbols>();
    if (msa.selector == NULL)
	//msa.selector = new HighSelector<GeneticSymbols>();
	msa.selector = new StochasticSelector<GeneticSymbols>();
    if (msa.generator == NULL)
	msa.generator = new RandomPicker<GeneticSymbols>();
    if (msa.terminator == NULL)
	msa.terminator = new IterationsTerminator<GeneticSymbols>(10);
    if (msa.mutator == NULL)
	msa.mutator = new TotallyRandomMutator<GeneticSymbols>(25, 50);

    {Timer a("Execution");
	msa.execute();
    }

    StarScore<GeneticSymbols> ss;

    {Timer a("Score compute");
	cout << "Star score of best alignment found: " << 
	    ss.score(*msa.best().second, msa) << endl;
    }
    
    return 0;
}
