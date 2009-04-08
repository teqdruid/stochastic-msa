#include <string>
#include <fstream>
#include <algorithm>
#include <msa.h>
#include <algo.h>

template<class T>
void MSA<T>::read(istream& is) {
    //Do input reading stuff here.
}

template<class T>
void MSA<T>::execute() {
    {
	set<ImmutableSequence<T>*> iSet = generator->generateSet(*this);
	this->scoreAddProfiles(iSet);
    }

    while (!terminator->terminate(*this)) {
	set<ImmutableSequence<T>*> newSet =
	    this->mutator->mutate(*this);

	this->scoreAddProfiles(newSet);

	vector< pair<double, ImmutableSequence<T>* > > selset = 
	    this->selector->select(*this, this->K);

	this->profiles.swap(selset);

	for (size_t i=0; i<selset.size(); i++) {
	    if (find(this->profiles.begin(), this->profiles.end(), selset[i])
		== this->profiles.end()) {
		//This profile was not selected
		delete selset[i].second;
	    }
	}
    }
}

template<class T>
void MSA<T>::scoreAddProfiles(set<ImmutableSequence<T>*>& profs) {
    typeof(profs.begin()) iter;

    for (iter = profs.begin(); iter != profs.end(); iter++) {
	pair<double, ImmutableSequence<T>*> p;
	p.second = *iter;
	p.first = this->scorer->score(**iter, *this);

	this->profiles.push_back(p);
    }
}

template<class T>
pair<double, ImmutableSequence<T>* > MSA<T>::best() {
    pair<double, ImmutableSequence<T>* > high;

    high.first = -99999999.9;

    for (size_t i=0; i<profiles.size(); ++i) {
	if (profiles[i].first > high.first)
	    high = profiles[i];
    }

    return high;
}

template<class T>
MSA<T>::MSA():
    scores(NULL), scorer(NULL), selector(NULL), mutator(NULL),
    terminator(NULL), generator(NULL), K(0)
{}

/*************
 *   Plugin implementations
 ****************/

template<class T>
class StarScore: public Scorer<T> {
    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa) {
	double total = 0.0;

	for (size_t i = 0; i<msa.sequences.size(); ++i) {
	    total += alignmentScore(*msa.scores, profile, *msa.sequences[i]);
	}

	return total;
    }
};



//Used to sort profiles by their score
template<class T>
bool pairCmp( pair<double, ImmutableSequence<T>* > a,
	  pair<double, ImmutableSequence<T>* > b ) {
    return a.first > b.first;
}


//This guys just selects the K best profiles
template<class T>
class HighSelector: public Selector<T> {
    vector< pair<double, ImmutableSequence<T>* > >
    select(MSA<T>& msa, size_t k) {
	vector< pair<double, ImmutableSequence<T>* > > ret;
	sort(msa.profiles.begin(), msa.profiles.end(), pairCmp<T>);

	for (size_t i=0; i<k; ++i) {
	    ret.push_back(msa.profiles[i]);
	}

	return ret;
    }
};

template<class T>
class RandomPicker: public Generator<T> {
public:
    RandomPicker() { srand(time(NULL)); }
    virtual set<ImmutableSequence<T>*> generateSet(MSA<T>& msa) {
	set<ImmutableSequence<T>*> ret;

	for (size_t i=0; i<msa.K; ++i) {
	    size_t s = rand() % msa.sequences.size();
	    ret.insert(msa.sequences[s]);
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
	return count++ >= max;
    };
};

/********
 *  Primary execution area
 ******/

int main(int argv, char** argc) {
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


    //Fill in defaults
    if (msa.K == 0)
	msa.K = 10;

    if (msa.scores == NULL)
	msa.scores = new GenScores(.8, .3, 1, .1);
    if (msa.scorer == NULL)
	msa.scorer = new StarScore<GeneticSymbols>();
    if (msa.selector == NULL)
	msa.selector = new HighSelector<GeneticSymbols>();
    if (msa.generator == NULL)
	msa.generator = new RandomPicker<GeneticSymbols>();
    if (msa.terminator == NULL)
	msa.terminator = new IterationsTerminator<GeneticSymbols>(10);

    ifstream inp(filename.c_str(), ios::in);
    msa.read(inp);

    if (msa.sequences.size() < 3) {
	cerr << "Error: only found " <<
	    msa.sequences.size() << " sequences in file." << endl;
	return 1;
    }

    msa.execute();

    cout << "Score of best alignment found: " << msa.best().first << endl;
    
    return 0;
}
