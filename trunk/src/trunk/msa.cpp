#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <msa.h>
#include <algo.h>
#include "util.h"

static size_t SIZE_DIV = 10;

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

    avgSize = count / ((double)this->sequences.size());

    printf("-----------\nRead %u sequences of average length %lf\n",
	   this->sequences.size(), avgSize);
}

template<class T>
static void eliminateDuplicates(vector<ImmutableSequence<T>*>& vec, bool del = true) {
    bool dupes[vec.size()];
    memset(dupes, false, sizeof(bool) * vec.size());

    for (size_t i=0; i<vec.size(); ++i) {
	ImmutableSequence<T>* a = vec[i];

	if (dupes[i])
	    continue;

	for (size_t j=i+1; j<vec.size(); ++j) {
	    ImmutableSequence<T>* b = vec[j];

	    if (dupes[j])
		continue;

	    if (a->hash() != b->hash())
		continue;

	    if (a->isEqual(*b))
		dupes[j] = true;
		//dupes.insert(j);
	}
    }

    size_t removed = 0;
    for (long i=vec.size()-1; i >= 0; --i)  {
	if (dupes[i]) {
	    ImmutableSequence<T>* dup = vec[i];
	    vec.erase(vec.begin() + i);
	    if (del) delete dup;
	    removed++;
	}
    }

    cout << "\tRemoved " << removed << " duplicates" << endl;
}

template<class T>
void MSA<T>::execute() {
    size_t bigCount = 0;
    do { {
	    Timer a("Iteration");
	    
	    if (profiles.size() > 0)
		if (this->scorer->rescore()) {
		    cout << "\tRescoring profiles..." << endl;
		    int size = this->profiles.size();
		    
#pragma omp parallel for
		    for (int i=0; i<size; i++) {
			pair<double, ImmutableSequence<T>*>& p = this->profiles[i];
			
#pragma omp critical
			if (profileSiteInfo.count(p.second)) {
			    delete profileSiteInfo[p.second];
			    profileSiteInfo.erase(p.second);
			}
			
			SiteInformation* si = 
			    new SiteInformation(p.second->length(),
						sequences.size() / SIZE_DIV);
			p.first = this->scorer->score(*p.second, *this, si);
		    
#pragma omp critical
			profileSiteInfo[p.second] = si;
		    }
		}
 
	    vector<ImmutableSequence<T>*> newSet;
	    if (profiles.size() > 0) {
		cout << "\tRunning successor function..." << endl;
		newSet = this->mutator->mutate(*this);
	    } else {
		cout << "\tRunning generation function..." << endl;
		newSet = generator->generateSet(*this);
	    }

	    eliminateDuplicates(newSet);

	    if (B > 0.0) {
		//Use our bag of tricks to eliminate a few
		trickyEliminate(newSet);
	    }

	    cout << "\tScoring new profiles..." << endl;
	    this->scoreAddProfiles(newSet);

	    cout << "\tSelecting new profile set..." << endl;
	    vector< pair<double, ImmutableSequence<T>* > > selset = 
		this->selector->select(*this, this->K);

	    this->profiles.swap(selset);

	    vector<ImmutableSequence<T>*> delSet;	    

	    for (size_t i=0; i<selset.size(); i++) {
		ImmutableSequence<T>* p = selset[i].second;
		delSet.push_back(p);
	    }

	    eliminateDuplicates(delSet, false);

	    for (size_t i=0; i<delSet.size(); i++) {
		ImmutableSequence<T>* p = delSet[i];
		//This profile was not selected

		if (profileSiteInfo.count(p)) {
		    delete profileSiteInfo[p];
		    profileSiteInfo.erase(p);
		}
		delete p;
	    }

	    cout << endl;

	    double bestScore = best().first;
	    double worstScore = worst().first;
	    cout << "Ending iteration: " << bigCount++
		 << " best score: " << bestScore
		 << " range: " << bestScore - worstScore
		 << ".  Sizeof(best): " << best().second->length()
		 << endl;
	    cout << "\tScores: ";
	    for (size_t i=0; i<profiles.size(); i++) {
		cout << profiles[i].second->length() << ":" 
		     << profiles[i].first << " ";
	    }
	    cout << endl;
	}
	cout << "----------------------" << endl << endl;
    } while (!terminator->terminate(*this));
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

	SiteInformation* si =
	    new SiteInformation(p.second->length(),
				sequences.size() / SIZE_DIV);
	p.first = this->scorer->score(*p.second, *this, si);

#pragma omp critical
	profileSiteInfo[p.second] = si;

	cout << p.second->length() << ":" << p.first << " ";
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
    long seqSize = sequences.size();

#pragma omp parallel for
    for (long i=0; i<seqSize; i++) {
	ImmutableSequence<T>* seq = sequences[i];
	vector<T>* align = nwAlignment(*this->scores, *profile, *seq);

#pragma omp critical 
	{
	    cout << "Writing results for " << seq->identifier << endl;
	    os << ">" << seq->identifier << endl;
	    for (size_t i=0; i<align->size(); ++i) {
		if ((i % 70) == 0)
		    os << endl;
		os << toChar(align->at(i));
	    }
	    os << endl << endl;
	    delete align;
	}
    }
}

/*************
 *   Plugin implementations
 ****************/

template<class T>
class StarScore: public Scorer<T> {
public:
    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa, SiteInformation* si) {
	double total = 0.0;

	long seqSize = msa.sequences.size();
	for (long i = 0; i<seqSize; ++i) {
	    total += alignmentScore(*msa.scores, profile, *msa.sequences[i], si);
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

    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa, SiteInformation* si) {
	double total = 0.0;

	long seqSize = msa.sequences.size();
	for (long i = 0; i<samples; ++i) {
	    long rIndx = rand() % seqSize;
	    total += alignmentScore(*msa.scores, profile, *msa.sequences[rIndx], si);
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

	for (size_t i=0; i<(msa.K*msa.A); ++i) {
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

	for (size_t i=0; i<(msa.K*msa.A); ++i) {
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
	return count++ >= max;
    };

    virtual void print() { cout << "Using iterations (" << max << ") terminator" << endl; }
};

template<class T>
class ProgressTerminator: public Terminator<T> {
public:
    double percent;
    double lastBest;
    size_t max;
    size_t count;
    size_t bigCount;

    ProgressTerminator(size_t max, double percent):
	percent(percent), lastBest(-9999999.9), max(max), count(0), bigCount(0) {}
    virtual bool terminate(MSA<T>& msa) {
	double diff = msa.best().first - lastBest;
	diff = diff / lastBest;
	if (diff < percent)
	    count++;
	else
	    count = 0;
	
	lastBest = msa.best().first;
	return count >= max;
    };

    virtual void print() { cout << "Using progress (" << max << "x, "
				<< percent*100 << "%) terminator" << endl; }
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
		msa.profiles[rand() % msa.profiles.size()].second;
	    MutableSequence<T> m(*is);

	    //Allow no more than 10*expected mutations
	    for (long j=0; j<(mutations*10); j++) {
		size_t loc = rand() % is->length();
		size_t insDel = rand() % 3;
		GeneticSymbols t = (GeneticSymbols) (rand() & 3);

		switch(insDel) {
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

template<class T>
class StochasticMutator : public Mutator<T> {
public:
    long outputs;
    double mutations;//Expected number of mutations
    StochasticMutator(size_t outputs, size_t mutations) :
	outputs(outputs), mutations(mutations) {
	srand(time(NULL));
	srandom(time(NULL));
    }

    virtual vector<ImmutableSequence<T>*> mutate(MSA<T>& msa) {
	vector<ImmutableSequence<T>*> ret;

	double muts = 0.0;

#pragma omp parallel for shared(ret)
	for (long i=0; i<outputs; ++i) {
	    ImmutableSequence<T>* is =
		msa.profiles[rand() % msa.profiles.size()].second;
	    MutableSequence<T> m(*is);

	    SiteInformation* si;
#pragma omp critical
	    {
		assert(msa.profileSiteInfo.count(is) > 0);
		si = msa.profileSiteInfo[is];
	    }

/*	    if (i == 0) {
		cout << endl;
		for (size_t j=0; j<is->length(); ++j) {
		    cout << j << ": " 
			         << si->subst[j]
			 << ", " << si->ins[j] 
			 << ", " << si->dels[j] 
			 << endl;
		}
		} */

	    size_t iLen = is->length() * 3;
	    size_t total = 0;
	    for (size_t j=0; j<iLen; ++j) {
		total += si->subst[j];
	    }

	    double totalSites = msa.sequences.size() * msa.avgSize;

	    //Choose a number of mutations proportional to the number
	    //of "bad" sites
	    double mFrac = (((double)total) / totalSites);
	    long mut = round(mutations * mFrac);

	    for (long j=0; j<mut; j++) {
		size_t aLoc = stochasticSelect(si->subst, iLen);
		size_t loc = aLoc / 3;
		size_t insDel = aLoc - (loc * 3);

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

	    ImmutableSequence<T>* cm = m.commit();

#pragma omp critical
	    {
		ret.push_back(cm);
		muts += mut;
	    }
	}

	cout << "\t\tAverage mutations per profile: " << muts / outputs << endl;
	return ret;
    }

    virtual void print() { cout << "Using stochastic mutator, outputs = "
				<< outputs << ", mutations < " << mutations << endl; }
};




/*********************************
 *   Done with plugins
 *********************************/

//Implement our bag of tricks down here
template<class T>
void MSA<T>::trickyEliminate(vector<ImmutableSequence<T>*>& profs) {
    vector< pair<double, ImmutableSequence<T>* > > scores;

    SampledStarScore<T> samp(0.05, *this);
    long pSize = profs.size();
    long total = 0;

    cout << "\tQuick scoring: ";
    cout.flush();

#pragma omp parallel for shared(total)
    for (long i=0; i<pSize; ++i) {
	double sizeDiff = profs[i]->length() - avgSize;
	sizeDiff /= sequences.size();

	double score = samp.score(*profs[i], *this, NULL);
	score += score * sizeDiff;

	//cout << profs[i]->length() << ":" << score << " ";
	//cout.flush();

	pair<double, ImmutableSequence<T>* > p(score, profs[i]);
#pragma omp critical
	{
	    scores.push_back(p);
	    total ++;
	    if ((total % 100) == 0) {
		cout << total << ", ";
		cout.flush();
	    }
	}
    }

    cout << endl;

    sort(scores.begin(), scores.end(), pairCmp<T>);

    size_t save = K * A * B;
    while (profs.size() > save) {
	ImmutableSequence<T>* del = scores.back().second;
	scores.pop_back();
	profs.erase(find(profs.begin(), profs.end(), del));

	if (profileSiteInfo.count(del)) {
	    delete profileSiteInfo[del];
	    profileSiteInfo.erase(del);
	}
	delete del;
    }
}

//MSA's primary run function
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

    msa.A = 2.0;
    msa.B = 0.0;
    double m = 300.0;
    double hGap = 250;

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
	    msa.A = atof(opt.c_str());
	}

	if (opt == "-b") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -b" << endl;
		return 1;
	    }
	    opt = argc[i];
	    msa.B = atof(opt.c_str());
	}

	if (opt == "-hg") {
	    i++;
	    if (i >= argv) {
		cerr << "Need value after -hg" << endl;
		return 1;
	    }
	    opt = argc[i];
	    hGap = atof(opt.c_str());
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
	    SIZE_DIV = SIZE_DIV / r;
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

	if (opt == "-prog") {
	    i++;
	    if (i >= argv) {
		cerr << "Need 2 values after -prog" << endl;
		return 1;
	    }
	    opt = argc[i];
	    long max = atol(opt.c_str());

	    i++;
	    if (i >= argv) {
		cerr << "Need 2 values after -prog" << endl;
		return 1;
	    }
	    opt = argc[i];
	    double p = atof(opt.c_str());

	    msa.terminator = new ProgressTerminator<GeneticSymbols>(max, p/100);
	}


	if (opt == "-highsel") {
	    msa.selector = new HighSelector<GeneticSymbols>();
	}

	if (opt == "-randgen") {
	    msa.generator = new RandomGenerator<GeneticSymbols>();
	}

	if (opt == "-randmut") {
	    msa.mutator = new TotallyRandomMutator<GeneticSymbols>(msa.K*msa.A, m);
	}

	if (opt == "-stomut") {
	    msa.mutator = new StochasticMutator<GeneticSymbols>(msa.K*msa.A, m);
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
	msa.scores = new GenScores(1.9, 0, 10, 0.2, hGap, hGap/10);

    if (msa.scorer == NULL)
	msa.scorer = new StarScore<GeneticSymbols>();
    if (msa.selector == NULL)
	msa.selector = new HighSelector<GeneticSymbols>();
    if (msa.generator == NULL)
	msa.generator = new RandomPicker<GeneticSymbols>();
    if (msa.terminator == NULL)
	msa.terminator = new ProgressTerminator<GeneticSymbols>(3, .005);
    if (msa.mutator == NULL)
	msa.mutator = new StochasticMutator<GeneticSymbols>(msa.K*msa.A, m);

    cout << "K = " << msa.K << endl;
    cout << "A = " << msa.A << endl;
    cout << "Size_Div = " << SIZE_DIV << endl;
    if (msa.B > 0.0) {
	cout << "B = " << msa.B 
	     << " ... Using bag of tricks for quick filtering" << endl;
    }
    msa.scorer->print();
    msa.selector->print();
    msa.generator->print();
    msa.mutator->print();
    msa.terminator->print();

    cout << endl
	 << "Starting execution:" << endl
	 << "-----------------------" << endl;

    {Timer a("Execution");
	msa.execute();
    }

    StarScore<GeneticSymbols> ss;

    {Timer a("Score compute");
	cout << "Star score of best alignment found: " << 
	    ss.score(*msa.best().second, msa, NULL) << endl;
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

