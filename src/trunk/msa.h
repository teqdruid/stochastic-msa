/*
** msa.h
** 
*/

#ifndef   	MSA_H_
# define   	MSA_H_

#include <storage.h>

#include <set>
#include <vector>

#define Ktup 4

template<class T> class MSA;

template<class T>
class Scorer {
public:
    virtual bool rescore() = 0; //Do we have to re-score each round?
    virtual double score(ImmutableSequence<T>& profile,
			 MSA<T>& msa,
			 SiteInformation* si) = 0;
    virtual void print() { cout << "Using unknown scorer" << endl; }
};

template<class T>
class Selector {
public:
    virtual vector< pair<double, ImmutableSequence<T>* > >
	select(MSA<T>& msa, size_t k) = 0;
    virtual void print() { cout << "Using unknown selector" << endl; }
};

template<class T>
class Mutator {
public:
    virtual vector<ImmutableSequence<T>*>
	mutate(MSA<T>& msa) = 0;
    virtual void print() { cout << "Using unknown mutator" << endl; }
};

template<class T>
class Terminator {
public:
    virtual bool terminate(MSA<T>& msa) = 0;
    virtual void print() { cout << "Using unknown terminator" << endl; }
};

template<class T>
class Generator {
public:
    virtual vector<ImmutableSequence<T>*> generateSet(MSA<T>& msa) = 0;
    virtual void print() { cout << "Using unknown generator" << endl; }
};

template<class T>
class MSA {
public:
    ScoringMatrix<T, 4>* scores;

    Scorer<T>* scorer;
    Selector<T>* selector;
    Mutator<T>* mutator;
    Terminator<T>* terminator;
    Generator<T>* generator;

    size_t K;
    double avgSize;

    map<T[Ktup], size_t> seqIndex;

    vector< ImmutableSequence<T>* > sequences;
    vector< pair<double, ImmutableSequence<T>* > > profiles;
    map<ImmutableSequence<T>*, SiteInformation*> profileSiteInfo;

public:
    MSA() :
    scores(NULL), scorer(NULL), selector(NULL), mutator(NULL),
	terminator(NULL), generator(NULL), K(0) {}
    ~MSA();

    void read(istream& is);
    void index();
    void execute();
    void output(ostream& os, ImmutableSequence<T>* profile);

    pair<double, ImmutableSequence<T>* > best();
    pair<double, ImmutableSequence<T>* > worst();
    void printFinalAlign(ostream* os);

private:
    void scoreAddProfiles(vector<ImmutableSequence<T>*>& profs);
};


int msa_main(int argv, char** argc);


#endif 	    /* !MSA_H_ */
