/*
** msa.h
** 
*/

#ifndef   	MSA_H_
# define   	MSA_H_

#include <storage.h>

#include <set>
#include <vector>

template<class T> class MSA;

template<class T>
class Scorer {
public:
    virtual bool rescore() = 0; //Do we have to re-score each round?
    virtual double score(ImmutableSequence<T>& profile, MSA<T>& msa) = 0;
};

template<class T>
class Selector {
public:
    virtual vector< pair<double, ImmutableSequence<T>* > >
	select(MSA<T>& msa, size_t k) = 0;
};

template<class T>
class Mutator {
public:
    virtual vector<ImmutableSequence<T>*>
	mutate(MSA<T>& msa) = 0;
};

template<class T>
class Terminator {
public:
    virtual bool terminate(MSA<T>& msa) = 0;
};

template<class T>
class Generator {
public:
    virtual vector<ImmutableSequence<T>*> generateSet(MSA<T>& msa) = 0;
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
    
    vector< ImmutableSequence<T>* > sequences;
    vector< pair<double, ImmutableSequence<T>* > > profiles;

public:
    MSA();
    ~MSA();

    void read(istream& is);
    void execute();

    pair<double, ImmutableSequence<T>* > best();
    void printFinalAlign(ostream* os);

private:
    void scoreAddProfiles(vector<ImmutableSequence<T>*>& profs);
};




#endif 	    /* !MSA_H_ */
