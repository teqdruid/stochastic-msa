#include "storage.h"

#include <set>
#include <sys/time.h>
#include <stdlib.h>

using namespace std;

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

	printf("%20s Timer (ms):\t%lf\n", name.c_str(), 1000*(end-start));
    }
};

static string longRndSeq() {
    size_t size = 5000000; //Is 1/2 a million big enough?
    char* buffer = (char*)malloc(size+1);
    
    for (size_t i=0; i<size; i++) {
	int r = rand() % 4;
	switch (r) {
	case 0:
	    buffer[i] = 'A';
	    break;
	case 1:
	    buffer[i] = 'C';
	    break;
	case 2:
	    buffer[i] = 'G';
	    break;
	case 3:
	    buffer[i] = 'T';
	    break;
	}
    }
    
    buffer[size] = 0;
    return buffer;
}

int main(void) {
    string seq = longRndSeq();
    size_t totalMutations = 100000;

    printf("Making %u mutations on a sequence of length %u\n",
	   totalMutations, seq.length());

    srand(time(NULL));
    
    GISeq is(seq);
    MGISeq ms(is);

    cout << endl;


    int change = 0;
    {Timer a("Mutation");
	size_t len = is.length();
	set<int> touched;
	for (size_t i=0; i<totalMutations; i++) {
	    size_t loc = rand() % len;
	    size_t insDel = rand() % 3;
	    GeneticSymbols t = (GeneticSymbols) (rand() & 0b11);

	    if (touched.count(loc)) {
		--i;
		continue;
	    }

	    touched.insert(loc);

	    switch(insDel) {
	    case 0:
		//printf("Set %u to %c\n", loc, toChar(t));
		ms.set(loc, t);
		break;
	    case 1:
		//printf("Delete %u\n", loc);
		ms.del(loc);
		change--;
		break;
	    case 2:
		//printf("Insert %u, %c\n", loc, toChar(t));
		ms.insert(loc, t);
		change++;
		break;
	    }
	}
    }

    size_t newLength;

    {Timer a("Length calc");
	newLength = ms.length();
    }

    GISeq* ns;
    {Timer a("Commit");
	ns = ms.commit();
    }

    cout << endl << endl;
    cout << "New length:\t\t" << newLength << endl;
    cout << "\tCalc:\t\t" << change + is.length() << endl;
    cout << "Committed length:\t" << ns->length() << endl;
}
