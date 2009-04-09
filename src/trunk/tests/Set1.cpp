#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <ctype.h>
#include <set>
#include "storage.h"

static string longRndSeq() {
    size_t size = 500000; //Is 1/2 a million big enough?
    char buffer[size+1];
    
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

class ImmutableSequenceTest: public CppUnit::TestFixture  {

    CPPUNIT_TEST_SUITE( ImmutableSequenceTest );

    CPPUNIT_TEST( testParse );
    CPPUNIT_TEST( testMulti );

    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {

    }

    void tearDown() {

    }

    string normalizeStr(string& s) {
	char buf[s.length()+1];
	size_t size = 0;
	for (size_t i=0; i<s.length(); i++) {
	    char c = s[i];
	    if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		continue;
	    if (c == '>')
		break;
	    buf[size++] = toupper(c);
	}
	buf[size] = 0;
	return buf;
    }

#define parseTest(T, S)  {						\
	string s = S;							\
	T a(s);								\
	string u = normalizeStr(s);					\
	CPPUNIT_ASSERT(a.toString() == u);				\
	T b(s);								\
	CPPUNIT_ASSERT(a == b);						\
    }
   
    void testParse() {
	parseTest(GISeq, "AAACTTGGTC");
	parseTest(GISeq, "AA  AAA");
	parseTest(GISeq, "");
	parseTest(GISeq, "ttttttt");
	parseTest(GISeq, "AAaatctggGGctACc");
	parseTest(GISeq, "AcTg \n \tAATCct\r GGTC");

	try {
	    GISeq s("This will fail");
	    CPPUNIT_ASSERT(false);
	} catch (MsaException* e) {
	    //Success!
	}
    }

    void testMulti() {
	parseTest(GISeq, "ACGT  ATCGTA \t ATCGTA > ACGGAT");

	GISeq* seq1 = new GISeq("ACGT > AGGT");
	CPPUNIT_ASSERT_EQUAL(string("ACGT"), seq1->toString());

	GISeq* seq2 = new GISeq(">Hello\nACGT > AGGT");
	CPPUNIT_ASSERT_EQUAL(string("ACGT"), seq2->toString());
	CPPUNIT_ASSERT_EQUAL(string("Hello"), seq2->identifier);

	delete seq1;
	delete seq2;
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION( ImmutableSequenceTest );


class MutableSequenceTest: public CppUnit::TestFixture  {

    CPPUNIT_TEST_SUITE( MutableSequenceTest );

    //CPPUNIT_TEST( bigTest );
    CPPUNIT_TEST( changeTest );

    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {
	srand(time(NULL));
    }

    void tearDown() {

    }

    string normalizeStr(string& s) {
	char buf[s.length()+1];
	size_t size = 0;
	for (size_t i=0; i<s.length(); i++) {
	    char c = s[i];
	    if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		continue;
	    buf[size++] = toupper(c);
	}
	buf[size] = 0;
	return buf;
    }

    void changeTest() {
	string sSeq = longRndSeq();
	GISeq is(sSeq);

	MGISeq ms(is);

	int change = 0;
	size_t len = is.length();

	set<int> touched;
	for (size_t i=0; i<1000; i++) {
	    size_t loc = rand() % len;
	    size_t insDel = rand() % 3;
	    GeneticSymbols t = (GeneticSymbols) (rand() & 3);

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

	GISeq* ns = ms.commit();

	if (ns == NULL) {
	    CPPUNIT_ASSERT(false);
	    return;
	}

	//cout << "Compare: " << endl << sSeq << endl << ns->toString() << endl;

	CPPUNIT_ASSERT_EQUAL(is.length() + change, ns->length());
	CPPUNIT_ASSERT_EQUAL(sSeq.length() + change, ns->toString().length());

	CPPUNIT_ASSERT(sSeq != ns->toString());
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MutableSequenceTest );
