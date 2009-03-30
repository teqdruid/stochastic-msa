#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <ctype.h>
#include "storage.h"

class ImmutableSequenceTest: public CppUnit::TestFixture  {

    CPPUNIT_TEST_SUITE( ImmutableSequenceTest );

    CPPUNIT_TEST( testParse );

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
};


CPPUNIT_TEST_SUITE_REGISTRATION( ImmutableSequenceTest );
