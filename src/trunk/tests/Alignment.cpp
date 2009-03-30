#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <ctype.h>
#include "storage.h"
#include "algo.h"

class PairwiseAlignmentTest: public CppUnit::TestFixture  {

    CPPUNIT_TEST_SUITE( PairwiseAlignmentTest );

    CPPUNIT_TEST( identAlignment );

    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {

    }

    void tearDown() {

    }

    void identAlignment() {
	const double alpha = .1;
	const double diag = 1.0 - alpha * 4;
	GenScores sm(alpha, -.3);
	
	GISeq seq("AAAAAAAA");

	double score = alignmentScore(sm, seq, seq);

	CPPUNIT_ASSERT_EQUAL(diag * seq.length(), score);
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION( PairwiseAlignmentTest );
