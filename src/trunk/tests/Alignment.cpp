#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <ctype.h>
#include "storage.h"
#include "algo.h"

class PairwiseAlignmentTest: public CppUnit::TestFixture  {

    CPPUNIT_TEST_SUITE( PairwiseAlignmentTest );

    CPPUNIT_TEST( scoreMatrix );
    CPPUNIT_TEST( identAlignment );

    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {

    }

    void tearDown() {

    }

    void scoreMatrix() {
	GenScores sm(.8, .3, 1, .1);
	CPPUNIT_ASSERT_EQUAL(0.8, sm.score(A, A));
	CPPUNIT_ASSERT_EQUAL(0.3, sm.score(A, T));
    }

#define FEQ(A, B, EP) ((A < (B + EP)) && (A > (B - EP)))

    void identAlignment() {
	GenScores sm(.8, .3, 1, .1);
	
	GISeq seq("AACATAGAACTTAGCTTAGCTATCGCGGTTTTACGATTCGATCGATTCGAT");

	double score = alignmentScore(sm, seq, seq);

	CPPUNIT_ASSERT(FEQ(score, (.8*seq.length()), .0000001));
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION( PairwiseAlignmentTest );
