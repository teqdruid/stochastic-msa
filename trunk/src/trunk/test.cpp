//#include <cppunit/extensions/TestFactoryRegistry.h>
//#include <cppunit/ui/text/TestRunner.h>
#include "storage.h"
#include "pairwisealign.cpp"
/*int main( int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );
  bool wasSuccessful = runner.run( "", false );
  return wasSuccessful;
}
*/

int main()
{
	GenScores* scoreMatrix = new GenScores(0.3, 0.02);
	GISeq* seqA = new GISeq("ACGT");
	GISeq* seqB = new GISeq("GCAT");
	const char*** directions;
	double score = alignmentScore <GenScores, GISeq, GISeq> (*scoreMatrix, *seqA, *seqB);
	directions = getAlignment <GenScores, GISeq, GISeq> (*scoreMatrix, *seqA, *seqB);
	reconstructAlignment(*seqA, *seqB, directions);
	cout << "Maximum score is " << score << "\n";
	for(size_t i = 0; i < a.length() + 1; i++)
	{
		delete[] directions[i];
	}
	delete[] directions;
}
