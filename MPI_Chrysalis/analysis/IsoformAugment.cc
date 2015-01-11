 #include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "base/FileParser.h"

#include "aligns/KmerAlignCore.h"




int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-t","fasta transcript file");
  commandArg<string> bStringCmmd("-i","inchworm file");
  commandArg<string> oStringCmmd("-o","output fasta file");
  commandLineParser P(argc,argv);
  P.SetDescription("Maps chrysalis components to a reference.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(oStringCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);
  string outFile = P.GetStringValueFor(oStringCmmd);
  
  vecDNAVector ref;
  vecDNAVector iworm;
  vecDNAVector out;
  
  cout << "Reading files..." << endl;
  ref.Read(aString);
  iworm.Read(bString);
  cout << "done!" << endl;
  //test.ReverseComplement();


  TranslateBasesToNumberExact trans;
  KmerAlignCore core;
  core.SetTranslator(&trans);
  core.SetNumTables(2);
  core.AddData(ref);
  //core.SortAll();

  int i, j;

  int min = 24;
  int k = 24;


  for (i=0; i<iworm.isize(); i++) {
    const DNAVector & d = iworm[i];
    //if  (d.isize() > 300)
    //  continue;

    vec<KmerAlignCoreRecord> matchesLeft, matchesRight, matchesMiddle;
    int l = d.isize();
    core.GetMatches(matchesMiddle, d, l/2-k/2);

    if (matchesMiddle.isize() != 0) {
      //cout << "Middle match." << endl;
      continue;
    }

    core.GetMatches(matchesLeft, d, 0);
    core.GetMatches(matchesRight, d, l-k);

    //cout << "left=" << matchesLeft.isize() << " right=" << matchesRight.isize() << endl;

    if (matchesLeft.isize() == 0 || matchesRight.isize() == 0)
      continue;
    if (matchesLeft.isize() > 8 || matchesRight.isize() > 8)
      continue;

    bool bBad = false;

    svec<int> matches;
    svec<int> posLeft;
    svec<int> posRight;

    int diff = d.isize()-k;

    for (j=0; j<matchesLeft.isize(); j++) {
      int c = matchesLeft[j].GetContig();
      int p = matchesLeft[j].GetPosition();
      for (int x = 0; x<matchesRight.isize(); x++) {
	int cR = matchesRight[x].GetContig();
	int pR = matchesRight[x].GetPosition();
	if (c != cR)
	  continue;

	//if (iworm.Name(i) == ref.Name(c))
	//continue;

	int cDiff = pR-p;
	int dist = cDiff - diff;
	if (dist < 0)
	  dist = -dist;
	

	if (dist < min) {
	  bBad = true;
	  break;
	} else {
	  //cout << "Found dist=" << dist << endl;
	  matches.push_back(c);
	  posLeft.push_back(p);
	  posRight.push_back(pR);
	}

      }
    }
    if (bBad)
      continue;
    if (matches.isize() == 0)
      continue;
    
    for (j=0; j<matches.isize(); j++) {
      //if (j == 1)
      //break;

      cout << ref.Name(matches[j]) << "_alt" << endl;
      const DNAVector & orig = ref[matches[j]];
      DNAVector nn;
      nn.resize(orig.isize() + d.isize() - 2*k);
      int c = 0;
      for (int y=0; y<posLeft[j]; y++) {
	cout << orig[y];
	//nn[c] = orig[y];
	//c++;
      }
      for (int y=0; y<d.isize(); y++) {
	cout << d[y];
	//nn[c] = d[y];
	//c++;
      }
      for (int y=posRight[j]+k; y<orig.isize(); y++) {
	cout << orig[y];
	//nn[c] = orig[y];
	//c++;
      }
      cout << endl;
      //out.push_back(nn, ref.Name(matches[j])+"_alt");
    }
  }

  out.Write(outFile);

  return 0;

}
  
