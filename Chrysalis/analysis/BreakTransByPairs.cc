#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>
#include "base/FileParser.h"
#include "analysis/KmerTable.h"


  

//===================================================================
//========================================================================
//========================================================================


bool Irregular(char l)
{
  if (l == 'A' || l == 'C' || l == 'G' || l == 'T')
    return false;
  //cout << "Irregular char: " << l << endl;
  return true;
}


int FindPartner(svec<IDS> & ids, const string & name) 
{
  char tmp[256];
  strcpy(tmp, name.c_str());
  tmp[strlen(tmp)-1] = '1';

  IDS dummy;
  dummy.SetName(tmp);

  int index = BinSearch(ids, dummy);
  return index;
}


int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-i","read fasta file");
  commandArg<string> fStringCmmd("-f","fasta file");
  commandArg<string> oStringCmmd("-o","fasta output");
  commandArg<int> kCmmd("-k","kmer size", 24);
  commandArg<bool> strandCmmd("-doublestrand","not strand-specific", false);
  //commandArg<bool> cCmmd("-nc","do not fully connected graph", false);

  commandLineParser P(argc,argv);
  P.SetDescription("Assembles k-mer sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(fStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(kCmmd);
  P.registerArg(strandCmmd);
  //P.registerArg(cCmmd);


  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string fString = P.GetStringValueFor(fStringCmmd);
  string oString = P.GetStringValueFor(oStringCmmd);
  int k = P.GetIntValueFor(kCmmd);
  bool bStrand = P.GetBoolValueFor(strandCmmd);




  int i, j;
  //vecbasevector contigBases;
 
  vecDNAVector dna;
  dna.Read(fString);

  vecDNAVector seq;
  seq.Read(aString);

  
  KmerSequence kmers(k, &seq);;
  kmers.Add(seq);
  
  long long m = kmers.GetBoundValue();


  vecDNAVector out;
  //dna.Read(fString);

  svec<int> used;
  used.resize(seq.size(), 0);

  int offset = seq[0].isize();

  int broken = 0;

  int extra = 25;

  for (i=0; i<dna.size(); i++) {
    DNAVector one = dna[i];
    double all = 0.;
    broken = 0;
    //cout << "Sequence length: " << one.isize() << endl;
    svec<IDS> reads;
    for (j=0; j<=one.isize()-k; j++) {
      DNAVector sub;
      sub.SetToSubOf(one, j, k);
      long long n1, n2;
      
      svec<IDS> ids;
     
      int edge = 0;
      kmers.BasesToNumberCountPlus(ids, n1, sub, edge);

      for (int x=0; x<ids.isize(); x++) {
	if (used[ids[x].ID()] > 0)
	  continue;

	double ff = FracAlign(one, j, seq[ids[x].ID()], ids[x].Start());
	if (ff > 0.95) {
	  used[ids[x].ID()] = 1;
	  ids[x].SetEdge(j);
	  ids[x].SetName(seq.Name(ids[x].ID()));
	  reads.push_back(ids[x]);
	}
      }


      if (bStrand) {
	sub.ReverseComplement();
	kmers.BasesToNumberCountPlus(ids, n2, sub, edge);
      }
      all += n1 + n2;
    }

    Sort(reads);

    for (j=0; j<reads.isize(); j++) {
      const char * p = reads[j].Name().c_str();
      //cout << "Considering " << p << endl;
      if (p[strlen(p)-1] != '2')
	continue;
      int idx = FindPartner(reads, p);
      if (idx == -1)
	continue;
      //cout << "Not found." << endl;
      //cout << "Found partner for " << p << endl;

      offset = seq[j].isize();
      int start = reads[j].Edge()+offset+extra;
      int end = reads[idx].Edge()-extra;
      //cout << "Span from " << start-offset << " to " << end << "\t" << end - start + offset << "\t" << p << "\t" << reads[idx].Name() << endl;
      for (int x = start; x<end; x++) {
	if (one.Qual(x) < 100)
	  one.SetQual(x, one.Qual(x)+1);
      }

    }
    int baldStart = -2;
    int lastBaldStart = 0;
    char tmp[256];
    //cout << "Printing" << endl;
    for (int x=0; x<one.isize(); x++) {
      //cout << x << "\t" << one.Qual(x) << endl;
      if (one.Qual(x) > 0) {
	if (baldStart == -2)
	  baldStart = -1;
	//cout << "baldStart=" << baldStart << endl;
	if (baldStart >= 0) {
	  int middle = (x+baldStart)/2;
	  
	  sprintf(tmp, "%s_part_%d", dna.Name(i).c_str(), broken);
	  broken++;
	  DNAVector piece;
	  cout << "Found Break: start=" << lastBaldStart << " to=" <<  middle-lastBaldStart+k/2 << " middle=" << middle << endl;
	  piece.SetToSubOf(one, lastBaldStart, middle-lastBaldStart+k/2);
	  out.push_back(piece, tmp);
	  baldStart = -1;
	  lastBaldStart = middle-k/2;
	}
      } else {
	if (baldStart == -1)
	  baldStart = x;
      }
    }
    sprintf(tmp, "%s_part_%d", dna.Name(i).c_str(), broken);
    broken++;
    DNAVector piece2;
    piece2.SetToSubOf(one, lastBaldStart, one.isize()-lastBaldStart);
    out.push_back(piece2, tmp);
    
    

    /*
    double cov = all/(double)one.isize();
    char tmp[128];
    sprintf(tmp, "%f", cov); 
    string name = dna.Name(i);
    name += "_cov_";
    name += tmp;
    dna.SetName(i, name);
    */
  }

  out.Write(oString);
 
  return 0;

}
