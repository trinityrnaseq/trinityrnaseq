#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>

#include "analysis/KmerTable.h"




//========================================================================
//========================================================================
//========================================================================


void Right(vecDNAVector & out, const DNAVector & in, int K)
{
  DNAVector tmp;
  tmp.resize(K);
  int len = in.isize();
  for (int i=0; i<K-1; i++) {
    tmp[i] = in[len-K+i+1];
  }
  tmp[K-1] = 'A';
  out.push_back(tmp);
  tmp[K-1] = 'C';
  out.push_back(tmp);
  tmp[K-1] = 'G';
  out.push_back(tmp);
  tmp[K-1] = 'T';
  out.push_back(tmp);

}


void Extend(KmerSequence & seq, svec<int> & counts, DNAVector & sub, bool rc, bool bStrand, int k) 
{
  //static int deep = 0;
  //deep++;
  //cout << "Deep=" << deep << endl;

  size_t i, j;
  
  vecDNAVector ex;
  Right(ex, sub, k);
  int max = 0;
  long long indexFW = -1;
  long long indexRC = -1;
  long long index = -1;

  int max2 = 0;

  bool doFwd = false;
  bool doRC = false;

  if (!bStrand) {
    doFwd = true;
    doRC = true;
  } else {
    if (rc) {
      doRC = true;
    } else {
      doFwd = true;
    }
  }


  for (i=0; i<ex.size(); i++) {
    int nFW = 0;
    int nRC = 0;
    long long iFW = -1;

    if (doFwd) {
      iFW = seq.BasesToNumber(ex[i], 0);
      if (iFW >=0)
	nFW = counts[iFW];
    }
    long long iRC = -1;
    if (doRC) {
      ex[i].ReverseComplement();
      iRC = seq.BasesToNumber(ex[i], 0);
      //ex[i].ReverseComplement();
      if (iRC >=0)
	nRC = counts[iRC];
    }

    if (nFW + nRC > max) {
      max2 = max;
      max = nFW + nRC;
      indexFW = iFW;
      indexRC = iRC;
      index = i;
    } else {
      if (nFW + nRC > max2)
	max2 = nFW + nRC;
    }
  }
  if (max < 2)
    return;

  int diff = 0; //max2/2;
  //cout << "max=" << max << " diff=" << diff << endl;

  if (indexFW >= 0)
    counts[indexFW] = diff;
  if (indexRC >= 0)
    counts[indexRC] = diff;
  
  sub.resize(sub.isize()+1);
  sub[sub.isize()-1] = NucLetter(index);
  Extend(seq, counts, sub, rc, bStrand, k);
}


bool Irregular(char l)
{
  if (l == 'A' || l == 'C' || l == 'G' || l == 'T')
    return false;
  //cout << "Irregular char: " << l << endl;
  return true;
}


int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-i","fasta file");
  commandArg<string> oStringCmmd("-o","fasta output");
  commandArg<int> kCmmd("-k","kmer size", 25);
  commandArg<bool> strandCmmd("-strand","strand-specific data", false);
  //commandArg<bool> cCmmd("-nc","do not fully connected graph", false);

  commandLineParser P(argc,argv);
  P.SetDescription("Assembles k-mer sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(kCmmd);
  P.registerArg(strandCmmd);
  //P.registerArg(cCmmd);


  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string oString = P.GetStringValueFor(oStringCmmd);
  int k = P.GetIntValueFor(kCmmd);
  bool bStrand = P.GetBoolValueFor(strandCmmd);


  size_t i, j;
  //vecbasevector contigBases;
 
  vecDNAVector seq;
  seq.Read(aString);

  /*
  if (bStrand) {
    int rev = 0;
    cout << "Reversing all /2 reads." << endl;
    for (i=0; i<seq.isize(); i++) {
      const char *p = seq.Name(i).c_str();
      int n = strlen(p);
      if (p[n-2] == '/' && p[n-1] == '2') {
	seq[i].ReverseComplement();
	rev++;
      }
    }
    
    cout << "done, reversed " << rev << " sequences" << endl;
    }*/


  
  KmerSequence kmers(k, &seq);
  kmers.Add(seq);
  
  long long m = kmers.GetBoundValue();


  svec<int> counts;
  counts.resize(m, 0);

  cout << "Counting." << endl;

  kmers.Count(counts);
  /*
  for (i=0; i<seq.lsize(); i++) {
    DNAVector d = seq[i];
    for (j=0; j<=d.lsize()-k; j++) {
     
      int num = 0;
      //cout << "j=" << j << endl;
      int n = kmers.BasesToNumberCount(num, d, j);
      //cout << "n=" << n << " num=" << num << endl;
      if (n < 0 || counts[n] > 0)
	continue;
    
      counts[n] = num;
    }  
    } */   

  cout << "Done." << endl;

  FILE * pOut = fopen(oString.c_str(), "w");
  int cc = 0;
  for (i=0; i<seq.size(); i++) {
    DNAVector d = seq[i];
    
    for (j=0; j<=d.lsize()-k; j++) {
      bool bNs = false;
      for (int x=j; x<j+k; x++) {
	if (d[x] != 'A' && d[x] != 'C' && d[x] != 'G' && d[x] != 'T') 
	  bNs = true;
      }
      if (bNs)
	continue;

      long long n = kmers.BasesToNumber(d, j);
      if (n < 0)
	continue;
      if (counts[n] == 0)
	continue;
      DNAVector sub;
      sub.SetToSubOf(d, j, k);
      
      Extend(kmers, counts, sub, false, bStrand, k);
      sub.ReverseComplement();
      Extend(kmers, counts, sub, true, bStrand, k);

      counts[n] = 0;

      if (sub.isize() >= 2*k-2) {
	//if (sub.isize() >= k) {
	sub.ReverseComplement();
	//cout << ">Sequence_" << cc << endl;
	fprintf(pOut, ">Sequence_%d\n", cc);
	cc++;
	for (int x=0; x<sub.isize(); x++) {
	  if (x>0 && x%80 == 0)
	    fprintf(pOut, "\n");
	  fprintf(pOut, "%c", sub[x]);
	    //cout << endl;
	    //cout << sub[x];
	  
	}
	fprintf(pOut, "\n");
	//cout << endl;
      }
      
    }  
  }    


  fclose(pOut);
  return 0;

}
