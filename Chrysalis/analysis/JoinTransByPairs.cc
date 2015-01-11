#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>
#include "base/FileParser.h"

#include "util/mutil.h"
#include "analysis/KmerTable.h"


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

int FindPartnerRead(vecDNAVector & dna, const string & name) 
{
  char tmp[256];
  strcpy(tmp, name.c_str());
  tmp[strlen(tmp)-1] = '1';


  int index = dna.NameIndex(tmp);
  return index;
}

bool IsFW(const string & name) {
  if (name[name.size()-1] == '2') {
    //cout << "Add fw " << name << endl;
    return true;
  }
  return false;
}


class ReadPlaces
{
public:
  ReadPlaces() {
  }

  void AddFW(int i, int pos) {
    m_fw.push_back(i);
    m_fwPos.push_back(pos);
  }
  void AddRC(int i, int pos) {
    m_rc.push_back(i);
    m_rcPos.push_back(pos);
  }

  int GetFWCount() const {return m_fw.isize();}
  int GetFW(int i) const {return m_fw[i];}
  int GetFWPos(int i) const {return m_fwPos[i];}
  int GetRCCount() const {return m_rc.isize();}
  int GetRC(int i) const {return m_rc[i];}
  int GetRCPos(int i) const {return m_rcPos[i];}

  int GetRCPosByRead(int read) const {
    for (int i=0; i<m_rc.isize(); i++) {
      if (m_rc[i] == read)
	return m_rcPos[i];
    }
    cout << "Error: read not found!" << endl;
    return 10000;
  }

  void Clear() {
    m_fw.clear();
    m_rc.clear();
    m_fwPos.clear();
    m_rcPos.clear();
  }

  void Absorb(const ReadPlaces & pl, int len1, int len2) {
    int i;
    for (i=0; i<pl.GetFWCount(); i++) {
      m_fw.push_back(pl.GetFW(i));
      m_fwPos.push_back(pl.GetFWPos(i)+len2);
    }
    
    for (i=0; i<pl.GetRCCount(); i++) {
      m_rc.push_back(pl.GetRC(i));
      m_rcPos.push_back(pl.GetRCPos(i)+len1);
    }
  }

  void Write(CMWriteFileStream & f) const {
    int i;
    f.Write(m_fw.isize());
    for (i=0; i<m_fw.isize(); i++) {
      f.Write(m_fw[i]);
      f.Write(m_fwPos[i]);
    }
    f.Write(m_rc.isize());
    for (i=0; i<m_rc.isize(); i++) {
      f.Write(m_rc[i]);
      f.Write(m_rcPos[i]);
    }
  }

  void Read(CMReadFileStream & f) {
    int i, n;
    f.Read(n);
    m_fw.resize(n);
    m_fwPos.resize(n);
    for (i=0; i<m_fw.isize(); i++) {
      f.Read(m_fw[i]);
      f.Read(m_fwPos[i]);
    }
    f.Read(n);
    m_rc.resize(n);
    m_rcPos.resize(n);
    for (i=0; i<m_rc.isize(); i++) {
      f.Read(m_rc[i]);
      f.Read(m_rcPos[i]);
    }
    
  }
  
  
private:
  svec<int> m_fw;
  svec<int> m_fwPos;
  
  svec<int> m_rc;
  svec<int> m_rcPos;
};


void ReadThePlaces(svec<ReadPlaces> & plac, const string & file) {
  CMReadFileStream f;
  f.Open(file.c_str());
  int n;  
  f.Read(n);
  plac.resize(n);
  for (int i=0; i<n; i++)
    plac[i].Read(f);
  f.Close();

}

void WriteThePlaces(const svec<ReadPlaces> & plac, const string & file) {
  CMWriteFileStream f;
  f.Open(file.c_str());
  
  f.Write(plac.isize());

  for (int i=0; i<plac.isize(); i++)
    plac[i].Write(f);
  f.Close();

}

bool DevOK(double a, double a_dev, double b, double b_dev) 
{
  if ((a-a_dev)/(b+b_dev) < 2.5 && (b-b_dev)/(a+a_dev) < 2.5)
    return true;
  
  return false;
}

bool IsJoinOK(double a, double b, double c) 
{
  double a_dev = sqrt(a);
  double b_dev = sqrt(b);
  double c_dev = sqrt(c);

  cout << "Try this - contig a: mean=" << a << " +/-" << a_dev;
  cout << "  contig b: mean=" << b << " +/-" << b_dev;
  cout << "  joins: mean=" << c << " +/-" << c_dev << endl;
  
  return DevOK(a, a_dev, b, b_dev) && DevOK((a+b)/2, a_dev, c, c_dev);

}



int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-i","read fasta file");
  commandArg<string> fStringCmmd("-f","fasta file");
  commandArg<string> oStringCmmd("-o","fasta output");
  commandArg<int> kCmmd("-k","kmer size", 24);
  commandArg<bool> strandCmmd("-doublestrand","not strand-specific", false);
  commandArg<bool> padCmmd("-pad","pad witn N's if no merge", false);
  commandArg<bool> contCmmd("-cont","continue w/ previous run", false);
  //commandArg<bool> cCmmd("-nc","do not fully connected graph", false);

  commandLineParser P(argc,argv);
  P.SetDescription("Assembles k-mer sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(fStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(kCmmd);
  P.registerArg(strandCmmd);
  P.registerArg(padCmmd);
  P.registerArg(contCmmd);
  //P.registerArg(cCmmd);


  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string fString = P.GetStringValueFor(fStringCmmd);
  string oString = P.GetStringValueFor(oStringCmmd);
  int k = P.GetIntValueFor(kCmmd);
  bool bStrand = P.GetBoolValueFor(strandCmmd);

  bool bCont = P.GetBoolValueFor(contCmmd);

  bool bPad = P.GetBoolValueFor(padCmmd);


  string placementFile = fString + ".placements";

  int i, j;
  //vecbasevector contigBases;
 
  vecDNAVector seq;
  seq.Read(aString);

  vecDNAVector dna;
  dna.Read(fString);

  int readlen = seq[0].isize();
 
  KmerSequence kmers(k, &seq);

  if (!bCont) {
    kmers.AddRestrict(seq, dna);
    //kmers.Add(seq);
  } else {
    cout << "Continuing..." << endl;
  }

  long long m = kmers.GetBoundValue();


  vecDNAVector out;
  //dna.Read(fString);

  //svec<int> used;
  //used.resize(seq.isize(), 0);

  int offset = seq[0].isize();

  int broken = 0;

  int extra = 100;

  svec<int> contigs;
  contigs.resize(seq.isize(), -1);

  svec<ReadPlaces> places;

  if (!bCont) {
    places.resize(dna.isize());
  } else {
    cout << "Reading placements..." << endl;
    ReadThePlaces(places, placementFile);
    cout << "done." << endl;
  }
   
  int maxDist = 400;


  cout << "Assigning contigs to reads." << endl;
  for (i=0; i<dna.isize(); i++) {
    const DNAVector & one = dna[i];
    double all = 0.;
    broken = 0;
    //cout << "Sequence length: " << one.isize() << " seq: " << i << endl;
   
    ReadPlaces & pl = places[i];

    if (bCont) {
      for (j=0; j<pl.GetFWCount(); j++) {
	contigs[pl.GetFW(j)] = i;
      }
      for (j=0; j<pl.GetRCCount(); j++) {
	contigs[pl.GetRC(j)] = i;
      }

    } else {
      for (j=0; j<=one.isize()-k; j++) {
	DNAVector sub;
	sub.SetToSubOf(one, j, k);
	long long n1, n2;
	
	svec<IDS> ids;
	
	int edge = 0;
	kmers.BasesToNumberCountPlus(ids, n1, sub, edge);
	
	for (int x=0; x<ids.isize(); x++) {
	  if (contigs[ids[x].ID()] > 0)
	    continue;
	  
	  double ff = FracAlign(one, j, seq[ids[x].ID()], ids[x].Start());
	  if (ff > 0.98) {
	    if (IsFW(seq.Name(ids[x].ID()))) {
	      pl.AddFW(ids[x].ID(), one.isize()-j);
	    } else {
	      pl.AddRC(ids[x].ID(), j+readlen);
	    }
	    contigs[ids[x].ID()] = i;
	  }
	}
      }
    }
  }

  //if (!bCont)
  //WriteThePlaces(places, placementFile);
  
  for (i=0; i<dna.isize(); i++) {

    DNAVector & d = dna[i];
    if (d.isize() == 0)
      continue;
    if (d.isize() < 100)
      continue;

    svec<int> possible;
    possible.resize(dna.isize(), 0);

    //cout << "Checking contig " << i << endl;
    ReadPlaces & pl = places[i];

    double currCov = (double)(pl.GetFWCount() + pl.GetRCCount() * readlen)/(double)d.isize();

    int best = 0;
    int max = 0;
    for (j=0; j<pl.GetFWCount(); j++) {
      int o = pl.GetFW(j);
      int dist = pl.GetFWPos(j); 
      const string & name = seq.Name(o);
      int r = FindPartnerRead(seq, name);
      if (r == -1)
	continue;
      
      int cont = contigs[r];
      if (cont == i || cont == -1)
	continue;
      //cout << "Found possible: " << cont << endl;
      
      const ReadPlaces & plOther = places[cont];
      dist += plOther.GetRCPosByRead(r);

      if (dist > maxDist)
	continue;

      possible[cont]++;
      if (possible[cont] > max) {
	max = possible[cont];
	best = cont;
      }
    }
    if (max > 2) {
      DNAVector & second = dna[best];
      double otherCov = (double)(places[best].GetFWCount() + places[best].GetRCCount() * readlen)/(double)second.isize();

      if (IsJoinOK(currCov, otherCov, max)) {

	cout << "Joining " << dna.Name(i) << " to " << dna.Name(best) << " support: " << max << endl;
	int firstLen = d.isize();
	int secondLen = second.isize();

	bool bMerged = false;

	if (d.Append(second, 12, 25, 0.95)) {
	  cout << "Merged sequences w/o N's." << endl;
	  second.resize(0);
	  bMerged = true;
	} else {
	  if (bPad) {
	    cout << "Padding with N's" << endl;
	    DNAVector nn = d;	
	    
	    nn.resize(d.isize() + 8);
	    for (j=d.isize(); j<nn.isize(); j++)
	      nn[j] = 'N';
	    int off = nn.isize();
	    nn.resize(off+second.isize());
	    for (j=off; j<nn.isize(); j++)
	      nn[j] = second[j-off];
	    second.resize(0);
	    dna[i] = nn;
	    bMerged = true;
	  }
	}

	if (bMerged) {
	  string newName = dna.Name(i);
	  newName += "_";
	  newName += dna.Name(best);
	  dna.SetName(i, newName);
	  
	  for (j=0; j<places[best].GetFWCount(); j++) {
	    contigs[places[best].GetFW(j)] = i;
	  }
	  for (j=0; j<places[best].GetRCCount(); j++) {
	    contigs[places[best].GetRC(j)] = i;
	  }
	  
	  
	  pl.Absorb(places[best], firstLen, secondLen);
	  places[best].Clear();
	  i--;
	}
      }
    }
  }


  

  dna.Write(oString, true);
 
  return 0;

}
