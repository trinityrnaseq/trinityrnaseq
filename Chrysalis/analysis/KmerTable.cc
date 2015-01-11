#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include "analysis/KmerTable.h"
#include "util/mutil.h"
#include "analysis/NonRedKmerTable.h" 


int K = 12;
const vecDNAVector * pAll = NULL;
DNAVector dummy;


KmerSequence::KmerSequence(int k, const vecDNAVector * p) 
{
  m_count = 0;
  K = k;
  pAll = p;
}

int KmerSequence::GetK() const
{
  return K;
}


long long KmerSequence::BasesToNumber(const DNAVector & d, int off) {
  long long i;
  
  dummy = d;
  KmerEntry dummyKmer;
  dummyKmer.SetIndex(-1, off);
  //cout << "Searching " << tmp << endl;
  
  long long ret = BinSearch(m_kmers, dummyKmer);
  
  return ret;
} 

long long KmerSequence::BasesToNumberCount(int & count, const DNAVector & d, int off) {
  long long i;
  
  dummy = d;
  KmerEntry dummyKmer;
  dummyKmer.SetIndex(-1, off);
  //cout << "Searching " << tmp << endl;
  
  
  long long ret = BinSearch(m_kmers, dummyKmer);
  //cout << "ret=" << ret << endl;
  
  count = 0;
  if (ret < 0)
    return -1;
  
  for (i=ret; i<m_kmers.lsize(); i++) {
    //cout << "ret=" << ret << " i=" << i << endl;
    if (m_kmers[i] < m_kmers[ret])
      break;
    if (m_kmers[ret] < m_kmers[i])
      break;
    count++;
  }
  
  return ret;
} 

long long KmerSequence::BasesToNumberCountPlus(svec<IDS> & ids, long long & count, const DNAVector & d, int edge) {
  long long i;
  
  dummy = d;
  KmerEntry dummyKmer;
  dummyKmer.SetIndex(-1, 0);
  //cout << "Searching " << tmp << endl;
  
  //FIXME
  //------------------------------------------------------------
  // note, this is very confusing as written. global 'dummy', as assigned to parameter kmer 'd' is actually being searched. 
  //------------------------------------------------------------
  
  

  long long ret = BinSearch(m_kmers, dummyKmer);
  //cerr << "ret=" << ret << endl;
  
  count = 0;
  if (ret < 0) {
      //cerr << "Searched: " << d.AsString() << "\tMISSING" << endl;
      return -1;
  }
  for (i=ret; i<m_kmers.lsize(); i++) {
      //cerr << "ret=" << ret << " i=" << i << endl;
    if (m_kmers[i] < m_kmers[ret])
      break;
    if (m_kmers[ret] < m_kmers[i])
      break;
    ids.push_back(IDS(m_kmers[i].Index(), m_kmers[i].Pos(), edge));
    count++;
  }
  

  // cerr << "Searched: " << d.AsString() << "\tcount: " << count << endl;
  
  return ret;
} 






void KmerSequence::Add(const vecDNAVector & all) 
{
  pAll = &all;

  //long long est = (long long)(all[0].isize()-K+1)*(long long)all.isize();

  long long est = 0;
  for (int i=0; i<all.isize(); i++) {
    est += (long long)(all[i].isize()-K+1);
  }
  // cout << "Estimate space: " << est << endl;


  m_kmers.resize(est);
  m_count = 0;
  for (int j=0; j<all.isize(); j++) {
    const DNAVector & d = all[j];
    for (int i=0; i<=d.isize()-K; i++) {
      if (m_kmers.lsize() <= m_count)
	m_kmers.resize(m_count + 1000000);
      m_kmers[m_count].SetIndex(j, i);
      m_count++;
    }
  }
  // cout << "Effective space: " << m_count << endl;
  m_kmers.resize(m_count);
  // cout << "Sorting..." << endl;
  Sort(m_kmers);
  // cout << "done!" << endl;
  
}

void KmerSequence::AddRestrict(const vecDNAVector & all, const vecDNAVector & restrict) 
{
  pAll = &all;


  cout << "Restriction table..." << endl;
  NonRedKmerTable nr(K);
  nr.SetUp(restrict);


  long long est = (long long)(all[0].isize()-K+1)*(long long)all.isize() / 4;
  cout << "Estimate space: " << est << " (" << est * 4 << ")" << endl;
  m_kmers.resize(est);
  m_count = 0;
  int kk = 0;
  for (int j=0; j<all.isize(); j++) {
    const DNAVector & d = all[j];
    for (int i=0; i<=d.isize()-K; i++) {
      if (m_kmers.lsize() <= m_count)
	m_kmers.resize(m_count + 1000000);

      kk++;
      if (kk % 1000000 == 0) {
	cout << kk << " (" << m_count << ")" << endl;
      }

      if (!nr.IsPresent(d, i))
	continue;

      m_kmers[m_count].SetIndex(j, i);
      m_count++;
    }
  }
  cout << "Effective space: " << m_count << endl;
  m_kmers.resize(m_count);
  cout << "Sorting..." << endl;
  Sort(m_kmers);
  cout << "done!" << endl;
  
}



bool KmerEntry::operator < (const KmerEntry & k) const {
  int i;

  //const DNAVector & me = (*pAll)[m_index];
  const DNAVector * me = NULL;
  if (m_index < 0) {
    //cout << "Using dummy." << endl;
    me = &dummy;
  } else {
    me = &((*pAll)[m_index]);
  }
  
  const DNAVector * you = NULL;
  if (k.Index() < 0) {
    //cout << "Using dummy." << endl;
    you = &dummy;
  } else {
    you = &((*pAll)[k.Index()]);
  }

  //ML: perform the actual comparison on the raw char array
  const char* me_str  = ( &((*me)[0])  ) + m_pos;
  const char* you_str = ( &((*you)[0]) ) + k.Pos();
  for (i=0; i<K; i++) {
    if (me_str[i] > you_str[i])
      return false;
    if (me_str[i] < you_str[i])
      return true;
  }
  return false;
}
  

void KmerSequence::Read(const string & fileName)
{

  int ver;
  CMReadFileStream f;
  f.Open(fileName.c_str());
  long long s;

  f.Read(ver);
  f.Read(s);
  f.Read(K);

  m_kmers.resize(s);
  for (int i=0; i<s; i++) {
    int index, pos;
    f.Read(index);
    f.Read(pos);
    
    m_kmers[i].SetIndex(index, pos);
  }

}

void KmerSequence::Write(const string & fileName)
{
  int ver = 1;
  
  CMWriteFileStream f;
  f.Open(fileName.c_str());
  f.Write(ver);

  f.Write(m_count);
  f.Write(m_kmers.isize());
  f.Write(K);

  for (int i=0; i<m_kmers.isize(); i++) {
    f.Write(m_kmers[i].Index());
    f.Write(m_kmers[i].Pos());
  }
}
