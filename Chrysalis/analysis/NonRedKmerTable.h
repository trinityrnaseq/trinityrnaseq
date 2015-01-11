#ifndef NONREDKMERTABLE_H
#define NONREDKMERTABLE_H


#include "analysis/DNAVector.h"


class NonRedKmerTable
{
 public:
  NonRedKmerTable(int k) {
    m_k = k;
  }

  // Restricts all k-mers to what's in here.
  void SetUp(const vecDNAVector & templ, bool noNs = false);

  void AddData(const vecDNAVector & d);

  void AddData(vecDNAVectorStream & d);

  //ML: openmp parallel version using DNAStringStreamFast instead vecDNAVectorStream
  void AddData(DNAStringStreamFast & d);
 
  bool IsPresent(const DNAVector & d, int pos) {
    
    int i = Index(d, pos);
    //cout << s << " " << i << endl;
    
    return (i >= 0);
  }

  void SetCount(const DNAVector & d, int pos, int count) {
    int i = Index(d, pos);
    if (i < 0)
      return;
    m_counts[i] = count;
  }

  int GetCount(const DNAVector & d, int pos) {
    int i = Index(d, pos);;
    if (i < 0)
      return 0;
    return m_counts[i];
  }

  int GetCountReal(const DNAVector & d, int pos) {
    int i = Index(d, pos);;
    if (i < 0)
      return -1;
    return m_counts[i];
  }

  int GetCountReal(const string & seq, int pos) {
    string s = seq.substr(pos, m_k);
    int i = BinSearch(m_data, s);
    if (i < 0)
      return -1;
    return m_counts[i];
  }

  void SetAllCounts(int c) {
    for (int i=0; i<m_counts.isize(); i++)
      m_counts[i] = c;
  }
  
 private:
  int Index(const DNAVector & d, int pos) {
    string s; 
    d.Substring(s, pos, m_k);
    
    int i = BinSearch(m_data, s);
    return i;
  }

  int m_k;
  svec<string> m_data;
  svec<int> m_counts;


};



#endif //NONREDKMERTABLE_H

