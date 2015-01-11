#ifndef KMERTABLE_H
#define KMERTABLE_H




 
//class KmerSequence;


//KmerSequence * GetGlobalSeq();

const string defaultKmerDBName = "kmers.db";

inline void TempFile(string & out, const string & in) 
{
  char tmp[1024];
  strcpy(tmp, in.c_str());

  for (int i=0; i>=0; i--) {
    if (tmp[i] == '/') {
      tmp[i+1] = 0;
      break;
    }
  }
  out = tmp;
  out += defaultKmerDBName;
}


class KmerEntry
{
public:
  KmerEntry() {
    m_index = -1;
    m_pos = -1;
  }

  KmerEntry(int i, int pos) {
    m_index = i;
    m_pos = pos;
  }


  void SetIndex(int index, int pos) {
    m_index = index;
    m_pos = pos;
  }
  
  int Index() const {return m_index;}
  int Pos() const {return m_pos;}

  bool operator < (const KmerEntry & k) const;  


private:
  int m_index;
  int m_pos;
};



class IDS
{
public:
  IDS() {
    m_id = -1;
    m_start = -1;
    m_edge = -1;
    m_ori = 1;
  }

  IDS(int id, int start, int edge) {
    m_id = id;
    m_start = start;
    m_edge = edge;
    m_ori = 1;
  }

  int ID() const {return m_id;}
  int Start() const {return m_start;}
  int Edge() const {return m_edge;}
  int Ori() const {return m_ori;}
  void SetStart(int s) {m_start = s;}

  const string & Name() const {return m_name;}

  void SetEdge(int i) {m_edge = i;}
  void SetName(const string & n) {
    m_name = n;
  }


  bool operator < (const IDS & d) const {
#ifndef NO_REVERSE_OUT   
    if (m_ori < d.m_ori)
      return true;
    if (m_ori > d.m_ori)
      return false;
#endif    
    
    if (m_id < d.m_id)
      return true;
    if (m_id == d.m_id) {
      if (m_start < d.m_start)
	return true;
    }
    return false;
  }
  


//  bool operator < (const IDS & d) const {
//    /*if (m_id < d.m_id)
//      return true;
//    if (m_id == d.m_id) {
//      if (m_start < d.m_start)/
//	return true;
//    }
//    return false;*/
//    return m_name < d.m_name;
//  }

  void SetOri(int n) {
    m_ori = n;
  }


private:
  int m_id;
  int m_start;
  int m_edge;
  int m_ori;
  string m_name;
};


inline double FracAlign(const DNAVector & d, int startSeq, const DNAVector & read, int startRead)
{
  int i;
  int fit = 0;
  int k = startRead;
  for (i=startSeq; i<d.isize(); i++) {
    if (k >= read.isize())
      break;
    if (d[i] == read[k])
      fit++;
  }
  return (double)fit/(double)read.isize();
}



class KmerSequence
{
public:
  KmerSequence(int k, const vecDNAVector * p);

  //void SetK(long long k) {m_k = k;}
  //long long K() const {return m_k;}

  void Read(const string & fileName); 
  void Write(const string & fileName); 


  long long GetBoundValue() const {return m_count;}

  void Add(const vecDNAVector & all);

  void AddRestrict(const vecDNAVector & all, const vecDNAVector & restrict);

  void Count(svec<int> & counts) {
    long long i, j;
    
    for (i=0; i<m_kmers.lsize(); i++) {
      int c = 0;
      for (j=i; j<m_kmers.lsize(); j++) {
	if (m_kmers[j] < m_kmers[i])
	  break;
	if (m_kmers[i] < m_kmers[j])
	  break;
	c++;
      }
      counts[i] = c;
      i = j-1;
    }
  }

  long long BasesToNumber(const DNAVector & d, int off);

  long long BasesToNumberCount(int & count, const DNAVector & d, int off);

  long long BasesToNumberCountPlus(svec<IDS> & ids, long long & count, const DNAVector & d, int edge);    

  int GetK() const;

private:
  svec<KmerEntry> m_kmers;
  long long m_count;

};




//===================================================================
//========================================================================
//========================================================================






#endif 



