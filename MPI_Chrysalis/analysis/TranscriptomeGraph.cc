//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>

#include "analysis/TranscriptomeGraph.h"

//int K = 12;

class SeqChunk
{
public:
    SeqChunk() {
        m_size = 100000;
        m_counter = 0;
        // cerr << "Initial resizing m_seq to " << m_size << endl;
        m_seq.resize(m_size, 0);
    }
    
    long long Offset() const {return m_counter;}
    
    long long lsize() const {return m_size;}
    
    const char * Seq(long long i) const {
        if (i >= m_seq.isize()) {
            cerr << "Error! i=" << i << " size=" << m_seq.isize() << endl;
            return NULL;
        }
        return &m_seq[i];
    }
    
    bool Add(long long & from, long long & to, const string & s) {
        //cout << "Adding " << s << endl;
        //cout << "m_counter=" << m_counter << endl;
        
        long long len = strlen(s.c_str());
        if (len + m_counter >= m_size) {
            // cerr << "Overflow: " << len + m_counter << endl;
            return false;
        }
        
        from = m_counter;
        
        if (m_seq.lsize() == 0) {
            //cerr << "Resizing m_seq to " << m_size << endl;
            m_seq.resize(m_size, 0);
        }
        
        for (long long i=0; i<len; i++) {
            m_seq[m_counter] = s[i];
            m_counter++;
        }
        
        to = m_counter;
        
        return true;
    }
    

private:
    svec<char> m_seq; 
    long long m_counter;
    long long m_size;
};

class KmerSequence;


KmerSequence * GetGlobalSeq();

class KmerEntry
{
public:
    KmerEntry() {m_index = -1;}
    KmerEntry(long long i) {m_index = i;}
    
    
    void SetIndex(long long index) {
        m_index = index;
    }
    
    long long Index() const {return m_index;}
    
    bool operator < (const KmerEntry & k) const;  
    
private:
    long long m_index;
    
};



class KmerSequence
{
public:
    KmerSequence() {
        m_k = 24;
        m_seq.resize(256);
        m_count = 0;
        m_kmerIndex = 0;
    }
    
    void Clear() {
        m_count = 0;
        m_kmerIndex = 0;
        m_seq.clear();
        m_seq.resize(256);
        m_kmers.clear();
        m_current = "";  
    }
    
    
    void SetK(long long k) {m_k = k;}
    long long K() const {return m_k;}
    
    long long GetBoundValue() const {return m_kmerIndex+2;}
    
    const char * Sequence(long long i) const {
        if (i == -1) {
            //cout << "Called 'Sequence', returning " << m_current << endl;
            return m_current.c_str();
        }
        long long index = i / m_seq[0].lsize();
        long long i2 = i - index * m_seq[0].lsize();
        
        //cout << "index=" << index << " offset=" << i2 << endl;
        
        if (index > m_seq.isize()) {
            cerr << "ERROR! " << index << " " << m_seq.isize() << endl;
        }
        const SeqChunk & c = m_seq[index];
        const char * pRet = c.Seq(i2);
        if (pRet == NULL) {
            cerr << "Err: index=" << index << " i2=" << i2 << endl;
        }
        return pRet;
        
    }
    
    
    void Add(const string & line) {
        long long from = -1;
        long long to = -1;
        
        if (!m_seq[m_count].Add(from, to, line)) {
            m_count++;
            // cerr << "Added chunk " << m_count << endl;
            m_seq[m_count].Add(from, to, line);
        }
        
        long long index = from + m_count * m_seq[m_count].lsize();
        long long len = strlen(line.c_str());
        for (long long i=0; i<=len-m_k; i++) {
            
            if (m_kmerIndex >= m_kmers.lsize()) {
                //cerr << "Resizing to " << m_kmerIndex + 100000 << endl;
                m_kmers.resize(m_kmerIndex + 100000);
                //cerr << "size=" << m_kmers.lsize() << endl;
            }
            m_kmers[m_kmerIndex].SetIndex(index + i);
            m_kmerIndex++;
        }
        
        
    }
    
    
    void Setup() {
        //cerr << "Sorting k-mers" << endl;
        long long i;
        //for (i=0; i<m_kmers.isize(); i++) {
        //cout << "kmer=" << i << "\tindex=" << m_kmers[i].Index() << endl;
        //}
        
        m_kmers.resize(m_kmerIndex);
        
        //cerr << "Final size=" << m_kmers.lsize() << endl;
        Sort(m_kmers);
        //cerr << "done" << endl;
        
        /*
          cout << "Printing" << endl;
          for (int i=0; i<m_kmers.isize(); i++) {
          const char * p = Sequence(m_kmers[i].Index());
          for (long long j=0; j<K(); j++)
          cout << p[j];
          cout << endl;
          }
          cout << "Done" << endl;
        */
    }
    
    
    
    long long BasesToNumber(const DNAVector & d, long long off) {
        long long i;
        char tmp[256];
        
        for (i=0; i<m_k; i++) {
            tmp[i] = d[i+off];
        }
        tmp[i] = 0;
        m_current = tmp;
        
        KmerEntry dummy;
        
        //cout << "Searching " << tmp << endl;
        
        long long ret = BinSearch(m_kmers, dummy);
        //cout << "Done" << endl;
        
        //if (ret >= m_kmers.isize())
        //cout << "ERROR: " << ret << endl;
        
        //if (ret < 0)
        //ret = m_kmerIndex+1;
        
        return ret;
    } 
    
    
private:
    svec<SeqChunk> m_seq;
    svec<KmerEntry> m_kmers;
    long long m_kmerIndex;
    long long m_count;
    long long m_k;
    string m_current;
};


KmerSequence * GetGlobalSeq()
{
    static KmerSequence seq;
    return &seq;
}


bool KmerEntry::operator < (const KmerEntry & k) const {
    long long i;
    
    //cout << "my index=" << m_index << " comp=" << k.Index() << endl;
    
    const char * pOne = GetGlobalSeq()->Sequence(m_index);
    const char * pTwo = GetGlobalSeq()->Sequence(k.Index());
    int n = GetGlobalSeq()->K();
    for (i=0; i<n; i++) {
        //cout << i << "\t" << pOne[i] << "\t" << pTwo[i] << endl;
        if (pOne[i]> pTwo[i])
            return false;
        if (pOne[i] < pTwo[i])
            return true;
        
    }
    return false;
}





//========================================================================
//========================================================================
//========================================================================

class KmerNode
{
public:
    KmerNode(long long id, long long prev, long long mult, const DNAVector & d) {
        m_id = id;
        m_prev = prev;
        m_kmer = d;
        m_mult = mult;
    }
    
    const DNAVector & Kmer() const {return m_kmer;}
    long long ID() const {return m_id;}
    long long Prev() const {return m_prev;}
    long long Mult() const {return m_mult;}
    
private:
    DNAVector m_kmer;
    long long m_id;
    long long m_prev;
    long long m_mult;
    
};


class KmerGraph
{
public:
    KmerGraph() {
        m_high = 0;
    }
    
    long long Add(long long id, long long last, long long supp, const DNAVector & d) {
        //cout << "Adding node " << id << " last=" << last << endl;
        m_nodes.push_back(KmerNode(id, last, supp, d));
        return m_high;
    }
    long long Add(long long last, long long supp, const DNAVector & d) {    
        //cout << "Adding node last=" << last << endl;
        m_nodes.push_back(KmerNode(m_high, last, supp, d));
        m_high++;
        return m_high-1;
    }
    
    
    void Print(FILE * pOut, const svec<long long> & ids) const {
        if (m_nodes.lsize() < 2)
            return;
        //cout << "Prlong longing graph (verbose) " << endl;
        //fprlong longf(pOut, "Partial graph\n");
        long long i, j;
        bool bCont = false;
        
        if (m_nodes.lsize() >= 10 || ids.lsize() > 1) {
            bCont = true;
        }
        
        if (!bCont)
            return;
        
        for (i=0; i<m_nodes.lsize(); i++) {
            //cout << m_nodes[i].ID() << "\t" << m_nodes[i].Prev() << "\t" << m_nodes[i].Mult() << "\t";
            fprintf(pOut, "%d\t%d\t%d\t", (int)m_nodes[i].ID(), (int)m_nodes[i].Prev(),  (int)m_nodes[i].Mult());
            for (j=0; j<m_nodes[i].Kmer().lsize(); j++) {
                //cout << (m_nodes[i].Kmer())[j];
                fprintf(pOut, "%c", (m_nodes[i].Kmer())[j]);
            }
            for (j=0; j<ids.lsize(); j++) {
                fprintf(pOut, "\t%d", (int)ids[j]);
            }
            fprintf(pOut, "\n");
            //cout << endl;
        }
        return;
    }
    
    bool Valid() const {
        return (m_nodes.lsize() >= 10);
    }
    
    
    void Clear() {
        m_nodes.clear();
        //m_high = 0;
    }
    
private:
    svec<KmerNode> m_nodes;
    long long m_high;
};




class KmerSearch
{
public:
    KmerSearch(long long k, bool connect) {
        m_k = k;
        //TranslateBasesToNumberExact trans;
        
        //trans.SetSize(k);
        
        long long m = GetGlobalSeq()->GetBoundValue();
        m_used.resize(m, 0);
        m_usedLocal.resize(m, 0);
        m_graphID.resize(m, -1);
        
        m_connect = connect;
        if (!connect)
            m_local.reserve(100000);
        
        m_counter = 0;
    }
    
    
    
    void Extend(long long last, 
                DNAVector & seq, 
                const svec<long long> & count,
                DNAVector & longest);
    
    void Clear();
    
    long long Used(long long i) const {return m_used[i];}
    void SetUsed(long long i) {m_used[i]++;}
    void ShiftLeft(vecDNAVector & out, const DNAVector & in);
    
    void Print(FILE * p) {
        m_branch.push_back(m_counter);
        //cout << "Before: " << m_branch.isize() << endl;
        UniqueSort(m_branch);
        //cout << "After: " << m_branch.isize() << endl;
        m_graph.Print(p, m_branch);
    }
    
private:
    void Shift(vecDNAVector & out, const DNAVector & in);
    void Right(DNAVector & out, const DNAVector & in);
    
    
    
    long long m_k;
    svec<long long> m_used;
    svec<long long> m_usedLocal;
    svec<long long> m_local;
    svec<long long> m_graphID;
    svec<long long> m_branch;
    
    KmerGraph m_graph;
    bool m_connect;
    
    long long m_counter;
    bool m_bPrinted;
    
};




bool Irregular(char l)
{
    if (l == 'A' || l == 'C' || l == 'G' || l == 'T')
        return false;
    //cout << "Irregular char: " << l << endl;
    return true;
}


void KmerSearch::Clear()
{
    if (m_graph.Valid())
        m_counter++;
    
    long long i;
    
    if (!m_connect) {
        for (i=0; i<m_local.lsize(); i++) {
            m_usedLocal[m_local[i]] = 0;
        }
        
        m_local.clear();
        m_local.reserve(100000);
    }
    m_graph.Clear();
    m_branch.clear();
}

void KmerSearch::Right(DNAVector & out, const DNAVector & in)
{
    long long plus = in.lsize() - m_k;
    out.resize(m_k);
    
    for (long long j=0; j<m_k; j++) {
        out[j]=in[j+plus];
    }
    
}

void KmerSearch::Shift(vecDNAVector & out, const DNAVector & in)
{
    long long plus = in.lsize() - m_k;
    for (long long i=0; i<5; i++) {
        DNAVector tmp;
        tmp.resize(m_k);
        for (long long j=1; j<m_k; j++) {
            tmp[j-1]=in[j+plus];
        }
        tmp[tmp.lsize()-1] = NucLetter(i);
        out.push_back(tmp);
    }
}

void KmerSearch::ShiftLeft(vecDNAVector & out, const DNAVector & in)
{
    
    for (long long i=0; i<5; i++) {
        DNAVector tmp;
        tmp.resize(m_k);
        for (long long j=1; j<m_k; j++) {
            tmp[j]=in[j-1];
        }
        tmp[0] = NucLetter(i);
        out.push_back(tmp);
    }
}


long long GetFullCount(const svec<long long> & count, const DNAVector & r, int from)
{
    long long cc = 0;
    long long num = GetGlobalSeq()->BasesToNumber(r, from);
    if (num >= 0)
        cc += count[num];
    //DNAVector tmp;
    //tmp.SetToSubOf(r, from, GetGlobalSeq()->K());
    //tmp.ReverseComplement();
    
    //num = GetGlobalSeq()->BasesToNumber(tmp, 0);
    //if (num >= 0)
    //cc += count[num];
    
    return cc;
}

void KmerSearch::Extend(long long last,
                        DNAVector & seq, // kmer
                        const svec<long long> & count, // counts of kmers
                        DNAVector & longest)
{ 
    
    //TranslateBasesToNumberExact trans;
    //trans.SetSize(m_k);
    DNAVector curr;
    Right(curr, seq);
    
    long long currNum = GetGlobalSeq()->BasesToNumber(curr, 0);
    if (currNum < 0)
        return;
    long long support = count[currNum];
    
    //if (support < 3)
    //return;
    
    last = m_graph.Add(last, support, curr);
    
    m_used[currNum]++;
    m_graphID[currNum] = m_counter;
    m_usedLocal[currNum] = last;
    if (!m_connect)
        m_local.push_back(currNum);
    
    
    long long j;
    
    vecDNAVector right;
    Shift(right, seq);
    
    static int plusminus = 0;
    
    //cout << "Entering Extend" << endl;
    
    /*cout << "+";
      plusminus++;
      if (plusminus > 80) {
      cout << endl;
      plusminus = 0;
      }*/
    
    long long valid = 0;
    
    long long hi = 0;
    long long hi_index = -1;
    for (j=0; j<right.lsize(); j++) {
        const DNAVector & r = right[j];
        long long num = GetGlobalSeq()->BasesToNumber(r, 0);
        //long long num = GetFullCount(count, r, 0);
        long long cc = 0;
        if (num >= 0)
            cc = count[num];
        if (cc > hi) {
            hi = cc;
            hi_index = j;
        }
    }
    
    if (hi_index > 0) {
        //cout << "!!!!!" << endl;
        DNAVector tmp = right[0];
        right[0] = right[hi_index];
        right[hi_index] = tmp;
    }
    
    
    for (j=0; j<right.lsize(); j++) {
        const DNAVector & r = right[j];
        //cout << "r size=" << r.lsize() << endl;
        long long num = GetGlobalSeq()->BasesToNumber(r, 0);
        if (num < 0)
            continue;
        
        if (m_usedLocal[num] > 0) {
            long long mergePos = m_usedLocal[num];
            long long supp = count[num];
            long long graphID = m_graphID[num];
            m_branch.push_back(graphID);
            //cout << "Merge path, seq len=" << seq.isize() << " merge pos=" << mergePos << endl;
            m_graph.Add(mergePos, last, supp, r);
            continue;
        }
        
        if (count[num] > 0 && m_used[num] == 0 && seq.isize() < 400000) {
            
            seq.resize(seq.lsize()+1);
            seq[seq.lsize()-1] = r[r.lsize()-1];
            //cout << "Call extend." << endl;
            Extend(last, seq, count, longest);
            //cout << "Leave extend." << endl;
            seq.resize(seq.lsize()-1);
            valid++;
        }
        
    }
    if (seq.lsize() > longest.lsize()) {
        //cout << "Got seq, l=" << seq.isize()  << endl;
        longest = seq;
    }
    
    /*
      cout << "-";
      plusminus++;
      if (plusminus > 80) {
      cout << endl;
      plusminus = 0;
      }
    */
    
    
    /*
      if (valid == 0 && seq.isize() > 100) {
      cout << ">Hypothesis"  << endl;
      for (j=0; j<seq.isize(); j++) {
      if (j>0 && j % 80 == 0)
      cout << endl;
      cout << seq[j];
      }
      cout << endl;
      }*/
}


//====================================================
int TranscriptomeGraph(vecDNAVector & seq,
                       FILE * pOut,
                       int k,
                       bool connect)
{
    
    bool bAppend = true;
    
    //cerr << "Sequences: " << seq.isize() << endl;
    
    //vecbasevector contigBases;
    
    //vecDNAVector seq;
    //seq.Read(aString);
    
    int i;
    
    if (seq.isize() == 1) {    // only one sequence, the graph is simple:  linear set of overlapping kmers
        /*FILE * pOut = NULL;
          
        if (bAppend)
        pOut = fopen(oString.c_str(), "a");
        else
        pOut = fopen(oString.c_str(), "w");
        */
        
        DNAVector & d = seq[0];
        for (i=0; i<=d.isize()-k; i++) {
            fprintf(pOut, "%d\t%d\t1\t", i, i-1);
            //cout << i << "\t" << i-1 << "\t1\t";
            for (int x=i; x<i+k; x++)
                fprintf(pOut, "%c", d[x]);
            //cout << d[x];
            fprintf(pOut, "\t0\n");
            //cout << "\t0" << endl;
        }
        
        //fclose(pOut);
        return 0;
    }
    
    
    //TranslateBasesToNumberExact trans;
    
    long long j;
    // K = k;
    //cout << "Using k=" << k << endl;
    GetGlobalSeq()->Clear();
    GetGlobalSeq()->SetK(k);
    
    
    
    for (i=0; i<seq.lsize(); i++) {
        DNAVector d = seq[i];
        //cout << "i=" << i << endl;
        char * tmp = new char[d.lsize()+1];
        for (j=0; j<d.lsize(); j++) {
            //cout << j << "\t" << d[j] << endl;
            tmp[j] = d[j];
        }
        tmp[j] = 0;
        
        GetGlobalSeq()->Add(tmp);
        delete [] tmp;
        
        // why is below disabled?  Do we need it if not strand-specific data? (or do we accommodate this elsewhere?)
        /*
          d.ReverseComplement();
          for (j=0; j<d.lsize(); j++) {
          //cout << j << "\t" << d[j] << endl;
          tmp[j] = d[j];
          }
          tmp[j] = 0;
          
          GetGlobalSeq()->Add(tmp);
        */
    }
    // trans.SetSize(k);
    GetGlobalSeq()->Setup();
    
    
    long long m = GetGlobalSeq()->GetBoundValue();
    //cout << "Bound value: " << m << endl;
    
    
    svec<long long> counts;
    counts.resize(m, 0);
    
    //cout << "WARNING: no check for left seed extensions!!" << endl;
    
        
    for (i=0; i<seq.lsize(); i++) {
        DNAVector d = seq[i];
        //cout << i << endl;
        for (j=0; j<=d.lsize()-k; j++) {
            //cout << "j=" << j << endl;
            long long num = GetGlobalSeq()->BasesToNumber(d, j);
            if (num < 0)
                continue;
            //cout << "num=" << num << endl;
            counts[num]++;
        }
        /*
          d.ReverseComplement();
          for (j=0; j<=d.lsize()-k; j++) {
          //cout << "j=" << j << endl;
          long long num = GetGlobalSeq()->BasesToNumber(d, j);
          //cout << "num=" << num << endl;
          counts[num]++;
          }*/
    }    
    
    //cerr << "Done here." << endl;
    
    KmerSearch search(k, connect);
    
    /*FILE * pOut = NULL;
      
    if (bAppend)
    pOut = fopen(oString.c_str(), "a");
    else
    pOut = fopen(oString.c_str(), "w");
    */
    
    
    // examine kmers
    for (i=0; i<seq.lsize(); i++) {
        const DNAVector & v = seq[i];
        
        for (j=0; j<=v.lsize()-k; j++) {
            DNAVector d; // kmer
            d.SetToSubOf(v, j, k);
            
            long long num2 = GetGlobalSeq()->BasesToNumber(d, 0);
            
            //cerr << "iworm(" << i << "," << j << ") = " << d.AsString() << " set to num2: " << num2 << endl;
            

            if (num2 < 0 || search.Used(num2) > 0)
                continue;
            /*
              vecDNAVector left;
              search.ShiftLeft(left, d);
              
              bool bLeft = false;
              for (long long y=0; y<left.isize(); y++) {
              long long num2 = trans.BasesToNumber(left[y], 0);
              if (counts[num2] > 0 && search.Used(num2) == 0)
              bLeft = true;
              }
              if (bLeft)
              continue;
            */
            
            
            search.Clear();
            
            //cout << "Trying..." << endl;
            
            DNAVector longest;
            search.Extend(-1, d, counts, longest);     
            
            /*
              if (longest.lsize() > 100) {
              cout << ">Sequence_" << i << endl;
              for (long long x=0; x<longest.lsize(); x++) {
              if (x>0 && x % 80 == 0)
              cout << endl;
              cout << longest[x];
              }
              cout << endl;
              }*/
            
            //cout << "Printing." << endl;
            search.Print(pOut);
            search.SetUsed(num2);
        }
    }
    
    
    //cout << "All done!" << endl;
    //fclose(pOut);
    
    return 0;
    
}
