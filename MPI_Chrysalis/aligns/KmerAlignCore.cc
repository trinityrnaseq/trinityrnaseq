#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "aligns/KmerAlignCore.h"

KmerAlignCore::KmerAlignCore() 
{
    m_numTables = 2;    
    m_pTrans = NULL;
    m_lookAhead = 0;
    m_lookAheadMaxFreq = 50000;
    
    m_max12 = 0x7FFFFFFF;
}

void KmerAlignCore::AddData(const vecDNAVector & bases)
{
    vecNumVector dummy;
    dummy.resize(bases.size());
    AddData(bases, dummy, 1);
}


bool IsRepeat(const NumVector & t, int i, int size, int min)
{
    if (t.size() == 0)
        return false;
    
    int r = 0; 
    int n = 0; 
    for (int j=i; j<i+size; j+=4) {
        if (t[j] >= min)
            r++;
        n++;
    }
    if (r >= n / 2) {
        //cout << "Repeat, skipping!" << endl;
        return true;
    } else {
        return false;
    }
}

void KmerAlignCore::AddData(const vecDNAVector & bases, const vecNumVector & tags, int min)
{
    // First, count k-mers
    int i, j, k;
    
    int size = m_pTrans->GetSize(); // set to 12
    
    
    svec<int> counts;
    cerr << "Counting k-mers..." << endl;
    counts.resize(m_pTrans->GetBoundValue(), 0);
    
    for (j=0; j<(int)bases.size(); j++) {

        // iterating through contigs [j]

        const DNAVector & b = bases[j];
        const NumVector & t = tags[j];
        k=0;
        
        //cerr << "Read inchworm contig: " << b.AsString() << endl;
        
        if (j % 1000 == 0 || j == (int)bases.size()-1)
            cerr << "\rKmerAlignCore- Contigs: " << j << "   ";
        
        while (k <= (int)b.size()-size) {
            //cout << "k=" << k << endl;
            
            //cerr << "k=" << k << ", size=" << size << ", min=" << min << endl;
            
            // 

            if (!IsRepeat(t, k, size, min)) {
                KmerAlignCoreRecordStoreTable & t = m_table; // not used?
                int n = m_pTrans->BasesToNumber(b, k); // (bases, position_k): uniquely identifies a kmer
                if (n >= 0)
                    counts[n]++; // counting the kmer
            }
            k++;
        }
    }
    //cerr << "done, re-sizing arrays..." << endl;
    cerr << endl;
    
    // ensure can hold matching positions for each of the kmer positions.
    for (j=0; j<counts.isize(); j++) {
        KmerAlignCoreRecordStoreTable & t = m_table;
        KmerAlignCoreRecordStore & s = t[j];
        s.Resize(counts[j]);
        counts[j] = 0;    
    }
    
    cerr << "done, assigning k-mers..." << endl;
    
    for (j=0; j<(int)bases.size(); j++) {
        const DNAVector & b = bases[j];
        const NumVector & t = tags[j];
        k=0;
        if (j % 1000 == 0 || j == (int)bases.size()-1)
            cerr << "\rKmerAlignCore- Contigs: " << j << "   ";
        
        while (k <= (int)b.size()-size) {
            // walking position k of contig j
            if (!IsRepeat(t, k, size, min)) {
                KmerAlignCoreRecordStoreTable & t = m_table;
                int n = m_pTrans->BasesToNumber(b, k);
                if (n >= 0) {
                    KmerAlignCoreRecordStore & s = t[n];
                    s.Add(j, k, counts[n]); // for kmer (n), store contig (j) position (k) 'hit'
                    counts[n]++; // if > 0, then have more than the self-match.
                }
            }      
            k++;
        }
    }
    
    cerr << endl << "done!" << endl;
    
}

void KmerAlignCore::AddData(const DNAVector & b, int contig, int offset, bool bSort)
{
    int i, j, k;
    
    int size = m_pTrans->GetSize();
    k = 0;
    
    while (k < (int)b.size()-size) {
        KmerAlignCoreRecordStoreTable & t = m_table;
        int n = m_pTrans->BasesToNumber(b, k);
        if (n >= 0) {
            KmerAlignCoreRecordStore & s = t[n];
            s.Add(contig, k+offset);
        }
        k++;      
    }
    if (bSort)
        SortAll();
}

void KmerAlignCore::SortAll()
{
    int i, j;
    KmerAlignCoreRecordStoreTable & t = m_table;
    for (j=0; j<t.GetSize(); j++) {
        KmerAlignCoreRecordStore & s = t[j];
        s.Sort();    
    }
}

const svec<KmerAlignCoreRecord> & KmerAlignCore::GetMatchesDirectly(const DNAVector & b, int start)
{
    int size = m_pTrans->GetSize();
    int i, j;
    
    int n = m_pTrans->BasesToNumber(b, start);
    
    if (n < 0) {
        static svec<KmerAlignCoreRecord> dummy;
        return dummy;
    }
    
    KmerAlignCoreRecordStore & s = m_table[n];   
    return s.GetData();
}


bool KmerAlignCore::GetMatches(svec<KmerAlignCoreRecord> & matches, const DNAVector & b, int start) // start=0
{
    int size = m_pTrans->GetSize(); // set to 12
    int i, j, k, l;
    k = start;
    
    //cerr << "m_numTables=" << m_numTables << endl;
    
    if (start + m_numTables * size > (int)b.size()) {
        
        // kmer length > contig length?
        
        cerr << "Error: sequence length=" << b.size() << " and k-kmer end is " << start + m_numTables * size << endl; 
        return false;
    } 
    
    
    //cout << "Start position: " << start << endl;
    
    KmerAlignCoreRecordStoreTable hits;
    
    hits.SetSize(m_numTables);  // store first and second 12-mer hit
    
    for (i=0; i<m_numTables; i++) {
        KmerAlignCoreRecordStoreTable & t = m_table;
        KmerAlignCoreRecordStore & r = hits[i];
        
        if (m_lookAhead == 0) {

            int n = m_pTrans->BasesToNumber(b, k + i * size); // b=kmer_seq, k=start=0, size=12, i=[1,2]  (note only reads the first 12-mer of nucs)
            
            // 12-mer lookups.  Try first 12-mer, then second 12-mer of the 24-mer.??

            if (n >= 0) {
                KmerAlignCoreRecordStore & s = t[n];   
                if (s.GetNumRecords() > m_max12) {
                    matches.clear();
                    return false;
                }
                
                // Pre-size it!
                r.Resize(s.GetNumRecords());
                int count = 0;
                //cout << "size=" << s.GetNumRecords() << endl;

                // copy over hits from (s) to (r)
                for (j=0; j<s.GetNumRecords(); j++) {
                    //cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
                    r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size, count);
                    //cout << "Done" << endl;
                    count++;
                    //r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size);
                }
            }
        } else {
            
            int penalty = 0;
            for (l=0; l<m_lookAhead; l++) {
                if (k + i * size + l + size > (int)b.size())
                    break;
                //cout << "n=" << k + i * size + l << endl;
                int n = m_pTrans->BasesToNumber(b, k + i * size + l);
                
                if (n > 0) {
                    KmerAlignCoreRecordStore & s = t[n];   
                    
                    int count = r.GetSize();
                    if (l > 0)
                        penalty = 1;
                    
                    if (l > 0 && count > m_lookAheadMaxFreq) {
                        //cout << "***** Count: " << count << endl;
                        continue;
                    }
                    r.Resize(count + s.GetNumRecords());
                    for (j=0; j<s.GetNumRecords(); j++) {
                        //cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
                        r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size - l, count, penalty);
                        count++;
                        //r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size - l);
                    }
                }
                //cout << "Before sort: " << r.GetSize() << endl;
                r.UniqueSort();      
                //cout << "After sort:  " << r.GetSize() << endl;
                //for (l=0; l<r.GetSize(); l++) {
                //const KmerAlignCoreRecord & record = r.Get(l);
                //cout << "   c" << record.GetContig() << " @" << record.GetPosition() << endl;
                //}
            }
        }
    }
    
    matches.clear();
    svec<KmerAlignCoreRecord> tmp;
    
    if (hits.GetSize() == 0)
        return false;
    
    //cout << "Merging, hits size=" << hits.GetSize() << endl;
    tmp = hits[0].GetData();
    for (i=1; i<m_numTables; i++) {
        MergeSortFilter(matches, tmp, hits[i].GetData());
        tmp = matches;
        //cout << "Matches remaining: " << matches.isize() << endl;
    }
    
    if (m_numTables == 1)
        matches = tmp;
    //cout << "All done!" << endl;
    
    return (matches.isize() > 0);
}

void KmerAlignCore::MergeSortFilter(svec<KmerAlignCoreRecord> & result,
                                    const svec<KmerAlignCoreRecord> & one,
                                    const svec<KmerAlignCoreRecord> & two) {
    // Stupid for now...
    int i;
    
    result.clear();
    
    //cout << "one=" << one.isize() << "  two=" << two.isize() << endl;
    if (one.isize() == 0 || two.isize() == 0)
        return;
    
    svec<KmerAlignCoreRecord> tmp;
    result.resize(one.isize() + two.isize());
    tmp.resize(one.isize() + two.isize());
    
    int k = 0;
    
    /*
      for (i=0; i<one.isize(); i++)
      tmp[i] = one[i];
      for (i=0; i<two.isize(); i++)
      tmp[i+one.isize()] = two[i];
      ::Sort(tmp);
      
      
      int shift = m_pTrans->GetSize();
      for (i=0; i<tmp.isize()-1; i++) {
      if (tmp[i].GetContig() == tmp[i+1].GetContig()) {
      if (tmp[i].GetPosition() == tmp[i+1].GetPosition()) {
      result[k] = tmp[i];
      k++;
      }
      }
      } */
    
    
    int x = 0;
    int y = 0;
    for (x=0; x<one.isize(); x++) {
        while (y<two.isize() && two[y] < one[x]) {
            y++;
        }
        if (x >= one.isize() || y >= two.isize())
            break;
        if (one[x] == two[y]) {
            result[k] = one[x];
            k++;
        }
    }
  
    
    result.resize(k);
}
