/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */ 
#ifndef _AACODONS_H_
#define _AACODONS_H_

#include <vector>
#include "base/SVector.h"
#include "analysis/DNAVector.h"

// forward declare some classes
// http://stackoverflow.com/questions/2133250/does-not-name-a-type-error-in-c
class DNAVector;
class AAAmb;

class AACodons {
public:
    static void init();
    static const string & GetBases(char aa, int i);
    static char GetCodon(const DNAVector& d, int pos);
    static char GetCodon(const char * p);
    static int GetCodonCount();
    static const string & GetCodon(int i);
    static void AminoAcidToBases(char * out, char aa);
    static char BasesToAminoAcid(char * b);
private:
    static bool initialized;
    static svec<int> m_table;
    static svec<AAAmb> m_bases;
    static svec<char> m_aminoAcids;
    static svec<string> m_codonBases;
    static void Set(const string & aa, const string & codon, int index);
};

#endif //_AACODONS_H_
