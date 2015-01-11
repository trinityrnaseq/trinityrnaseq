/* -*- mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */ 

#include "analysis/DNAVector.h"
#include "analysis/AACodons.h"


// initialize variables for AACodons static object
// see [http://stackoverflow.com/questions/272900/c-undefined-reference-to-static-class-member]
bool AACodons::initialized = false;
svec<int> AACodons::m_table;
svec<AAAmb> AACodons::m_bases;
svec<char> AACodons::m_aminoAcids;
svec<string> AACodons::m_codonBases;

void AACodons::init()
{
    if(initialized){
		// don't initialize more than once
		return;
	}
  m_table.resize(256, -1);
  m_bases.resize(25);

  Set("A", "GCT", 0); 
  Set("A", "GCC", 0); 
  Set("A", "GCA", 0); 
  Set("A", "GCG", 0); 
  
  Set("L", "TTA", 1); 
  Set("L", "TTG", 1); 
  Set("L", "CTT", 1); 
  Set("L", "CTC", 1); 
  Set("L", "CTA", 1); 
  Set("L", "CTG", 1);
  
  Set("R", "CGT", 2); 
  Set("R", "CGC", 2); 
  Set("R", "CGA", 2); 
  Set("R", "CGG", 2); 
  Set("R", "AGA", 2); 
  Set("R", "AGG", 2);
  
  Set("K", "AAA", 3); 
  Set("K", "AAG", 3); 
  
  Set("N", "AAT", 4); 
  Set("N", "AAC", 4); 
  
  Set("M", "ATG", 5); 
  
  Set("D", "GAT", 6); 
  Set("D", "GAC", 6); 
  
  Set("F", "TTT", 7); 
  Set("F", "TTC", 7); 
  
  Set("C", "TGT", 8); 
  Set("C", "TGC", 8); 
  
  Set("P", "CCT", 9); 
  Set("P", "CCC", 9); 
  Set("P", "CCA", 9); 
  Set("P", "CCG", 9); 
  
  Set("Q", "CAA", 10); 
  Set("Q", "CAG", 10); 
  
  Set("S", "TCT", 11); 
  Set("S", "TCC", 11); 
  Set("S", "TCA", 11); 
  Set("S", "TCG", 11); 
  Set("S", "AGT", 11); 
  Set("S", "AGC", 11); 
  
  Set("E", "GAA", 12); 
  Set("E", "GAG", 12); 
  
  Set("T", "ACT", 13); 
  Set("T", "ACC", 13); 
  Set("T", "ACA", 13); 
  Set("T", "ACG", 13); 
  
  Set("G", "GGT", 14); 
  Set("G", "GGC", 14); 
  Set("G", "GGA", 14); 
  Set("G", "GGG", 14); 
  
  Set("W", "TGG", 15); 
  
  Set("H", "CAT", 16); 
  Set("H", "CAC", 16); 
  
  Set("Y", "TAT", 17); 
  Set("Y", "TAC", 17); 
  
  
  Set("I", "ATT", 18); 
  Set("I", "ATC", 18); 
  Set("I", "ATA", 18); 
  
  Set("V", "GTT", 19); 
  Set("V", "GTC", 19); 
  Set("V", "GTA", 19); 
  Set("V", "GTG", 19); 
  
  Set("*", "TAG", 20); 
  Set("*", "TGA", 20); 
  Set("*", "TAA", 20); 


  //Special/ambiguous 

  //B (N or D)
  Set("B", "AAT", 21); 
  Set("B", "AAC", 21); 
  Set("B", "GAT", 21); 
  Set("B", "GAC", 21); 


  //Z (Q or E)
  Set("Z", "GAA", 22); 
  Set("Z", "GAG", 22); 
  Set("Z", "CAA", 22); 
  Set("Z", "CAG", 22); 


  //J (L or I)
  Set("J", "ATT", 23); 
  Set("J", "ATC", 23); 
  Set("J", "ATA", 23); 
  Set("J", "TTA", 23); 
  Set("J", "TTG", 23); 
  Set("J", "CTT", 23); 
  Set("J", "CTC", 23); 
  Set("J", "CTA", 23); 
  Set("J", "CTG", 23);
  
 

  //X anything
  Set("X", "AAA", 24); 
  Set("X", "CCC", 24); 
  Set("X", "GGG", 24); 
  Set("X", "TTT", 24); 
	initialized = true;
}

const string& AACodons::GetBases(char aa, int i) {
	if(!initialized){init();}
	int index = 0;
	//cout << "Checking AA: " << aa << endl;
	if (aa < 0) {
		cout << "ERROR(1): Illegal amino acid: " << aa << endl;
	} else {
		index = m_table[aa];
	}
	if (index >= m_bases.isize() || index < 0) {
		cout << "ERROR(2): Illegal amino acid: " << aa << endl;
		index = 0;
	}
	const AAAmb & a = m_bases[index];
	return a.Base(i);
}


char AACodons::GetCodon(const DNAVector & d, int pos)
{
	if(!initialized){init();}
  char tmp[64];
  tmp[0] = d[pos];
  tmp[1] = d[pos+1];
  tmp[2] = d[pos+2];
  tmp[3] = 0;
  int i;
  
  //cout << "Asking for " << tmp << endl;
  
  for (i=0; i<m_codonBases.isize(); i++) {
    if (m_codonBases[i] == tmp)
      return m_aminoAcids[i];
  }
  cout << "WEIRD ERROR!" << endl;
  return (char)0xFF;
}

char AACodons::GetCodon(const char * p) {
	if(!initialized){init();}
	char tmp[64];
	tmp[0] = p[0];
	tmp[1] = p[1];
	tmp[2] = p[2];
	tmp[3] = 0;
	int i;
	
	//cout << "Asking for " << tmp << endl;
	
	for (i=0; i<m_codonBases.isize(); i++) {
		if (m_codonBases[i] == tmp)
			return m_aminoAcids[i];
	}
	cout << "WEIRD ERROR!" << endl;
	return (char)0xFF;
}

int AACodons::GetCodonCount() {
	return m_codonBases.isize();
}

const string& AACodons::GetCodon(int i) {
	return m_codonBases[i];
}

void AACodons::Set(const string & aa, const string & codon, int index) {
	char i = (aa.c_str())[0];
	m_table[i] = index;
	m_bases[index].Add(codon);
  
	m_aminoAcids.push_back(i);
	m_codonBases.push_back(codon);
}



void AACodons::AminoAcidToBases(char * out, char aa)
{
  //cout << "Query " << aa << endl;
  int i;
  for (i=0; i<3; i++) {
    const string & b = GetBases(aa, i);
    out[i] = GetAmbiguous(b);
    if (out[i] == 'N')
      out[i] = 'X';
  }
  out[4] = 0;
}



char AACodons::BasesToAminoAcid(char * b)
{
  return GetCodon(b);
}
