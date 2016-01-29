/* -*- mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include "analysis/DNAVector.h"
#include "analysis/AACodons.h"
//#include "analysis/CodonTranslate.h"
#include "base/FileParser.h"
#include "util/mutil.h"

#define ONE_THIRD 1./3.

DNACodec::DNACodec() 
{
  int i;
  for (i=0; i<256; i++) {
    m_A[i] = m_C[i] = m_G[i] = m_T[i] = 0;
    m_rc[i] = 0;
  }

  Set('A', 1., 0., 0., 0., 'T');
  Set('C', 0., 1., 0., 0., 'G');
  Set('G', 0., 0., 1., 0., 'C');
  Set('T', 0., 0., 0., 1., 'A');



  Set('K', 0., 0., 0.5, 0.5, 'M');
  Set('M', 0.5, 0.5, 0., 0., 'K');

  Set('R', 0.5, 0., 0.5, 0., 'Y');
  Set('Y', 0., 0.5, 0., 0.5, 'R');

  Set('S', 0., 0.5, 0.5, 0., 'S');
  Set('W', 0.5, 0., 0., 0.5, 'W');

  //------------------------------------------
  //Set('K', 0., 0., 0.75, 0.75, 'M');
  //Set('M', 0.75, 0.75, 0., 0., 'K');

  //Set('R', 0.75, 0., 0.75, 0., 'Y');
  //Set('Y', 0., 0.75, 0., 0.75, 'R');

  //Set('S', 0., 0.75, 0.75, 0., 'S');
  //Set('W', 0.75, 0., 0., 0.75, 'W');
  //-----------------------------------------

  Set('B', 0., ONE_THIRD, ONE_THIRD, ONE_THIRD, 'V');
  Set('V', ONE_THIRD,ONE_THIRD , ONE_THIRD, 0., 'B');

  Set('H', ONE_THIRD, ONE_THIRD, 0., ONE_THIRD, 'D');
  Set('D', ONE_THIRD, 0., ONE_THIRD, ONE_THIRD, 'H');

  Set('-', 0., 0., 0., 0., '-');
  Set('N', 0.25, 0.25, 0.25, 0.25, 'N');
  Set('X', 0.25, 0.25, 0.25, 0.25, 'X');

}

char ResolveAmbiguous(const string & bases)
{
  const char * p = (const char*)bases.c_str();
  int n = strlen(p);
  if (n != 2)
    return GetAmbiguous(bases);

  bool bNuke0 = false;
  if (p[0] == 'A' || p[0] == 'C' || p[0] == 'G' || p[0] == 'T') {
    bNuke0 = true;
  }
  bool bNuke1 = false;
  if (p[1] == 'A' || p[1] == 'C' || p[1] == 'G' || p[1] == 'T') {
    bNuke1 = true;
  }

  if (bNuke0 && bNuke1)
    return GetAmbiguous(bases);

  if (bNuke0)
    return p[0];
  
  if (bNuke1)
    return p[1];
  
  
  return p[0];
}



char GetAmbiguous(const string & bases)
{
  bool a = false;
  bool c = false;
  bool g = false;
  bool t = false;

  const char * p = (const char*)bases.c_str();
  int i;
  int n = (int)strlen(p);
  for (i=0; i<n; i++) {
    if (p[i] == 'A')
      a = true;
    if (p[i] == 'C')
      c = true;
    if (p[i] == 'G')
      g = true;
    if (p[i] == 'T')
      t = true;
  }


  if (a && c && g && t)
    return 'N';

  if (a && g && t)
    return 'D';
  if (a && c && t)
    return 'H';
  if (a && c && g)
    return 'V';
  if (c && g && t)
    return 'B';


  if (a && t)
    return 'W';
  if (c && g)
    return 'S';
  if (c && t)
    return 'Y';
  if (a && g)
    return 'R';
  if (a && c)
    return 'M';
  if (g && t)
    return 'K';


  if (a)
    return 'A';
  if (c)
    return 'C';
  if (g)
    return 'G';
  if (t)
    return 'T';

  return '?';
}


void DNACodec::Set(int letter, double a, double c,double g, double t, int rc)
{
  m_rc[letter] = rc;
  m_A[letter] = a;
  m_C[letter] = c;
  m_G[letter] = g;
  m_T[letter] = t;
}


DNACodec theCodec;

char GetRC(const char c)
{
  return theCodec.GetRC(c);
}

double DNA_A(char l)
{
  return theCodec.A_Letter(l);
}

double DNA_C(char l)
{
  return theCodec.C_Letter(l);
}

double DNA_G(char l)
{
  return theCodec.G_Letter(l);
}

double DNA_T(char l)
{
  return theCodec.T_Letter(l);
}


double DNA_EqualAmb(char letter1, char letter2)
{
  if (letter1 == letter2) {
    if (letter1 != 'N')
      return 1.;
  }
  double a = theCodec.A_Letter(letter1) * theCodec.A_Letter(letter2);  
  double c = theCodec.C_Letter(letter1) * theCodec.C_Letter(letter2);  
  double g = theCodec.G_Letter(letter1) * theCodec.G_Letter(letter2);  
  double t = theCodec.T_Letter(letter1) * theCodec.T_Letter(letter2);  

  double sum = a + c + g + t;
  double plus = 0.;

  return (sum + plus);
}


double DNA_Equal(char letter1, char letter2)
{
  double a = theCodec.A_Letter(letter1) * theCodec.A_Letter(letter2);  
  double c = theCodec.C_Letter(letter1) * theCodec.C_Letter(letter2);  
  double g = theCodec.G_Letter(letter1) * theCodec.G_Letter(letter2);  
  double t = theCodec.T_Letter(letter1) * theCodec.T_Letter(letter2);  

  double sum = a + c + g + t;
  double plus = 0.;
  //if (sum < 1. && letter1 == letter2)
  //plus = 0.25;

  return (sum + plus);
}


double DNA_EqualEmph(char letter1, char letter2)
{
  double v = DNA_Equal(letter1, letter2);
  if (v > 0.26)
    v += 0.2;
  if (v > 1.0)
    v = 1.;
  return v;
}


double DNA_Diff(char letter1, char letter2)
{
  return (1 - DNA_Equal(letter1,letter2));
}


//====================================================
void DNAVector::SetToSubOf(const DNAVector & v, int start, int len)
{
  m_data.resize(len);
  int i;
  for (i=start; i<start+len; i++)
    m_data[i-start] = v[i];

  if (v.QualSize() > 0) {
    m_qual.resize(len);
    for (i=start; i<start+len; i++)
      m_qual[i-start] = v.Qual(i);
    
  }

}


void DNAVector::ReverseComplement(string & seq) {
  int n = seq.size();
  int i = 0;
  int j = n-1;

  n = (n+1)/2;

  for (i=0; i<n; i++) {
    char one = seq[i];
    char two = seq[j];
    one = theCodec.GetRC(one);
    two = theCodec.GetRC(two);

    seq[i] = two;
    seq[j] = one;
    
    j--;
    
  }

}



void DNAVector::ReverseComplement() {
  int n = m_data.isize();
  int i = 0;
  int j = n-1;

  n = (n+1)/2;

  for (i=0; i<n; i++) {
    char one = m_data[i];
    char two = m_data[j];
    one = theCodec.GetRC(one);
    two = theCodec.GetRC(two);

    m_data[i] = two;
    m_data[j] = one;
    
    j--;
    
  }

  if (m_qual.isize() > 0) {
    n = m_data.isize();
    i = 0;
    j = n-1;

    n = (n+1)/2;

    for (i=0; i<n; i++) {
      char one = m_qual[i];
      char two = m_qual[j];           
      m_qual[i] = two;
      m_qual[j] = one;
      
      j--;
    
    }
  }


}

void DNAVector::SetFromBases(const string & s)
{
  int n = strlen(s.c_str());
  m_data.resize(n);

  const char * p = s.c_str();
  int i;
  for (i=0; i<n; i++)
    m_data[i] = p[i];

}

void DNAVector::Proteinize()
{
  char * p = new char[m_data.isize()+1];
  int i;
  for (i=0; i<m_data.isize(); i++)
    p[i] = m_data[i];
  p[m_data.isize()] = 0;
  SetFromProteins(p);
  delete [] p;
}


void DNAVector::ToUpper()
{
  int i;
  for (i=0; i<m_data.isize(); i++)
    m_data[i] = (char)toupper(m_data[i]);
}

void DNAVector::SetFromProteins(const string & s)
{
  int n = strlen(s.c_str());
  m_data.resize(n*3);

  char tmp[64];

  const char * p = s.c_str();

  int i;
  for (i=0; i<n; i++) {
		AACodons::AminoAcidToBases(tmp, p[i]);
    m_data[3*i] = tmp[0];
    m_data[3*i+1] = tmp[1];
    m_data[3*i+2] = tmp[2];
  }

}

void DNAVector::Write(FILE * p) const
{
  int i;
  for (i=0; i<isize(); i++) {
    if (i > 0 && i % 80 == 0)
      fprintf(p, "\n");
    fprintf(p, "%c", m_data[i]);
    
  }
  fprintf(p, "\n");
}


void DNAVector::Write(ostream &s) const
{
  int i;
  for (i=0; i<isize(); i++) {
    if (i > 0 && i % 80 == 0)
      s << endl;
    s << m_data[i];
  }
  s << endl;
}


void DNAVector::WriteOneLine(ostream &s) const
{
  int i;
  for (i=0; i<isize(); i++) {
    s << m_data[i];
  }
  s << endl;
}



void DNAVector::WriteQual(ostream &s) const
{
  int i;
  for (i=0; i<isize(); i++) {
    if (i > 0 && i % 80 == 0)
      s << endl;
    s << (int) m_qual[i]<<" ";
  }
  s << endl;
}
void DNAVector::WriteQual(FILE * p) const
{
  int i;
  for (i=0; i<isize(); i++) {
    if (i > 0 && i % 80 == 0)
      fprintf(p, "\n");;
    fprintf(p, "%d ", m_qual[i]);
  }
  fprintf(p, "\n");;
 
}



bool DNAVector::Append(const DNAVector & d, int min, int max, double ident)
{
  int i, j;
  for (i=max; i>=min; i--) {
    int mis = (int)((1. - ident) * (double)i + 0.5);
    int no = 0;
    for (j=0; j<i; j++) {
      int x = isize()-i+j;
      if (x < 0)
	return false;
      //cout << (*this)[x] << " " << d[j] << " " << x << " " << j << endl;
      if ((*this)[x] != d[j]) {
	no++;
	if (no > mis)
	  break;
      }
    }
    if (no > mis)
      continue;

    //cout << "Merge, lap=" << i << endl;
    // Merge and get out.
    for (j=0; j<i; j++) {
      int x = isize()-i+j;
      if ((*this)[x] != d[j]) {
	string bases;
	bases = (*this)[x];
	bases += d[j];
	(*this)[x] = GetAmbiguous(bases);
      }
    }
    int oldSize = isize();
    resize(isize()-i+d.isize());
    //cout << "other size=" << d.isize() << " my size=" << isize() << endl;
    for (j=i; j<d.isize(); j++) {
      int x = oldSize-i+j;
      //cout << "Assign " << j << " to " << x << endl;
      (*this)[x] = d[j];
    }
    return true;
  }
  return false;
}



//////////////////////////////////////////////////////

vecDNAVector::vecDNAVector() {
	default_name_index = 0;
    m_pParser = NULL;
    m_lastName = "";
    max_seqs_to_read = -1; // no limit
}

//references intentionally not copied
vecDNAVector::vecDNAVector(const vecDNAVector &other) : m_names(other.m_names), m_data(other.m_data), default_name_index(other.default_name_index), m_name2index(other.m_name2index) {}

vecDNAVector::~vecDNAVector() {
	for(svec<string>::iterator currName = m_names.begin(); currName != m_names.end(); currName++)
		invalidateReferences(*currName);
    if(m_pParser != NULL){
        delete m_pParser;
    }
}

void vecDNAVector::setMaxSeqsToRead(long max_seq_count) {
    this->max_seqs_to_read = max_seq_count; // set to -1 to disable
    return;
}


//references intentionally not copied
vecDNAVector &vecDNAVector::operator = (const vecDNAVector &other) {
	if(this != &other) {
		for(svec<string>::iterator currName = m_names.begin(); currName != m_names.end(); currName++)
			invalidateReferences(*currName);
		m_names = other.m_names;
		m_data = other.m_data;
		default_name_index = other.default_name_index;
		m_name2index = other.m_name2index;
	}
	return *this;
}

const DNAVector &vecDNAVector::operator [] (int i) const {
  return m_data[i];
}

DNAVector &vecDNAVector::operator [] (int i) {
  return m_data[i];
}

const DNAVector &vecDNAVector::operator () (const string &name) const {
  return ((vecDNAVector*)this)->operator() (name);
}

DNAVector &vecDNAVector::operator () (const string &name) {
	  int index = NameIndex(name);
	  if(index == -1) {
        cout << "invalid sequence name " << name << endl;
        exit(-1);
	  }
	  return m_data[index];
}

vecDNAVector::const_DNAVectorRef vecDNAVector::getDNAVectorRef(int i) const {
  return const_DNAVectorRef(this, m_names[i]);
}

vecDNAVector::DNAVectorRef vecDNAVector::getDNAVectorRef(int i) {
	return DNAVectorRef(this, m_names[i]);
}

vecDNAVector::const_DNAVectorRef vecDNAVector::getDNAVectorRef(const string &name) const {
	  int index = NameIndex(name);
	  if(index == -1) {
        cout << "invalid sequence name " << name << endl;
        exit(-1);
	  }
	  return const_DNAVectorRef(this, string(name));
}

vecDNAVector::DNAVectorRef vecDNAVector::getDNAVectorRef(const string &name) {
	  int index = NameIndex(name);
	  if(index == -1) {
        cout << "invalid sequence name " << name << endl;
        exit(-1);
	  }
	  return DNAVectorRef(this, string(name));
}

void vecDNAVector::resize(int n) {
	int prevSize = m_data.isize();
  //If the size is decreased, erase the map entries for
  //indices beyond the new size and invalidate their references.
	if(n < prevSize)
		for(int i = n; i < prevSize; i++) {
			m_name2index.erase(m_names[i]);
			invalidateReferences(m_names[i]);
		}

  m_data.resize(n);
  m_names.resize(n);

	//If the size is increased, make new names for all the
	//new elements and add them to the names array and the map.
	if(n > prevSize) {
		for(int i = prevSize; i < n; i++) {
			stringstream currIndex;
			currIndex << default_name_index++;
			string currName = ">s_" + currIndex.str();
			m_names[i] = currName;
			m_name2index[currName] = i;
		}
	}
}

void vecDNAVector::clear() {
	for(int i = 0; i < m_names.isize(); i++)
		invalidateReferences(m_names[i]);

  m_data.clear();
  m_names.clear();
  m_name2index.clear();
}

int vecDNAVector::NameIndex(const string & name) const {
	  map<string,int>::const_iterator mIter = m_name2index.find(name);
	  if(mIter == m_name2index.end()) {
		  string realName;
		  string::size_type foundIndex = name.find('>');
		  if(foundIndex == string::npos)
			  realName = ">" + name;
		  else if(foundIndex == 0)
			  realName = name.substr(1);

		  mIter = m_name2index.find(realName);
		  if ( mIter == m_name2index.end() )
		  {
			  return -1;
		  }
	  }

	  return mIter->second;
}

const string & vecDNAVector::Name(size_t i) const {
	return m_names[i];
}

const char * vecDNAVector::NameClean(size_t i) const {
  const char * p = m_names[i].c_str();
  return &p[1];
}

void vecDNAVector::SetName(size_t i, const string & s) {
	m_name2index.erase(m_names[i]);
	invalidateReferences(m_names[i]);
	m_name2index[s] = i;
	m_names[i] = s;
}

void vecDNAVector::push_back(const DNAVector & v) {
  m_data.push_back(v);
  char tmp[256];
  sprintf(tmp, ">s_%d", default_name_index++);
  m_names.push_back(tmp);

  m_name2index.insert(make_pair(tmp, m_names.size()-1));
}

void vecDNAVector::push_back(const DNAVector & v, const string & name) {
  m_data.push_back(v);
  m_names.push_back(name);
  m_name2index.insert(make_pair(name, m_names.size()-1));
}

void vecDNAVector::erase(size_t index) {
	  m_name2index.erase(m_names[index]);
	  invalidateReferences(m_names[index]);
	  m_names.erase(m_names.begin() + index);
	  m_data.erase(m_data.begin() + index);
	  for(int i = index; i < m_names.isize(); i++)
		  m_name2index[m_names[i]] = i;
}

bool vecDNAVector::erase(const string &name) {
	int index = NameIndex(name);
	if(index != -1) {
		erase(index);
		return true;
	}
	return false;
}

void vecDNAVector::fast_erase(size_t index) {
	  m_name2index.erase(m_names[index]);
	  invalidateReferences(m_names[index]);
	  swap(m_names[m_names.size()-1], m_names[index]);
	  swap(m_data[m_data.size()-1], m_data[index]);
	  m_names.pop_back();
	  m_data.pop_back();
	  if(index < m_names.size())
		  m_name2index[m_names[index]] = index;
}

bool vecDNAVector::fast_erase(const string &name) {
	  int index = NameIndex(name);
	  if(index != -1) {
		  fast_erase(index);
		  return true;
	  }
	  return false;
}

size_t vecDNAVector::size() const {
	return m_data.size();
}

long long vecDNAVector::totalBases() const {
  long long sum(0);
  for (int i=0; i<(int) m_data.size(); ++i)
    sum += m_data[i].size();

  return sum;
}


void vecDNAVector::DoProteins(bool b) 
{
  int i, j;
  for (i=0; i<m_data.isize(); i++) {
    DNAVector & d = m_data[i];
    int n = d.isize() / 3;
    if (n * 3 != d.isize())
      cout << "WARNING: size of sequence " << m_names[i] << " is not a multiple of 3 (protein conversion: " << d.isize() << " )" << endl;

    char * p = new char[n+1];

    for (j=0; j<n; j++) {
      p[j] = AACodons::GetCodon(&d[j*3]);
    }
    p[j] = 0;
    //cout << "Set " << p << endl;
    if (!b) 
      d.SetFromProteins(p);
    else
      d.SetFromBases(p);
    delete [] p;
  }

}


void vecDNAVector::ReadV(const string & file) 
{
  int ver = 1;
  size_t i, j;
  CMReadFileStream s;
  s.Open(file.c_str());
  
  s.Read(ver);
  int size = 0;

  s.Read(size);
  resize(size);

  for (i=0; i<this->size(); i++) {
    CMString n;
    s.Read(n);
    m_names[i] = (const char*)n;
    DNAVector & d = m_data[i];
    m_name2index[m_names[i]] = i;
    size_t len = 0;
    s.Read(len);
    d.resize(len);
    for (j=0; j<len; j++)
      s.Read(d[j]);
  }

  s.Close();
}

void vecDNAVector::WriteV(const string & file) const
{
  int ver = 1;
  size_t i, j;
  CMWriteFileStream s;
  s.Open(file.c_str());
  
  s.Write(ver);
  s.Write(size());

  for (i=0; i<size(); i++) {
    CMString n = m_names[i].c_str();
    s.Write(n);
    const DNAVector & d = m_data[i];
    size_t len = d.size();
    s.Write(len);
    for (j=0; j<len; j++)
      s.Write(d[j]);
  }

  s.Close();
}



void vecDNAVector::Write(const string & fileName, bool bSkipEmpty) const
{
  FILE * p = fopen(fileName.c_str(), "w");
  for (size_t i=0; i<size(); i++) {
    if (bSkipEmpty && m_data[i].size() == 0)
      continue;
    //cout << "size: " << m_data[i].isize() << " skip=" << bSkipEmpty << endl;
    fprintf(p, "%s\n", m_names[i].c_str());
    m_data[i].Write(p);
    
  }
  fclose(p);
}
void vecDNAVector::WriteQuals(const string & fileName) const
{
  FILE * p = fopen(fileName.c_str(), "w");
  for (size_t i=0; i<size(); i++) {
    fprintf(p, "%s\n", m_names[i].c_str());
    m_data[i].WriteQual(p);
  }
  fclose(p);
}



void vecDNAVector::ReadQuals(const string & fileName)
{
  if (m_data.isize() == 0) {
    cout << "vecDNAVector ERROR: you need to load the bases first!" << endl;
    return;
  }
    
  FlatFileParser parser;
  
  parser.Open(fileName);
  
  int i = 0;
  int k = 0;
  int j;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const char * p = parser.AsString(0).c_str();
    if (p[0] == '>') {

      string tmpName = parser.AsString(0);
      for (int x=1; x<parser.GetItemCount(); x++) {
	tmpName += "_";
	tmpName += parser.AsString(x);
      }

      if (tmpName != m_names[i]) {
	cout << "vecDNAVector ERROR: qual file is out of sync with fasta file!" << endl;	
	return;
      }
      k = 0;
      ++i;
      continue;
    }
    for (j=0; j<parser.GetItemCount(); j++) {
      m_data[i-1].SetQual(k, parser.AsInt(j));
      k++;
    }
  }
}



// Half-way efficient implementation...?

int vecDNAVector::Read(const string & fileName, bool bProteins, bool shortName, bool allUpper, int bigChunk)
{
    if((max_seqs_to_read < 0) || (m_pParser == NULL) || (m_pParser->IsEndOfFile())){
        // Note: if the fileName is different, this will not currently reset the stream
        if (m_pParser != NULL){
            delete m_pParser;
        }
        this->clear();
        m_pParser = new FlatFileParser;
        m_pParser->Open(fileName);
        m_lastName = "";
    }

	int k = 0;
	int i;

	clear();

	// reserve some space
	int chunk = 20000;

	m_data.resize(chunk);

	DNAVector * pVec = NULL;
	DNAVector tmpVec;

	int j = 0;

    int num_seqs_read = 0;
    

	while (m_pParser->ParseLine()) {

        //if (num_seqs_read % 100 == 0)
        //    cerr << "\r[" << num_seqs_read << "] seqs read.  ";
        

		if (m_pParser->GetItemCount() == 0)
			continue;
		const char * p = m_pParser->AsString(0).c_str();
		if (p[0] == '>') {
			m_lastName = m_pParser->AsString(0);
			if ( !shortName )  {
				for (int x=1; x<m_pParser->GetItemCount(); x++) {
					m_lastName += "_";
					m_lastName += m_pParser->AsString(x);
				}
			}

            if ((max_seqs_to_read > 0) && (num_seqs_read >= max_seqs_to_read)) {
                break;
            }

			if (pVec != NULL) {
				pVec->SetToSubOf(tmpVec, 0, j);
				if (allUpper)
					pVec->ToUpper();
				j = 0;
				if (bProteins)
					pVec->Proteinize();
			}


			//cerr << "In: " << tmpName << endl;

			m_names.push_back(m_lastName);
            m_lastName = "";


			if (k >= m_data.isize())
				m_data.resize(k + chunk);

			pVec = &m_data[k];
            num_seqs_read++;
			k++;
		} else {
            if(m_lastName.length() != 0){
                //cerr << "Last name length is non-zero (" << m_lastName << ", " << p << "), must be continuing a stream" << endl;
                m_names.push_back(m_lastName);
                num_seqs_read++;
                k++;
            }
            int n = strlen(p);

            //cerr << "Parsing: " << p << endl;
            for (i=0; i<n; i++) {
                if (j >= tmpVec.isize())
                    tmpVec.resize(j + bigChunk);
                tmpVec[j] = p[i];
                j++;
            }
        }

	}

    //cerr << endl;


	if (pVec != NULL) {
		pVec->SetToSubOf(tmpVec, 0, j);
		if (allUpper)
			pVec->ToUpper();
		if (bProteins)
			pVec->Proteinize();
	}

	m_data.resize(k);
	//cerr << "Read in " << m_data.isize() << " sequences" << endl;
	//cerr << "Found " << m_names.isize() << " sequence headers" << endl;
	setupMap();
    if(m_pParser->IsEndOfFile()){
        return(-1 * num_seqs_read);
    } else {
        return(num_seqs_read);
    }
}

void vecDNAVector::Read(const string & fileName, svec<string> & names)
{
  Read(fileName);
  int i;
  names.resize(m_names.isize());
  for (i=0; i<m_names.isize(); i++) {
    const char * p = m_names[i].c_str();
    //if (strlen(p) > 0)
    names[i] = &p[1];
  }
}

void vecDNAVector::ReadPartial(const string & fileName, unsigned int firstToRead, unsigned int numToRead, bool bProteins, bool shortName, bool allUpper)
{
	if(numToRead == 0)
		return;

	FlatFileParser parser;
	parser.Open(fileName);

	clear();

	// reserve some space
	const int chunk = 20000;

	// 200 Megs?
	const int bigChunk = 200000000;

	m_data.resize(chunk);

	DNAVector tmpVec;

	unsigned int currSequenceLength = 0;
	int currIndex = -1;
	unsigned int numRead = 0;

	string currName;
	//When this terminates currName will hold the name of the first sequence to read
	while(currIndex < (int)firstToRead && parser.ParseLine()) {
		if (parser.AsString(0)[0] == '>') {
			currIndex++;
			currName = parser.AsString(0);
			if ( !shortName )  {
				for (int x=1; x<parser.GetItemCount(); x++) {
					currName += "_";
					currName += parser.AsString(x);
				}
			}
		}
	}

	while (numRead < numToRead && parser.ParseLine()) {
		if (parser.GetItemCount() == 0)
			continue;
		const char * p = parser.AsString(0).c_str();
		if (p[0] == '>') {
			m_names.push_back(currName);
			m_data[numRead].SetToSubOf(tmpVec, 0, currSequenceLength);
			currSequenceLength = 0;

			if (allUpper)
				m_data[numRead].ToUpper();
			if (bProteins)
				m_data[numRead].Proteinize();
			numRead++;

			currName = parser.AsString(0);
			if ( !shortName )  {
				for (int x=1; x<parser.GetItemCount(); x++) {
					currName += "_";
					currName += parser.AsString(x);
				}
			}


			if (numRead >= m_data.size())
				m_data.resize(numRead + chunk);
		}
		else {
			int n = strlen(p);

			for (int i=0; i<n; i++) {
				if (currSequenceLength >= (unsigned int)tmpVec.size())
					tmpVec.resize(currSequenceLength + bigChunk);
				tmpVec[currSequenceLength] = p[i];
				currSequenceLength++;
			}
		}
	}

	if (numRead < numToRead) {
		m_names.push_back(currName);
		m_data[numRead].SetToSubOf(tmpVec, 0, currSequenceLength);
		currSequenceLength = 0;

		if (allUpper)
			m_data[numRead].ToUpper();
		if (bProteins)
			m_data[numRead].Proteinize();

		numRead++;
	}

	m_data.resize(numRead);
	setupMap();
}

void vecDNAVector::ReverseComplement() {
  for (int i=0; i < m_data.isize(); i++)
    m_data[i].ReverseComplement();
}

void vecDNAVector::UniqueSort() {
	map<DNAVector, string> tempNameMap;
	for(size_t i = 0; i < m_data.size(); i++)
		tempNameMap[m_data[i]] = m_names[i];

  m_name2index.clear();
  ::UniqueSort(m_data);
  m_names.resize(m_data.size());
  for(size_t i = 0; i < m_data.size(); i++) {
  	m_name2index[tempNameMap[m_data[i]]] = i;
  	tempNameMap.erase(m_data[i]);
  }
  for(map<DNAVector, string>::iterator currRemoved = tempNameMap.begin(); currRemoved != tempNameMap.end(); currRemoved++)
  	invalidateReferences(currRemoved->second);
}

void vecDNAVector::addReference(ReferenceListener *newReference) const {
	  references[newReference->name()].push_back(newReference);
}

void vecDNAVector::removeReference(ReferenceListener *toRemove) const {
	map<string, vector<ReferenceListener *> >::iterator containingVector = references.find(toRemove->name());
	if(containingVector != references.end()) {
		for(vector<ReferenceListener *>::iterator currReference = containingVector->second.begin(); currReference != containingVector->second.end(); currReference++) {
			if(*currReference == toRemove) {
				containingVector->second.erase(currReference);
				if(containingVector->second.empty())
					references.erase(toRemove->name());
				break;
			}
		}
	}
}

void vecDNAVector::invalidateReferences(string toInvalidate) {
	  string realName = toInvalidate;
	map<string, vector<ReferenceListener *> >::iterator containingVector = references.find(realName);
	if(containingVector == references.end()) {
		  string::size_type foundIndex = toInvalidate.find('>');
		  if(foundIndex == string::npos)
			  realName = ">" + toInvalidate;
		  else if(foundIndex == 0)
			  realName = toInvalidate.substr(1);

		  containingVector = references.find(realName);
	}

	if(containingVector != references.end())
		for(vector<ReferenceListener *>::iterator currReference = containingVector->second.begin(); currReference != containingVector->second.end(); currReference++)
			(*currReference)->invalidate();

	references.erase(realName);
}

void vecDNAVector::setupMap() {
	if ( m_data.empty() )
		return;

	for ( int i=0; i<(int) m_names.size(); ++i )
		m_name2index.insert(make_pair(m_names[i],i));
}

vecDNAVector::DNAVectorRef::DNAVectorRef(vecDNAVector * owner, string DNAName) : owner(owner), DNAName(DNAName) {
	owner->addReference(this);
}

vecDNAVector::DNAVectorRef::DNAVectorRef(const DNAVectorRef &toCopy) : ReferenceListener() {
	  owner = toCopy.owner;
	  DNAName = toCopy.DNAName;
	  if(owner != NULL)
		  owner->addReference(this);
}

vecDNAVector::DNAVectorRef::DNAVectorRef() {
	  owner = NULL;
	  DNAName = "";
}

vecDNAVector::DNAVectorRef::~DNAVectorRef() {
	if(owner != NULL)
		owner->removeReference(this);
}

vecDNAVector::DNAVectorRef &vecDNAVector::DNAVectorRef::operator = (const DNAVectorRef &other) {
	if(this != &other) {
		if(owner != NULL)
			owner->removeReference(this);

		owner = other.owner;
		DNAName = other.DNAName;
		if(owner != NULL)
			owner->addReference(this);
	}
	return *this;
}

bool vecDNAVector::DNAVectorRef::operator == (const DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::DNAVectorRef::operator == (const const_DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::DNAVectorRef::operator != (const DNAVectorRef & other) const {
	return !(*this == other);
}

bool vecDNAVector::DNAVectorRef::operator != (const const_DNAVectorRef & other) const {
	return !(*this == other);
}

DNAVector &vecDNAVector::DNAVectorRef::operator *() {
	  if(owner == NULL) {
		  cout << "Tried to dereference an invalid DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return (*owner)(DNAName);
}

DNAVector *vecDNAVector::DNAVectorRef::operator ->() {
	  if(owner == NULL) {
		  cout << "Tried to dereference an invalid DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return &(*owner)(DNAName);
}

string vecDNAVector::DNAVectorRef::name() const {
	return DNAName;
}

bool vecDNAVector::DNAVectorRef::isInvalid() const {
	return owner == NULL;
}

void vecDNAVector::DNAVectorRef::invalidate() {
	owner = NULL;
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const vecDNAVector * owner, string DNAName) : owner(owner), DNAName(DNAName) {
	owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const const_DNAVectorRef &toCopy) : ReferenceListener() {
	owner = toCopy.owner;
	DNAName = toCopy.DNAName;
	if(owner != NULL)
		owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const DNAVectorRef &toCopy) : ReferenceListener() {
	owner = toCopy.owner;
	DNAName = toCopy.DNAName;
	if(owner != NULL)
		owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef() {
	  owner = NULL;
	  DNAName = "";
}

vecDNAVector::const_DNAVectorRef::~const_DNAVectorRef() {
	if(owner != NULL)
		owner->removeReference(this);
}

vecDNAVector::const_DNAVectorRef &vecDNAVector::const_DNAVectorRef::operator = (const const_DNAVectorRef &other) {
	if(this != &other) {
		if(owner != NULL)
			owner->removeReference(this);

		owner = other.owner;
		DNAName = other.DNAName;
		if(owner != NULL)
			owner->addReference(this);
	}
	return *this;
}

vecDNAVector::const_DNAVectorRef &vecDNAVector::const_DNAVectorRef::operator = (const DNAVectorRef &other) {
	if(owner != NULL)
		owner->removeReference(this);

	owner = other.owner;
	DNAName = other.DNAName;
	if(owner != NULL)
		owner->addReference(this);
	return *this;
}

bool vecDNAVector::const_DNAVectorRef::operator == (const DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::const_DNAVectorRef::operator == (const const_DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::const_DNAVectorRef::operator != (const DNAVectorRef & other) const {
	return !(*this == other);
}

bool vecDNAVector::const_DNAVectorRef::operator != (const const_DNAVectorRef & other) const {
	return !(*this == other);
}

const DNAVector &vecDNAVector::const_DNAVectorRef::operator *() const {
	  if(owner == NULL) {
		  cout << "Tried to dereference an invalid const_DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return (*owner)(DNAName);
}

const DNAVector *vecDNAVector::const_DNAVectorRef::operator ->() const {
	  if(owner == NULL) {
		  cout << "Tried to dereference an invalid const_DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }

	  return &(*owner)(DNAName);
}

string vecDNAVector::const_DNAVectorRef::name() const {
	return DNAName;
}

bool vecDNAVector::const_DNAVectorRef::isInvalid() const {
	return owner == NULL;
}

void vecDNAVector::const_DNAVectorRef::invalidate() {
	owner = NULL;
}

int CountNs(const DNAVector & d)
{
  int n = 0;
  for (int i=0; i<d.isize(); i++) {
    if (d[i] == 'N' || d[i] == 'X')
      n++;
  }
  return n;
}

bool IsHomopolymer(const DNAVector &d, double threshold)
{
  map<char,int> basecounts;
  for (int i=0; i<d.isize(); ++i)
  {
    basecounts[d[i]]++;
  }

  map<char,int>::iterator mIter = basecounts.begin();
  unsigned int max(basecounts.count(mIter->first));
  for (; mIter != basecounts.end(); ++mIter)
    if (basecounts.count(mIter->first) > max)
      max = basecounts.count(mIter->first);
  
  return ((double) max/d.isize() >= threshold);
}




//=================================================================

void vecDNAVectorStream::ReadStream(const string & fileName)
{
  if (m_pParser != NULL)
    delete m_pParser;

  m_pParser = new FlatFileParser;
  
  m_pParser->Open(fileName);

  while (m_pParser->ParseLine()) // skip first header
  {
	  if (m_pParser->GetItemCount() > 0)
	  {
		  const char * p = m_pParser->AsString(0).c_str();
		      if (p[0] == '>') {
		        break;
		      }
	  }
  }
}


const DNAVector & vecDNAVectorStream::Next()
{
    m_seq.resize(0);
    //cerr << "m_seq.size: " << m_seq.isize() << endl;
    while (m_pParser->ParseLine()) {
        if (m_pParser->GetItemCount() == 0)
            continue;
        string tmpName = m_pParser->AsString(0);
        if (tmpName.at(0) == '>') {
            for (int x=1; x<m_pParser->GetItemCount(); x++) {
                tmpName += "_";
                tmpName += m_pParser->AsString(x);
            }
            m_lastName = tmpName;
            break;
        }
        DNAVector tmp;
        tmp.SetFromBases(m_pParser->Line());
        m_seq += tmp;
    }

    // force to be in uppercase:
    m_seq.ToUpper();

    return m_seq;
}

const string & vecDNAVectorStream::getLastName()
{
    return m_lastName;
}


vecDNAVectorStream::~vecDNAVectorStream() 
{
  if (m_pParser != NULL)
    delete m_pParser;
}




//=================================================================
//ML: fast reader for fasta files, reads sequences and name into std:string

void DNAStringStreamFast::ReadStream(const string & fileName)
{
  m_ifs.open(fileName.c_str(), ifstream::in);
  
  getline(m_ifs,m_buf);

  // seek to first sequence header
  while(m_ifs.good() && m_buf[0]!='>')
  {
     getline(m_ifs,m_buf);
  }
}

void DNAStringStreamFast::Next(string & seq)
{
  getline(m_ifs,m_buf);
  while(m_ifs.good() && m_buf[0]!='>')
  {
     seq += m_buf;
     getline(m_ifs,m_buf);
  }
}


bool DNAStringStreamFast::NextToVector(vector<string> & seqv, vector<string> & namev)
{
  if(!m_ifs.good()) return false;
  // current buffer content is the name
  namev.push_back(m_buf);
  // get first line of sequence
  getline(m_ifs,m_buf);
  if(!m_ifs.good())
  {
    namev.pop_back();  // ensures consistency of name and seq vector, this only happens with broken files
    return false;
  }
  seqv.push_back(m_buf);
  // get next lines until next header and append to current sequence 
  getline(m_ifs,m_buf);
  while(m_ifs.good() && m_buf[0]!='>')
  {
     (seqv.back()) += m_buf;
     getline(m_ifs,m_buf);
  }
  return true;
}


void DNAStringStreamFast::formatReadNameString(const string & name, string & outname)
{
  outname = name;
  // strip " " at the beginning
  while(outname.length() && outname[0]==' ') outname.erase(0,1);
  // convert " " to "_"
  size_t pos=0;
  while((pos=outname.find_first_of(' ',pos))!=string::npos) outname[pos]='_';
  // strip " " at the end
  while(outname.length() && outname[outname.length()-1]==' ') outname.erase(outname.length()-1,1);
}


