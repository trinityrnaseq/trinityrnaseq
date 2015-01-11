
#include "base/FileParser.h"
#include <stdlib.h>




StringParser::StringParser()
{
}


StringParser::~StringParser()
{
  m_file.Close();
}

void StringParser::SetLine(const string & line, const string & delim)
{
  CMTokenizer t;
  t.AddDelimiter(delim.c_str());

  CMPtrStringList tokens;
  t.Tokenize(tokens, line.c_str());

  m_items.clear();

  m_items.resize(tokens.length());
  for (int i=0; i<tokens.length(); i++)
    m_items[i] = (const char*)*tokens(i);
}

void StringParser::SetLine(const string & line)
{
  CMPtrStringList tokens;
  Tokenize(tokens, line.c_str());

  m_items.clear();

  m_items.resize(tokens.length());
  for (int i=0; i<tokens.length(); i++)
    m_items[i] = (const char*)*tokens(i);
}

int StringParser::GetItemCount()
{
  return (int) m_items.size();
}

bool StringParser::IsString(int index)
{
  //always true....
  return true;
}


bool StringParser::IsInt(int index)
{
  const char * p = m_items[index].c_str();
  int n = strlen(p);

  for (int i=0; i<n; i++) {
    if (p[i] < '0' || p[i] > '9')
      return false;

  }
  return true;
}


bool StringParser::IsFloat(int index)
{
  const char * p = m_items[index].c_str();
  int n = strlen(p);

  for (int i=0; i<n; i++) {
    if (p[i] < '0' || p[i] > '9') {
      if (p[i] != '.')
	return false;
    }

  }
  return true;
}


const string & StringParser::AsString(int i)
{
  return m_items[i];
}


int StringParser::AsInt(int i)
{
  return atol(m_items[i].c_str());
}


double StringParser::AsFloat(int i)
{
  return atof(m_items[i].c_str());

}



//=======================================================
FlatFileParser::FlatFileParser() : StringParser()
{
}

FlatFileParser::FlatFileParser(const string & fileName)
{
  m_file.Open(fileName.c_str());
}

FlatFileParser::~FlatFileParser()
{
  m_file.Close();
}
void FlatFileParser::Open(const string & fileName)
{
  m_file.Open(fileName.c_str());    
}

bool FlatFileParser::Exists(const string &fileName)
{
  ifstream istrm(fileName.c_str());
  bool isFile = istrm.is_open();
  istrm.close();

  return isFile;

}

bool FlatFileParser::ParseLine()
{
  CMString line;
  m_file.ReadLine(line);
  if (m_file.IsEnd())
    return false;

  SetLine((const char*)line);

  m_line = line;
  
  return true;
}

bool FlatFileParser::GetLine(string & line2)
{
  CMString line;
  m_file.ReadLine(line);
  if (m_file.IsEnd())
    return false;

  line2 = (const char*)line;

  return true;
}

bool FlatFileParser::IsEndOfFile()
{
  return m_file.IsEnd();
}

void FlatFileParser::LoadVector(string &filename,
				vector<string> &element)
{
  element.clear();
  string line("");
  this->Open(filename);
  while (this->ParseLine())
  {
    element.push_back(m_line);
  }
}

void FlatFileParser::LoadSet(string &filename,
			     set<string> &element)
{
  element.clear();
  this->Open(filename);
  while (this->ParseLine())
  {
    element.insert(m_line);
  }
}






