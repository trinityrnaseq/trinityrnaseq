#ifndef FORCE_DEBUG
#define NDEBUG
#endif

// Utilities based on an open source library...

using namespace std;

#include "util/mutil.h"
#include <stdlib.h>
#include <iostream>

#define PAGE_SIZE 131072

class CMMemPage
{
public:
    CMMemPage(long size = PAGE_SIZE);
    ~CMMemPage();

	void * New(long len) {
	  if (m_ptr + len < m_data.length()) {
	    m_ptr += len;
		m_objectCount++;
		return (void*)&m_data(m_ptr - len);
	  }
	  return NULL;
	}

	bool Delete(void * p) {
	  long mem = (long)p;
	  if (mem >= m_head && mem < m_tail) {
	  	m_objectCount--;
		return true;
	  }
	  return false;
	}

    long GetObjectCount() {return m_objectCount;}

private:
    CMCharList m_data;
	long m_ptr;
	long m_objectCount;
	long m_head;
	long m_tail;
};




class CMMemoryManager
{
public:
    CMMemoryManager();
    ~CMMemoryManager();

	void * New(long size);
	void Delete(void *);

	void NewPage();

	void Switch(bool bOn = true);

private:
	bool m_bActive;
	TMPtrList<CMMemPage> m_buffers;
};


CMMemoryManager & GetMemoryManager();


CMMemoryManager theMemMgr;


CMMemoryManager & GetMemoryManager()
{
  return theMemMgr;
}



CMMemPage::CMMemPage(long size)
{
  m_data.reshape(size);

  for (int i=0; i<m_data.length(); i++)
	m_data(i) = (char)0xAA;

  m_ptr = 0;
  m_objectCount = 0;
  m_head = (long)&m_data(0);
  m_tail = m_head + m_data.length();
}

CMMemPage::~CMMemPage()
{
}


//======================================================
//======================================================
//======================================================

CMMemoryManager::CMMemoryManager()
{
  m_bActive = false;
}

CMMemoryManager::~CMMemoryManager()
{
}

void CMMemoryManager::Switch(bool bOn)
{
  m_bActive = bOn;
  if (m_buffers.length() == 0) {
	m_buffers.add(new CMMemPage);
    MLog("MEMMGR: Adding page - count: ", m_buffers.length());
  }
}

void CMMemoryManager::NewPage()
{
  m_buffers.add(new CMMemPage);
  MLog("MEMMGR: Adding page - count: ", m_buffers.length());
}

void * CMMemoryManager::New(long size)
{
  if (!m_bActive) 
    return malloc(size);
  
  CMMemPage * pLast = m_buffers(m_buffers.length()-1);
  
  void * p = pLast->New(size);
  
  if (p == NULL) {
    pLast = new CMMemPage;
    m_buffers.add(pLast);  
    MLog("MEMMGR: Adding page - count: ", m_buffers.length());
    p = pLast->New(size);
    if (p == NULL)
      return malloc(size);
  }
  
  return p;
}

void CMMemoryManager::Delete(void * p)
{
  int i;
  for (i=0; i<m_buffers.length(); i++) {
    if (m_buffers(i)->Delete(p)) {
      if (m_buffers(i)->GetObjectCount() == 0) {
	m_buffers.remove(i);
        MLog("MEMMGR: Removing page - count: ", m_buffers.length());
      }
      return;	
    }
  }
  
  free(p);
}





MDLLEXPORT void ThrowException()
{
  MLog("*** EXCEPTION *** (NO REASON)");
#ifndef MCL_NO_EXCEPTIONS
  throw;
#endif //MCL_NO_EXCEPTIONS
}

MDLLEXPORT void ThrowException(const CMString & reason)
{
  MLog("*** EXCEPTION *** ", reason);
#ifndef MCL_NO_EXCEPTIONS
  throw CMException(reason);
#endif //MCL_NO_EXCEPTIONS
}

MDLLEXPORT void ThrowException(const CMString & comment, const CMString & reason)
{
  CMString text = comment;
  text += " ";
  text += reason;
  MLog("*** EXCEPTION *** ", text);
#ifndef MCL_NO_EXCEPTIONS
  throw CMException(text);
#endif //MCL_NO_EXCEPTIONS
}


CMException::CMException(const CMString & error) : m_text(error)
{
}

CMException::~CMException()
{
}

const CMString & CMException::GetErrorText()
{
  return m_text;
}

void CMException::Print()
{
  mlog() << "EXCEPTION: " << (const char*)m_text << ML_ENDL;
}

CMAsciiReadFileStream::CMAsciiReadFileStream()
{
  m_pFile = NULL;
  m_bytesProcessed = 0;
  m_bIsEof = false;
}

CMAsciiReadFileStream::~CMAsciiReadFileStream()
{
  Close();
}

IMReadStream * CMAsciiReadFileStream::CloneAndOpen(const CMString &name)
{
  if (!IsOpen())
	return NULL;

  CMAsciiReadFileStream * pNewStream = new CMAsciiReadFileStream;
  
  CMString newPath;
  CMString temp;
  const char * pTemp = (const char *) m_fileName;
  for (int i=0; i<(long)strlen(m_fileName); i++) {
    temp += pTemp[i];
	if (pTemp[i] == '\\' || pTemp[i] == '/') {
	  newPath = temp;
	}
  }

  newPath += name;

  MCL_TRY 
  {
	pNewStream->Open(newPath);

  }

#ifndef MCL_NO_EXCEPTIONS
  catch (CMException & ) {
    delete pNewStream;
	return NULL;
  }
#endif //MCL_NO_EXCEPTIONS

  return pNewStream;
}

bool CMAsciiReadFileStream::Open(const CMString & name)
{
  m_fileName = name;
  if (IsOpen())
    Close();
  m_pFile = fopen(m_fileName, "r");

  if (m_pFile != NULL) {
    m_bIsEof = false;
    return true;
  } else {
    mlog() << "Could not open file: " << (const char*)m_fileName << ML_ENDL;
	ThrowException("Could not open file for read: ", m_fileName);
    return false;
  }
}


bool CMAsciiReadFileStream::Close()
{
  m_bytesProcessed = 0;
  if (IsOpen()) {
    fclose(m_pFile);
    m_pFile = NULL;
  }
  m_bIsEof = true;
  return true;
}

bool CMAsciiReadFileStream::IsOpen()
{
  return (m_pFile != NULL);
}

bool CMAsciiReadFileStream::IsEnd()
{
  return m_bIsEof;
}

bool CMAsciiReadFileStream::ReadSimpleType(void * pData, long lenInBytes)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;


  bool bSuccess = true;

  char szText[2048 * 10];

  if (fscanf(m_pFile, szText, sizeof(szText), m_pFile) == EOF) {
    bSuccess = false;
  } else {
    //long val = atol(szText);
	//Do Something HERE!!
    //*pData = atol(szText);
  }
  m_bytesProcessed += strlen(szText);

  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  return bSuccess;

}


bool CMAsciiReadFileStream::ReadBlob(void * pData, long lenInElements, long elSize)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;

  bool bSuccess = true;


  if (fread(pData, lenInElements * elSize, 1, m_pFile) != 1)
    bSuccess = false;

  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  m_bytesProcessed += lenInElements * elSize;
  return bSuccess;

}

bool CMAsciiReadFileStream::ReadStringLine(CMString & out)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;


  char text[128];
  string all;
  size_t n;
  bool exit = false;
  do {
    text[sizeof(text)-2] = 1;
    if (fgets(text, sizeof(text), m_pFile) == NULL) {
      m_bIsEof = true;
      break;  
    }
    n = strlen(text);
    
    exit = false;
    if (n >=1 && text[n-1] == '\n') {
      text[n-1] = 0;
      exit = true;
    }
    
    all += text;
  } while(!exit);
  
  out = all.c_str();


  /*
  int bufSize = 16384;
  char * szText = new char[bufSize];
  for (int i=0; i<bufSize; i++)
    szText[i] = 0;

  if (fgets(szText, bufSize-1, m_pFile) == NULL) {
    m_bIsEof = true;
    delete szText;
    return false;  
  }

  if (IsUTF8(szText)) {
    MLog("WARNING: Switching to UTF-8 Mode!!");
    printf("WARNING: Switching to UTF-8 Mode!!\n");
    MCLSetUTF8Encode(true);
  }
  
  //szText[strlen(szText)-1] = 0;
  
  int m = strlen(szText);
  if (m > 0 && szText[m-1] == 0x0a) {
    szText[m-1] = 0;
    m--;
  }
  if (m > 0 && szText[m-1] == 0x0d) {
    szText[m-1] = 0;
    m--;
  }
  string = szText;
  m_bytesProcessed += m;

  delete [] szText;
  */
  return true;
}

bool CMAsciiReadFileStream::ReadString(CMString & string)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;

  long len;

  Read(len);

  char szSmallBuffer[512];

  char * pData = szSmallBuffer;

  //Dynamic only if really needed (it's slower!)
  if (len > (long)sizeof(szSmallBuffer)) {
    pData = new char[len];
  }

  bool bSuccess = true;

  if (fscanf(m_pFile, pData, len, m_pFile) == EOF) {
    bSuccess = false;
    string = "";
  } else {
    string = pData;
  }

  if (len > (long)sizeof(szSmallBuffer)) {
    delete [] pData;
  }
  m_bytesProcessed += len;
  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  return bSuccess;
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
CMAsciiWriteFileStream::CMAsciiWriteFileStream()
{
  m_pFile = NULL;
  m_bytesProcessed = 0;

  m_bIsOpen = false;
  m_bIsEof = false;
}

CMAsciiWriteFileStream::~CMAsciiWriteFileStream()
{
  Close();
}

bool CMAsciiWriteFileStream::Open(const CMString & name)
{
  m_fileName = name;
  if (IsOpen())
    Close();
  m_pFile = fopen(m_fileName, "w");

  if (m_pFile != NULL) {
    m_bIsEof = false;
    return true;
  } else {
    ThrowException("Could not open file for write: ", m_fileName);
    return false;
  }
}

bool CMAsciiWriteFileStream::Close()
{
  if (IsOpen()) {
    fclose(m_pFile);
    m_pFile = NULL;
  }
  m_bytesProcessed = 0;

  m_bIsEof = true;
  return true;
}

bool CMAsciiWriteFileStream::IsOpen()
{
  return (m_pFile != NULL);
}

bool CMAsciiWriteFileStream::IsEnd()
{
  return m_bIsEof;
}

bool CMAsciiWriteFileStream::WriteSimpleType(const void * pData, long lenInBytes)
{
  if (m_pFile == NULL)
    return false;
  //This is a DESIGN FLAW - that's why it is the way it is!
  if (lenInBytes <= 4)
    fprintf(m_pFile, "%ld ", *((long*)pData));
  else
    fprintf(m_pFile, "%f ", *((double*)pData));
  m_bytesProcessed += lenInBytes;
  return true;

}

bool CMAsciiWriteFileStream::WriteBlob(const void * pData, long lenInElements, long elSize)
{
  if (m_pFile == NULL)
    return false;

  //No BLOBs
  return false;
  //fprintf(m_pFile, "%lf ", *((double*)pData));

  //if (fwrite(pData, elSize * lenInElements, 1, m_pFile) != 1)
    //return false;

  return true;

}

bool CMAsciiWriteFileStream::WriteStringLine(const CMString & string)
{
  if (m_pFile == NULL)
    return false;

  fprintf(m_pFile, "%s\n", (const char*)string);
  m_bytesProcessed += strlen(string);
  return true;
}

bool CMAsciiWriteFileStream::WriteString(const CMString & string)
{
  if (m_pFile == NULL)
    return false;

  fprintf(m_pFile, "%s", (const char*)string);
  m_bytesProcessed += strlen(string);
  return true;
}

bool CMAsciiWriteFileStream::WriteLine(const CMString & line)
{
  if (m_pFile == NULL)
    return false;

  fprintf(m_pFile, "%s\n", (const char*)line);
  m_bytesProcessed += strlen(line);

  return true;
}






#define MLOGFILENAME "mccllog.txt"

static IMLogSink * g_pLogSink = NULL;

static CMLog g_mlog;


MDLLEXPORT CMLog & mlog()
{
  return g_mlog;
}



                                                             
CMLog::CMLog()
{
    mVerboseLevel = 1;
}

CMLog::~CMLog()
{
}
  


CMLog & CMLog::operator<<(unsigned char c) 
{
	(*this) << (char)c;
	return *this;
}

CMLog & CMLog::operator<<(signed char c) 
{
	(*this) << (char)c; 
	return *this;
}

CMLog & CMLog::operator<<(const unsigned char *s)
{ 
	(*this) << (const MCL_TCHAR*)s; 
    return *this;
}
CMLog & CMLog::operator<<(const signed char *s)
{ 
	(*this) << (const MCL_TCHAR*)s; 
	return *this;
}
 
CMLog & CMLog::operator<<(char c)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%ld", (long)c);
  MLog(szLog, false);
  return *this;

}

CMLog & CMLog::operator<<(const CMString & s)
{
  MLog((const MCL_TCHAR*)s, false);
  return *this;
}


CMLog & CMLog::operator<<(const MCL_TCHAR *s)
{
  MLog(s, false);
  return *this;
}

CMLog & CMLog::operator<<(const void *p)
{
  MLog((const MCL_TCHAR*)p, false);
  return *this;
}

CMLog & CMLog::operator<<(int n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%d", n);
  MLog(szLog, false);
  return *this;
}

CMLog & CMLog::operator<<(double n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%f", n);
  MLog(szLog, false);
  return *this;
}

CMLog & CMLog::operator<<(float n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%f", n);
  MLog(szLog, false);
  return *this;
}

CMLog & CMLog::operator<<(unsigned int n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%d", n);
  MLog(szLog, false);
  return *this;
}

CMLog & CMLog::operator<<(long n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%ld", n);
  MLog(szLog, false);
  return *this;
}

CMLog & CMLog::operator<<(unsigned long n)
{
  MCL_TCHAR szLog[32];
  __mccl_sprintf(szLog, "%lu", n);
  MLog(szLog, false);
  return *this;
}


CMLog & CMLog::operator<<(const VerbositySetting &rVs) 
{
    mVerboseLevel = rVs.getLevel();
    return *this;
}

//=====================================================================
//=====================================================================
//=====================================================================


#define NO_LOG

class CLogClear
{
public:
  CLogClear() 
  {
#ifndef NO_LOG
    FILE * pLog = fopen(MLOGFILENAME, "w");
    if (pLog)
      fclose (pLog);
#endif 
    
  }
};

CLogClear logClear;

MDLLEXPORT void SetLogSink(IMLogSink * p)
{
  g_pLogSink = p;
}

MDLLEXPORT void MLog(const MCL_TCHAR * szText, bool bLineFeed)
{
#ifndef NO_LOG
  if (g_pLogSink != NULL) {
    g_pLogSink->OnLog(szText, bLineFeed);
  } else {
    FILE * pLog = fopen(MLOGFILENAME, "a");
    if (pLog == NULL)
      return;
    if (bLineFeed) {
      __mccl_fprintf(pLog, "%s\n", szText);
      //printf("%s\n", szText);
    } else {
      __mccl_fprintf(pLog, "%s", szText);
      //printf("%s", szText);
    }
    fclose (pLog);
  }
#else
  cout << szText << endl;
#endif
}

MDLLEXPORT void MLog(const MCL_TCHAR * szText, long val, bool bLineFeed)
{
#ifndef NO_LOG
  if (g_pLogSink != NULL) {
    char szLog[2048];
    __mccl_sprintf(szLog, "%s %ld", szText, val);
    g_pLogSink->OnLog(szLog, bLineFeed);
  } else {
    FILE * pLog = fopen(MLOGFILENAME, "a");
    if (pLog == NULL)
      return;
    if (bLineFeed) {
      __mccl_fprintf(pLog, "%s %ld\n", szText, val);
      //printf("%s %ld\n", szText, val);
    } else {
      __mccl_fprintf(pLog, "%s %ld", szText, val);
      //printf("%s %ld", szText, val);
    }
    fclose (pLog);
  }
#else
  cout << szText << " " << val << endl;
#endif

}

MDLLEXPORT void MLog(const MCL_TCHAR * szText, const MCL_TCHAR * szText2, bool bLineFeed)
{
#ifndef NO_LOG
  if (g_pLogSink != NULL) {
    char szLog[2048];
    __mccl_sprintf(szLog, "%s %s", szText, szText2);
    g_pLogSink->OnLog(szLog, bLineFeed);
  } else {
    FILE * pLog = fopen(MLOGFILENAME, "a");
    if (pLog == NULL)
      return;
    if (bLineFeed) {
      __mccl_fprintf(pLog, "%s %s\n", szText, szText2);
      //printf("%s %s\n", szText, szText2);
    } else {
      __mccl_fprintf(pLog, "%s %s", szText, szText2);
      //printf("%s %s", szText, szText2);
    }
    fclose (pLog);
  }
#else
  cout << szText << " " << szText2 << endl;
#endif
}


CMReadFileStream::CMReadFileStream()
{
  m_pFile = NULL;
  m_bytesProcessed = 0;

  m_bIsEof = false;
}

CMReadFileStream::~CMReadFileStream()
{
  Close();
}

bool CMReadFileStream::Open(const CMString & name)
{
  m_fileName = name;
  if (IsOpen())
    Close();
  m_pFile = fopen(m_fileName, "rb");

  if (m_pFile != NULL) {
    m_bIsEof = false;
    return true;
  } else {
    ThrowException("Could not open file for read:", m_fileName);
    return false;
  }
}


bool CMReadFileStream::Close()
{
  if (IsOpen()) {
    fclose(m_pFile);
    m_pFile = NULL;
  }
  m_bIsEof = true;
  m_bytesProcessed = 0;
  return true;
}

bool CMReadFileStream::IsOpen()
{
  return (m_pFile != NULL);
}

bool CMReadFileStream::IsEnd()
{
  return m_bIsEof;
}

bool CMReadFileStream::ReadSimpleType(void * pData, long lenInBytes)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;


  bool bSuccess = true;
  if (fread(pData, lenInBytes, 1, m_pFile) != 1) {
    bSuccess = false;
  }

  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  m_bytesProcessed += lenInBytes;

  return bSuccess;

}


bool CMReadFileStream::ReadBlob(void * pData, long lenInElements, long elSize)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;

  bool bSuccess = true;
  if (fread(pData, lenInElements * elSize, 1, m_pFile) != 1)
    bSuccess = false;

  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  m_bytesProcessed += lenInElements * elSize;

  return bSuccess;

}

bool CMReadFileStream::ReadStringLine(CMString & string)
{
  return ReadString(string);
}

bool CMReadFileStream::ReadString(CMString & string)
{
  if (m_pFile == NULL || m_bIsEof)
    return false;

  long len;

  Read(len);

  char szSmallBuffer[512];

  char * pData = szSmallBuffer;

  //Dynamic only if really needed (it's slower!)
  if (len > (long)sizeof(szSmallBuffer)) {
    pData = new char[len];
  }

  bool bSuccess = true;
  if (fread((void *)pData, len, 1, m_pFile) != 1) {
    bSuccess = false;
    string = "";
  } else {
    string = pData;
  }

  if (len > (long)sizeof(szSmallBuffer)) {
    delete [] pData;
  }

  m_bytesProcessed += len;

  if (feof(m_pFile) != 0)
    m_bIsEof = true;

  return bSuccess;
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
CMWriteFileStream::CMWriteFileStream()
{
  m_pFile = NULL;

  m_bytesProcessed = 0;
  m_bIsOpen = false;
  m_bIsEof = false;
}

CMWriteFileStream::~CMWriteFileStream()
{
  Close();
}

bool CMWriteFileStream::Open(const CMString & name)
{
  m_fileName = name;
  if (IsOpen())
    Close();
  m_pFile = fopen(m_fileName, "wb");

  if (m_pFile != NULL) {
    m_bIsEof = false;
    return true;
  } else {
    ThrowException("Could not open file for write:", m_fileName);
    return false;
  }
}

bool CMWriteFileStream::Close()
{
  if (IsOpen()) {
    fclose(m_pFile);
    m_pFile = NULL;
  }
  m_bytesProcessed = 0;

  m_bIsEof = true;
  return true;
}

bool CMWriteFileStream::IsOpen()
{
  return (m_pFile != NULL);
}

bool CMWriteFileStream::IsEnd()
{
  return m_bIsEof;
}

bool CMWriteFileStream::WriteSimpleType(const void * pData, long lenInBytes)
{
  if (m_pFile == NULL)
    return false;

  if (fwrite(pData, lenInBytes, 1, m_pFile) != 1)
    return false;

  m_bytesProcessed += lenInBytes;

  return true;

}

bool CMWriteFileStream::WriteBlob(const void * pData, long lenInElements, long elSize)
{
  if (m_pFile == NULL)
    return false;

  if (fwrite(pData, elSize * lenInElements, 1, m_pFile) != 1)
    return false;

  m_bytesProcessed += elSize * lenInElements;

  return true;

}

bool CMWriteFileStream::WriteString(const CMString & string)
{
  if (m_pFile == NULL)
    return false;

  long len = strlen(string) + 1;

  Write(len);


  if (fwrite((void *)((const char *)string), len, 1, m_pFile) != 1)
    return false;

  m_bytesProcessed += len;

  return true;
}







bool g_UseCasing = true;

MDLLEXPORT void SetMemoryManage(bool b)
{
  GetMemoryManager().Switch(b);
}

#ifndef WIN32
void __mccl_tolwr(MCL_TCHAR * szText)
{
  long len = strlen(szText); 
  for (int i=0; i<len; i++) {
    if (szText[i] <= 'Z' && szText[i] >= 'A')
      szText[i] += 32;
  }
}

void __mccl_toupr(MCL_TCHAR * szText)
{
  long len = strlen(szText); 
  for (int i=0; i<len; i++) {
    if (szText[i] <= 'z' && szText[i] >= 'a')
      szText[i] -= 32;
  }
}
#endif //WIN32


MDLLEXPORT void MCLSetCasing(bool b)
{
  g_UseCasing = b;
}


//===============================================================
CMString::CMString()
{
  m_pData = NULL;
}

CMString::CMString(const CMString& S)
{
  m_pData = NULL;
  *this = S;
}


CMString::CMString(const MCL_TCHAR * S)
{
  m_pData = NULL;
  *this = S;
}

CMString::~CMString()
{
  DeleteMemory(m_pData);
}

MCL_TCHAR * CMString::GetMemory(long charLen) const
{
  char * p = (char*)GetMemoryManager().New(charLen * sizeof(MCL_TCHAR));
  return p;
  //MCL_TCHAR * pNew;
  //pNew = new MCL_TCHAR[charLen];
  //return pNew;
}

void CMString::DeleteMemory(MCL_TCHAR * pData) const
{
  if (pData == NULL)
	return;

  GetMemoryManager().Delete((void*)pData);

  //delete [] m_pData;
}

bool CMString::IsEmpty() const
{
  if (m_pData == NULL)
    return true;
  else
    return false;
}

CMString::operator const MCL_TCHAR*() const
{
  if (m_pData == NULL) {
    ((CMString *) this)->m_pData = GetMemory(1);
    __mccl_strcpy(m_pData, "");
  }
  return m_pData;
}

// Assignment:
CMString& CMString::operator=(const MCL_TCHAR* s)
{
  if (m_pData == s)
    return * this;

  if (m_pData != NULL)
    DeleteMemory(m_pData);

  if (s != NULL && __mccl_strcmp(s, _TMCL("")) != 0) {
    long newLen = __mccl_strlen(s);
    m_pData = GetMemory(newLen + 1);
    __mccl_strcpy(m_pData, s);
  } else {
    m_pData = NULL;
  }
  return * this;
}


CMString& CMString::operator=(const CMString & s)
{
  if (s.IsEmpty()) {
    if (m_pData != NULL) {
      DeleteMemory(m_pData);
      m_pData = NULL;
    }
    return *this;
  }
  
  return operator=((const MCL_TCHAR*)s);
}

CMString& CMString::operator+=(const MCL_TCHAR * s)
{
  long oldLen = 0;
  
  if (m_pData != NULL) {
    oldLen = __mccl_strlen(m_pData);
    //DeleteMemory(m_pData);
  }

  long newLen = oldLen + __mccl_strlen(s) + 1;

  MCL_TCHAR * pTmpString = GetMemory(newLen);
  if (newLen > 0)
    pTmpString[0] = 0;

  if (m_pData != NULL) {
    __mccl_strcpy(pTmpString, m_pData);
    DeleteMemory(m_pData);
  }

  __mccl_strcat(pTmpString, s);

  m_pData = pTmpString;

  return *this;
}


CMString& CMString::operator+=(const CMString& s)
{

  return operator += ((const MCL_TCHAR *)s);
}


// Indexing operators:
MCL_TCHAR & CMString::operator[](long i)
{
  if (m_pData == NULL) {

    ThrowException(_TMCL("Invalid string access (NULL)")); 

    return m_pData[0];
  }

  if (i > (long)__mccl_strlen(m_pData)) {
    ThrowException(_TMCL("Invalid string access..."));

    return m_pData[0];
  
  }
  
  return m_pData[i];
}

MCL_TCHAR & CMString::operator()(long i)
{
  if (m_pData == NULL) {
    ThrowException(_TMCL("Invalid string access (NULL)"));

    return m_pData[0];
  }
  return m_pData[i];
}

CMString& CMString::operator+=(const MCL_TCHAR c)
{
  MCL_TCHAR szApp[2];
  szApp[1] = 0;
  szApp[0] = c;
  return operator+=(szApp);
}


MCL_TCHAR CMString::operator[](long i) const
{
  if (m_pData == NULL) {
    ThrowException(_TMCL("Invalid string access (NULL)"));
  
    return m_pData[0];
  }

  if (i > (long)__mccl_strlen(m_pData)) {
    ThrowException(_TMCL("Invalid string access..."));

    return m_pData[0];
  
  }
  
  return m_pData[i];
}

MCL_TCHAR CMString::operator()(long i) const
{
  if (m_pData == NULL) {
    ThrowException(_TMCL("Invalid string access (NULL)"));

    return m_pData[0];
  }
  return m_pData[i];
}

static bool g_bUseUTF8 = false;

MDLLEXPORT void MCLSetUTF8Encode(bool b)
{
  g_bUseUTF8 = true;
}

MDLLEXPORT bool MCLIsUTF8Encode()
{
  return g_bUseUTF8;
}



void CMString::toLower()
{


  if (!g_UseCasing)
	return;


#ifdef UNIX
  ThrowException();
#else

  if (m_pData != NULL) {
	long len = strlen(m_pData);

	if (!g_bUseUTF8) {
	  for (int i=0; i<len; i++) {
	    if (m_pData[i] >= 'A' && m_pData[i] <= 'Z') {
	      m_pData[i] += 0x20;
		}	else {
	      if (m_pData[i] == 'Ä')
		    m_pData[i] = 'ä';
	      if (m_pData[i] == 'Ü')
		    m_pData[i] = 'ü';
	      if (m_pData[i] == 'Ö')
		    m_pData[i] = 'ö';
		}		 
	    //__mccl_tolwr(m_pData);
	  } 
	} else {
	  	//UTF-8!!!!!

	  for (int i=0; i<len; i++) {
		if ((m_pData[i] & 0x80) != 0) {
		  i++;
		}

	    if (m_pData[i] >= 'A' && m_pData[i] <= 'Z') {
	      m_pData[i] += 0x20;
		}
	  }
	}

  }
#endif
}

void CMString::toUpper()
{
  if (!g_UseCasing)
	return;
#ifdef UNIX
  ThrowException();
#else
  if (m_pData != NULL) {
	  __mccl_toupr(m_pData);

  }
#endif
}

bool CMString::operator == (const CMString& s1) const
{
  return operator == ((const MCL_TCHAR *) s1);
}

bool CMString::operator == (const MCL_TCHAR * s1) const
{
  if (m_pData == NULL) {
    if (s1 == NULL || __mccl_strcmp(s1, _TMCL("")) == 0)
      return true;
    else
      return false;
  }
  if (__mccl_strcmp(m_pData, s1) == 0)
    return true;
  else
    return false;


}

bool CMString::operator != (const CMString& s1) const
{
  return operator != ((const MCL_TCHAR *) s1);
}

bool CMString::operator != (const MCL_TCHAR * s1) const
{
  if (*this == s1)
    return false;
  else
    return true;
}


bool CMString::operator > (const CMString& s1) const
{
  return operator > ((const MCL_TCHAR *) s1);
}

bool CMString::operator > (const MCL_TCHAR * s1) const
{
  if(__mccl_strcmp(m_pData,s1) > 0)
    return true;
  else
    return false;
}
  
bool CMString::operator < (const CMString& s1) const
{
  return operator < ((const MCL_TCHAR *) s1);
}

bool CMString::operator < (const MCL_TCHAR * s1) const
{
  if(__mccl_strcmp(m_pData,s1) < 0)
    return true;
  else
    return false;
}



bool CMString::operator >= (const CMString& s1) const
{
  return operator >= ((const MCL_TCHAR *) s1);
}

bool CMString::operator >= (const MCL_TCHAR * s1) const
{
  if(__mccl_strcmp(m_pData,s1) >= 0)
    return true;
  else
    return false;
}
  

bool CMString::operator <= (const CMString& s1) const
{
  return operator <= ((const MCL_TCHAR *) s1);
}

bool CMString::operator <= (const MCL_TCHAR * s1) const
{
  if(__mccl_strcmp(m_pData,s1) <= 0)
    return true;
  else
    return false;
}

long CMString::length() const
{
  if (m_pData == NULL)
	return 0;
  return __mccl_strlen(m_pData);
}

void CMString::removeLeadingChars(char c)
{
  if (m_pData == NULL)
	return;
  long len = length();

  int i;
  for (i=0; i<len; i++) {
    if (m_pData[i] != c)
	  break;
  } 

  if (i > 0) {
    memmove((void*)m_pData, (void*)&m_pData[i], len + 1 - i);
  }
}

void CMString::removeTrailingChars(char c)
{
  if (m_pData == NULL)
	return;
  long i = length() - 1;
  while (i > 0 && m_pData[i] == c) {
    m_pData[i] = 0;
	i--;
  }
}

void CMString::removeSpaces()
{
  removeLeadingChars();
  removeTrailingChars();
}


MDLLEXPORT bool Tokenize(CMPtrStringList & result, const CMString & source, char delimiter, long limit)
{
  const MCL_TCHAR * pszBuffer = source;

  MCL_TCHAR * szTemp = new MCL_TCHAR[strlen(source)+4];
  //Clear out any previous results in the tokenizer
  result.removeAll();
  __mccl_strcpy(szTemp, _TMCL(""));

  int k = 0;

  long len = (long)__mccl_strlen(pszBuffer);

  for (int i=0; i<=len; i++) {
	if (pszBuffer[i] == delimiter || pszBuffer[i] == '\t' || i == len) {
	  k = 0;
	  if (__mccl_strcmp(szTemp, _TMCL("")) != 0 && __mccl_strcmp(szTemp, _TMCL(" ")) != 0) {
	  	CMString * pNewString = new CMString;
		*pNewString = szTemp;
		result.add(pNewString);
        __mccl_strcpy(szTemp, _TMCL(""));
		if (result.length() >= limit)
		  break;
	  }
	} else {
	  szTemp[k] = pszBuffer[i];
	  szTemp[k + 1] = 0;
	  k++;
	}
  }
  delete [] szTemp;

  return true;
}

MDLLEXPORT bool GetNextToken(CMString & result, CMString & source, char delimiter)
{
  MCL_TCHAR * pszBuffer = (MCL_TCHAR *)(const MCL_TCHAR *)source;
  MCL_TCHAR szTemp[2048];

  __mccl_strcpy(szTemp, _TMCL(""));

  int k = 0;

  for (int i=0; i<(int)__mccl_strlen(pszBuffer); i++) {
	if (pszBuffer[i] == delimiter || i+1 == (int)__mccl_strlen(pszBuffer)) {
	  k = 0;
	  if (__mccl_strcmp(szTemp, _TMCL("")) != 0) {
	  	result = szTemp;
		memmove((void*)pszBuffer, (void*)&pszBuffer[i+1], __mccl_strlen(pszBuffer) + 1 - i);
		return true;
	  }
	} else {
	  szTemp[k] = pszBuffer[i];
	  szTemp[k + 1] = 0;
	  k++;
	}
  }

  return false;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
CMTokenizer::CMTokenizer()
{
}

CMTokenizer::~CMTokenizer()
{
}

  
void CMTokenizer::AddDelimiter(const CMString & delimiter, const CMString & replacement)
{
  m_map.add(new CMStringMap(delimiter, replacement));
}


bool CMTokenizer::IsDelimiter(const char * pBuffer, long & inc, long & mapIndex)
{
  inc = 0;
  long cmpLen = strlen(pBuffer);
  for (int i=0; i<m_map.length(); i++) {
    long len = strlen(m_map(i)->GetString());
	if (len > cmpLen)
	  continue;
	if (__mccl_strncmp(pBuffer, m_map(i)->GetString(), len) == 0) {
	  if (len > inc) {
		inc = len;
		mapIndex = i;
	  }
	}
  }

  return (inc > 0);
}

bool CMTokenizer::Tokenize(CMPtrStringList & result, const CMString & source)
{
  const MCL_TCHAR * pszBuffer = source;

  MCL_TCHAR szTemp[4096 * 8];

  //Clear out any previous results in the tokenizer
  result.removeAll();
  __mccl_strcpy(szTemp, _TMCL(""));

  int k = 0;

  long len = (long)__mccl_strlen(pszBuffer);

  for (int i=0; i<=len; i++) {
    long inc = 0;
    long idx = 0;
	if (IsDelimiter(&pszBuffer[i], inc, idx) || i == len) {
	  
	  k = 0;
	  if (__mccl_strcmp(szTemp, _TMCL("")) != 0) {
	  	CMString * pNewString = new CMString;
		*pNewString = szTemp;
		result.add(pNewString);
        __mccl_strcpy(szTemp, _TMCL(""));
	  }
	  if (i < len) {
	    const CMString & replacement = m_map(idx)->GetMap();
	    if (replacement != "") {
	      result.add(new CMString(replacement));
		}

	    i += inc - 1;
	  }

	} else {
	  szTemp[k] = pszBuffer[i];
	  szTemp[k + 1] = 0;
	  k++;
	}
  }

  return true;
}
const long utf8SigLength = 3;
const char utf8Signature[] = {
  (char)0xef,
  (char)0xbb,
  (char)0xbf,
  (char)0x00,
};

MDLLEXPORT const char * GetUTF8Sig()
{
  return utf8Signature;
}

MDLLEXPORT bool IsUTF8(const CMString & string)
{
  if (strlen(string) < 3)
	return false;

  const char * p = (const char*)string;
  if (p[0] != utf8Signature[0])
	return false;
  if (p[1] != utf8Signature[1])
	return false;
  if (p[2] != utf8Signature[2])
	return false;

  return true;
}

MDLLEXPORT bool RemoveUTF8Sig(CMString & string)
{
  if (!IsUTF8(string))
	return false;

  const char * p = (const char*)string;

  if (p[3] == 0) {
	string = "";
	return true;
  }

  CMString tmp = &p[3];
  string = tmp;
  return true;
}

MDLLEXPORT bool AddUTF8Sig(CMString & string)
{
  if (IsUTF8(string))
	return false;

  const char * p = (const char*)string;

  CMString tmp = GetUTF8Sig();
  tmp += string;
  string = tmp;
  return true;
}



//===========================================================

CMStringDictionary::CMStringDictionary()
{
  m_wordCount = 0;
}

CMStringDictionary::~CMStringDictionary()
{
}

long CMStringDictionary::GetDictID(const CMString & word)
{
  long index;
  if (!m_searcher.Search(index, word))
    return INVALID_STRING_DICT_ID;
  else
    return index;
}

CMStringDictionary & CMStringDictionary::operator = (const CMStringDictionary & dict)
{
  m_dict = dict.m_dict;
  m_wordCount = dict.m_wordCount;
  m_searcher = dict.m_searcher;

  m_searcher.SetDataPtr(&m_dict);

  return *this;
}

long CMStringDictionary::GetDictIDList(CMInt32List & result, const CMString & word)
{
  long index;
  if (!m_searcher.InternalSearch(index, word)) {
    result.reshape(0);
    return 0;
  }
  
  long k = 0;

  //if (result.length() == 0)
  result.reshape(1);
 
  result(0) = m_searcher.IndexToIndex(index);
  k++;
  //down

  long indexStore = index;
  index--;
  while (index >= 0 && GetWord(m_searcher.IndexToIndex(index)) == word) {
    if (result.length() <= k)
      result.reshape(k+1);
    result(k) = m_searcher.IndexToIndex(index);
    index--;
    k++;
  }

  //up
  index = indexStore + 1;
  while (index < m_wordCount && GetWord(m_searcher.IndexToIndex(index)) == word) {
    if (result.length() <= k)
      result.reshape(k+1);
    result(k) = m_searcher.IndexToIndex(index);
    index++;
    k++;
  }
  if (result.length() > k)
    result.reshape(k);
  return k;
}

const CMString & CMStringDictionary::GetWord(long id)
{
  return m_dict(id);
}

long CMStringDictionary::GetWordCount()
{
  return m_wordCount;
}

long CMStringDictionary::GetDictIDByIndex(long index)
{
  return index;
}

const CMString & CMStringDictionary::GetWordByIndex(long index)
{
  return GetWord(GetDictIDByIndex(index));
}

long CMStringDictionary::AddWordDontCheck(const CMString & word)
{

  if (m_wordCount >= m_dict.length()) {
    m_dict.reshape(m_wordCount + 512);
	m_searcher.SetDataPtr(&m_dict);
  }

  m_dict(m_wordCount) = word;

  m_searcher.Add(m_dict(m_wordCount), m_wordCount);
  m_wordCount++;

  return m_wordCount-1;
}

long CMStringDictionary::AddWord(const CMString & word)
{
  long id = GetDictID(word);
  if (id != INVALID_STRING_DICT_ID)
    return id;

  return AddWordDontCheck(word);

}


bool CMStringDictionary::Read(IMReadStream & stream)
{
  CMFileHeader header;
  header.Read(stream);
  
  stream.Read(m_wordCount);

  m_dict.reshape(m_wordCount);

  SetMemoryManage(true);
  for (int i=0; i<m_wordCount; i++)
    stream.Read(m_dict(i));

  SetMemoryManage(false);

  m_searcher.Read(stream);
  m_searcher.SetDataPtr(&m_dict);

  return true;
}

bool CMStringDictionary::Write(IMWriteStream & stream)
{
  CMFileHeader header;
  header.Write(stream);

  stream.Write(m_wordCount);

  for (int i=0; i<m_wordCount; i++)
    stream.Write(m_dict(i));

  m_searcher.Write(stream);

  return true;
}


//================================================

CMFileHeader::CMFileHeader()
{
  m_version = 0;
  m_revision = 0;
}

CMFileHeader::~CMFileHeader()
{
}

long CMFileHeader::GetVersion()
{
  return m_version;
}

long CMFileHeader::GetRevision()
{
  return m_revision;
}

void CMFileHeader::SetVersion(long v)
{
  m_version = v;
}

void CMFileHeader::SetRevision(long r)
{
  m_revision = r;
}

bool CMFileHeader::Read(IMReadStream & stream)
{
  stream.Read(m_version);
  stream.Read(m_revision);
  return true;
}

bool CMFileHeader::Write(IMWriteStream & stream)
{
  stream.Write(m_version);
  stream.Write(m_revision);
  return true;
}




