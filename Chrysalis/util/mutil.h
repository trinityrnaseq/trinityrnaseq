// Utilities based on an open source library...


#ifndef _MUTIL_H_
#define _MUTIL_H_


#ifdef WIN32 
#define SOME_GENERIC_WINDOWS
#endif 


#ifdef SOME_GENERIC_WINDOWS
#define MDLLEXPORT  __declspec (dllexport)
#else /* SOME_GENERIC_WINDOWS */
#define MDLLEXPORT
#endif /* SOME_GENERIC_WINDOWS */

#include <string.h>
#include <stdio.h>

#ifdef SOME_GENERIC_WINDOWS
#ifdef MCL_UNICODE_ENABLE
#define MCL_USE_UNICODE_IF_DESIRED
#endif //MCL_UNICODE
#else //SOME_GENERIC_WINDOWS
#endif //SOME_GENERIC_WINDOWS



#ifdef MCL_USE_UNICODE_IF_DESIRED

#include <tchar.h>
//WIN 32

typedef TCHAR MCL_TCHAR;

/////////////////////////////////////////////////////////////////////
#define __mccl_strcpy  _tcscpy 
#define __mccl_strcat  _tcscat
#define __mccl_strlen  _tcslen
#define __mccl_strcmp  _tcscmp
#define __mccl_strncmp  strncmp
#define __mccl_printf  _tprintf
#define __mccl_sprintf _stprintf

#define __mccl_fopen  _tfopen
#define __mccl_fclose  fclose
#define __mccl_atoi   _ttoi
#define __mccl_fscanf _ftscanf
#define __mccl_fgets  _fgetts
#define __mccl_fprintf _ftprintf

#define __mccl_tolwr  _tcslwr
#define __mccl_toupr  _tcsupr

#define __mccl_atof   wcstod

#define _TMCL _T

/////////////////////////////////////////////////////////////////////
#else //MCL_USE_UNICODE_IF_DESIRED

// GENERIC
typedef char MCL_TCHAR;
#define _TMCL 


/////////////////////////////////////////////////////////////////////
// Generic string operation definitions for Unices/Linux (UNICODE? MBCS?)
#define __mccl_strcpy strcpy 
#define __mccl_strcat strcat
#define __mccl_strlen strlen
#define __mccl_strcmp strcmp
#define __mccl_strncmp strncmp
#define __mccl_printf printf
#define __mccl_sprintf sprintf

#define __mccl_fopen  fopen
#define __mccl_fclose fclose
#define __mccl_atoi   atoi

#define __mccl_fscanf fscanf
#define __mccl_fgets  fgets
#define __mccl_fprintf fprintf

#ifdef WIN32
#ifdef UNICODE
#define __mccl_tolwr  _strlwr
#define __mccl_toupr  _strupr
#else
#define __mccl_tolwr  strlwr
#define __mccl_toupr  strupr
#endif
#else //WIN32
void __mccl_tolwr(MCL_TCHAR * szText);
void __mccl_toupr(MCL_TCHAR * szText);
#endif /* WIN32 */

#define __mccl_atof   atof
/////////////////////////////////////////////////////////////////////
#endif //MCL_USE_UNICODE_IF_DESIRED


MDLLEXPORT void SetMemoryManage(bool b = true);
MDLLEXPORT void MCLSetUTF8Encode(bool b = true);
MDLLEXPORT bool MCLIsUTF8Encode();
MDLLEXPORT void MCLSetCasing(bool b = true);

class CMString
{
public:

  MDLLEXPORT CMString();                 
  MDLLEXPORT CMString(const CMString& S);
  MDLLEXPORT CMString(const MCL_TCHAR* a);   

  MDLLEXPORT ~CMString();

  MDLLEXPORT operator const MCL_TCHAR*() const;

  MDLLEXPORT bool IsEmpty() const; 
  MDLLEXPORT long length() const; 
  MDLLEXPORT long len() const {return length();} 

  // Assignment:
  MDLLEXPORT CMString&    operator=(const MCL_TCHAR*); 
  MDLLEXPORT CMString&    operator=(const CMString&);
  MDLLEXPORT CMString&    operator+=(const MCL_TCHAR*);  
  MDLLEXPORT CMString&    operator+=(const CMString& s);
  MDLLEXPORT CMString&    operator+=(const MCL_TCHAR);   


  // Indexing operators:
  //() is bounds checking in DEBUG, [] does bounds check ALWAYS
  MDLLEXPORT MCL_TCHAR&         operator[](long);
  MDLLEXPORT MCL_TCHAR&         operator()(long);
  MDLLEXPORT MCL_TCHAR          operator[](long) const;
  MDLLEXPORT MCL_TCHAR          operator()(long) const;


  MDLLEXPORT bool operator == (const CMString& s1) const;
  MDLLEXPORT bool operator == (const MCL_TCHAR * s1) const;

  MDLLEXPORT bool operator != (const CMString& s1) const;
  MDLLEXPORT bool operator != (const MCL_TCHAR * s1) const;

  MDLLEXPORT bool operator > (const CMString& s1) const;
  MDLLEXPORT bool operator > (const MCL_TCHAR * s1) const;
  
  MDLLEXPORT bool operator < (const CMString& s1) const;
  MDLLEXPORT bool operator < (const MCL_TCHAR * s1) const;

  MDLLEXPORT bool operator >= (const CMString& s1) const;
  MDLLEXPORT bool operator >= (const MCL_TCHAR * s1) const;
  
  MDLLEXPORT bool operator <= (const CMString& s1) const;
  MDLLEXPORT bool operator <= (const MCL_TCHAR * s1) const;

  MDLLEXPORT void          toLower();
  MDLLEXPORT void          toUpper();
  MDLLEXPORT void          removeLeadingChars(char c = ' ');
  MDLLEXPORT void          removeTrailingChars(char c = ' ');
  MDLLEXPORT void          removeSpaces();

private:
  //We might want to hook in a mem mgr at some point...
  MCL_TCHAR * GetMemory(long charLen) const;
  void DeleteMemory(MCL_TCHAR *) const;

  MCL_TCHAR * m_pData;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class IMStream
{
public:
    virtual ~IMStream() {}

    virtual bool Open(const CMString &) = 0;
    virtual bool Close() = 0;
    virtual bool IsOpen() = 0;
    virtual bool IsEnd() = 0;

    virtual long BytesProcessed() = 0;
};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class IMReadStream : public IMStream
{
public:

    virtual ~IMReadStream() {}

    bool Read(long & d)             { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(unsigned long & d)    { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(int & d)             { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(unsigned int & d)    { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(short & d)            { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(unsigned short & d)   { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(char & d)             { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(unsigned char & d)    { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(signed char & d)      { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(float & d)            { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(double & d)           { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(long long & d)             { return ReadSimpleType((void*)&d, sizeof(d)); }
    bool Read(unsigned long long & d)    { return ReadSimpleType((void*)&d, sizeof(d)); }

    bool Read(CMString & d)         { return ReadString(d); }
    bool ReadLine(CMString & d)         { return ReadStringLine(d); }

    bool Read(void * p, long lenInElements, long elSize = 1)    { return ReadBlob(p, lenInElements, elSize); }

    virtual IMReadStream * CloneAndOpen(const CMString &name) = 0;

protected:
    virtual bool ReadSimpleType(void * pData, long lenInBytes) = 0;
    virtual bool ReadBlob(void * pData, long lenInElements, long elSize) = 0;
    virtual bool ReadString(CMString & string) = 0;
    virtual bool ReadStringLine(CMString & string) = 0;

};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class IMWriteStream : public IMStream
{
public:
    virtual ~IMWriteStream() {}

    bool Write(const long & d)             { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const unsigned long & d)    { return WriteSimpleType((const void*)&d, sizeof(d)); }

    bool Write(const int & d)             { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const unsigned int & d)    { return WriteSimpleType((const void*)&d, sizeof(d)); }

    bool Write(const short & d)            { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const unsigned short & d)   { return WriteSimpleType((const void*)&d, sizeof(d)); }

    bool Write(const char & d)             { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const unsigned char & d)    { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const signed char & d)      { return WriteSimpleType((const void*)&d, sizeof(d)); }

    bool Write(const float & d)             { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const double & d)    { return WriteSimpleType((const void*)&d, sizeof(d)); }
 
    bool Write(const long long & d)             { return WriteSimpleType((const void*)&d, sizeof(d)); }
    bool Write(const unsigned long long & d)    { return WriteSimpleType((const void*)&d, sizeof(d)); }

 
    bool Write(const CMString & d)         { return WriteString(d); }
    bool WriteLine(const CMString & d)         { return WriteStringLine(d); }

    bool Write(const void * p, long lenInElements, long elSize = 1)    { return WriteBlob(p, lenInElements, elSize); }

    virtual IMWriteStream * CloneAndOpen(const CMString &name) = 0;

protected:
    virtual bool WriteSimpleType(const void * pData, long lenInBytes) = 0;
    virtual bool WriteBlob(const void * pData, long lenInElements, long elSize) = 0;
    virtual bool WriteString(const CMString & string) = 0;
    virtual bool WriteStringLine(const CMString & string) = 0;
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class IMLogSink
{
public:
  virtual ~IMLogSink() {}
  virtual void OnLog(const MCL_TCHAR * szText, bool bLineFeed) = 0;
};

MDLLEXPORT void MLog(const MCL_TCHAR * szText, bool bLineFeed = true);
MDLLEXPORT void MLog(const MCL_TCHAR * szText, long val, bool bLineFeed = true);
MDLLEXPORT void MLog(const MCL_TCHAR * szText, const MCL_TCHAR * szText2, bool bLineFeed = true);

MDLLEXPORT void SetLogSink(IMLogSink * p);


class CMLog;
class VerbositySetting
{
 public:
    MDLLEXPORT VerbositySetting(int verboseLevel)
        {
            setLevel(verboseLevel);
        };
    MDLLEXPORT void setLevel(int verboseLevel)
        {
            mVerboseLevel = verboseLevel;
        };
    MDLLEXPORT int getLevel() const
        {
            return mVerboseLevel;
        };
 protected:
    int mVerboseLevel;
};


#define ML_ENDL "\n"


class CMLog
{
 public:
  MDLLEXPORT CMLog();
  MDLLEXPORT ~CMLog();
  
  MDLLEXPORT CMLog & operator<<(char c);
  MDLLEXPORT CMLog & operator<<(unsigned char c);
  MDLLEXPORT CMLog & operator<<(signed char c);
  MDLLEXPORT CMLog & operator<<(const MCL_TCHAR *s);
  MDLLEXPORT CMLog & operator<<(const CMString & s);
  MDLLEXPORT CMLog & operator<<(const unsigned char *s);
  MDLLEXPORT CMLog & operator<<(const signed char *s);
  MDLLEXPORT CMLog & operator<<(const void *p);
  MDLLEXPORT CMLog & operator<<(int n);
  MDLLEXPORT CMLog & operator<<(unsigned int n);
  MDLLEXPORT CMLog & operator<<(long n);
  MDLLEXPORT CMLog & operator<<(unsigned long n);                                                            
  MDLLEXPORT CMLog & operator<<(double n);                                                            
  MDLLEXPORT CMLog & operator<<(float n);                                                            

  MDLLEXPORT CMLog & operator<<(const VerbositySetting &rVs);                                                            

  //CMLog & operator<<(ostream & (*op)(ostream&));  
 protected:
  int mVerboseLevel;

  //void print(const MCL_TCHAR * szText, bool bLineFeed);
    

};

MDLLEXPORT CMLog & mlog();


#ifdef _DEBUG
#ifndef DEBUG
#define DEBUG
#endif
#endif /*  DEBUG */

#ifdef DEBUG
#ifndef MLIST_BOUNDSCHECK
// #define MLIST_BOUNDSCHECK
#endif
#endif /*  DEBUG */

template <class Item> class TMValueVector;
template <class Item> class TMSTValueVector;

typedef TMSTValueVector<unsigned char> CMUInt8List;
typedef TMSTValueVector<signed char> CMInt8List;
typedef TMSTValueVector<char> CMCharList;  
typedef TMSTValueVector<unsigned short> CMUInt16List;
typedef TMSTValueVector<unsigned long> CMUInt32List;
typedef TMSTValueVector<short> CMInt16List;
typedef TMSTValueVector<long> CMInt32List;
typedef TMSTValueVector<int> CMIntList;
typedef TMSTValueVector<double> CMDoubleList;



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item, long size>
class TMFixedValueVector
{
public:
	TMFixedValueVector(const TMFixedValueVector &);
	TMFixedValueVector();
	~TMFixedValueVector();

	inline const Item & operator[](const long index) const;
    inline const Item & operator()(const long index) const;

	inline Item & operator[](long index);
    inline Item & operator()(long index);

	inline long length() const;

    inline Item * data();

	inline const TMFixedValueVector<Item, size> & operator = (const TMFixedValueVector<Item, size> &);
	inline const TMFixedValueVector<Item, size> & operator = (const Item item);

protected:

    inline void BoundsCheck(const long index) const;
	Item m_array[size];

private:

};

template <class Item, long size>
TMFixedValueVector<Item, size>::TMFixedValueVector()
{
}


template <class Item, long size>
TMFixedValueVector<Item, size>::TMFixedValueVector(const TMFixedValueVector & vector)
{
  *this = vector;
}

template <class Item, long size>
TMFixedValueVector<Item, size>::~TMFixedValueVector()
{
}

template <class Item, long size>
inline const TMFixedValueVector<Item, size> & TMFixedValueVector<Item, size>::operator = (const Item val)
{
  for (long i=0; i<size; i++)
    m_array[i] = val;

  return *this;
}

template <class Item, long size>
inline const TMFixedValueVector<Item, size> & TMFixedValueVector<Item, size>::operator = (const TMFixedValueVector & vector)
{
  for (long i=0; i<size; i++)
    m_array[i] = vector(i);

  return *this;
}



template <class Item, long size>
inline void TMFixedValueVector<Item, size>::BoundsCheck(const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  if (index < size)
	return;

  ThrowException(_TMCL("TMValueVector<Item>::BoundsCheck"));

#endif
}

template <class Item, long size>
inline const Item & TMFixedValueVector<Item, size>::operator[](const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_array[index];
}

template <class Item, long size>
inline const Item & TMFixedValueVector<Item, size>::operator()(const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_array[index];
}

template <class Item, long size>
inline Item & TMFixedValueVector<Item, size>::operator[](long index)
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_array[index];
}

template <class Item, long size>
inline Item & TMFixedValueVector<Item, size>::operator()(long index)
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_array[index];
}


template <class Item, long size>
inline Item * TMFixedValueVector<Item, size>::data()
{
  return m_array;
}

template <class Item, long size>
inline long TMFixedValueVector<Item, size>::length() const
{
  return size;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item, long size>
class TMIOFixedValueVector  : public TMFixedValueVector<Item, size>
{
public:
    bool Read(IMReadStream & stream) {

      for (int i=0; i<size; i++)
        this->m_array[i].Read(stream);

      return true;
    }

    bool Write(IMWriteStream & stream) {

      for (int i=0; i<size; i++)
        this->m_array[i].Write(stream);

      return true;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item, long size>
class TMSTFixedValueVector  : public TMFixedValueVector<Item, size>
{
public:
    bool Read(IMReadStream & stream) 
    {
      return stream.Read(this->m_array,
			 size,
			 sizeof(Item));
    }

    bool Write(IMWriteStream & stream) 
    {
      return stream.Write(this->m_array,
			  size,
			  sizeof(Item));
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item>
class TMValueVector
{
public:
	TMValueVector(const TMValueVector &);
	TMValueVector();
	~TMValueVector();

	inline const Item & operator[](const long index) const;
    inline const Item & operator()(const long index) const;

	inline Item & operator[](long index);
    inline Item & operator()(long index);

	inline long length() const;

	inline void reshape(long);

    inline const Item * data() const;

	inline const TMValueVector<Item> & operator = (const TMValueVector &);
protected:

    inline void BoundsCheck(const long index) const;
	long m_length;
	Item * m_pArray;

private:

};

template <class Item>
TMValueVector<Item>::TMValueVector()
{
  m_length = 0;
  m_pArray = NULL;
}


template <class Item>
TMValueVector<Item>::TMValueVector(const TMValueVector & vector)
{
  m_length = 0;
  m_pArray = NULL;
  *this = vector;
}

template <class Item>
TMValueVector<Item>::~TMValueVector()
{
  if (m_pArray != NULL)
    delete [] m_pArray;
}


template <class Item>
inline const TMValueVector<Item> & TMValueVector<Item>::operator = (const TMValueVector & vector)
{
  reshape(vector.length());
  for (long i=0; i<m_length; i++)
    m_pArray[i] = vector(i);

  return *this;
}



template <class Item>
inline void TMValueVector<Item>::BoundsCheck(const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  if (index < m_length)
	return;

  ThrowException(_TMCL("TMValueVector<Item>::BoundsCheck"));

#endif
}

template <class Item>
inline const Item & TMValueVector<Item>::operator[](const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_pArray[index];
}

template <class Item>
inline const Item & TMValueVector<Item>::operator()(const long index) const
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_pArray[index];
}

template <class Item>
inline Item & TMValueVector<Item>::operator[](long index)
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_pArray[index];
}

template <class Item>
inline Item & TMValueVector<Item>::operator()(long index)
{
#ifdef MLIST_BOUNDSCHECK
  BoundsCheck(index);
#endif
  return m_pArray[index];
}


template <class Item>
inline const Item * TMValueVector<Item>::data() const
{
  return m_pArray;
}

template <class Item>
inline long TMValueVector<Item>::length() const
{
  return m_length;
}

template <class Item>
inline void TMValueVector<Item>::reshape(long length)
{

  Item * pNewArray;
  if (length > 0)
	pNewArray = new Item[length];
  else
	pNewArray = NULL;

  long realLength = m_length;
  if (m_length > length)
	realLength = length;

  m_length = length;

  for (long i=0; i<realLength; i++) 
    pNewArray[i] = m_pArray[i];

  Item * pTempArray = m_pArray;
  m_pArray = pNewArray;

  if (pTempArray != NULL) 
    delete [] pTempArray;
}




/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item>
class TMIOValueVector : public TMValueVector<Item>
{
public:
    bool Read(IMReadStream & stream) 
    {

      long length;
      stream.Read(length);
      
      this->reshape(length);

      for (int i=0; i<length; i++)
        (this->operator())(i).Read(stream);

      return true;
    }

    bool Write(IMWriteStream & stream) 
    {

      stream.Write(this->m_length);
      
      for (int i = 0; i < this->length(); i++)
        (this->operator())(i).Write(stream);

      return true;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item>
class TMSTValueVector : public TMValueVector<Item>
{
public:
    bool Read(IMReadStream & stream) 
    {
      long length;
      stream.Read(length);
      this->reshape(length);
      
      return stream.Read(this->m_pArray,
			 length,
			 sizeof(Item));
    }

    bool Write(IMWriteStream & stream) 
    {
      stream.Write(this->m_length);
      
      return stream.Write(this->m_pArray,
			  TMValueVector<Item>::m_length,
			  sizeof(Item));
    }

};


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item, unsigned long ChunkSize = 64>
class TMPtrList
{
public:
	TMPtrList();
	~TMPtrList();

	inline const Item * operator[](const long index) const;
    inline const Item * operator()(const long index) const;

	inline Item * operator[](long index);
    inline Item * operator()(long index);

	inline long length() const;

	//References the pointer
	inline void add(Item *);
	//Makes a copy
	inline void add(const Item &);

	//References the pointer
	inline void insert(Item *, long);
    //Makes a copy
	inline void insert(const Item &, long);

	//Creates an element
	inline long add();
	inline void remove(const Item *);
	inline void remove(long index);
	inline void removeNoDelete(const Item *);
	inline void removeNoDelete(long index);

    inline void removeNoDeleteAll();
	inline void removeAll();

    TMPtrList<Item> & operator = (const TMPtrList<Item> & toCopy);

    //Deletes the object and sets it to NULL
    //- but it does NOT remove it!!
    inline void nullify(long i);

	//Copies the last object over the current; does NOT delete it
	inline void replaceWithLast(long i);

protected:

	long m_realLength;
    TMValueVector<Item *> m_list;

private:

};


template <class Item, unsigned long ChunkSize>
TMPtrList<Item, ChunkSize>::TMPtrList()
{
  m_realLength = 0;
}

template <class Item, unsigned long ChunkSize>
TMPtrList<Item, ChunkSize>::~TMPtrList()
{
  removeAll();
}



template <class Item, unsigned long ChunkSize>
inline const Item * TMPtrList<Item, ChunkSize>::operator[](const long index) const
{
  return m_list(index);
}

template <class Item, unsigned long ChunkSize>
inline const Item * TMPtrList<Item, ChunkSize>::operator()(const long index) const
{
  return m_list(index);
}


template <class Item, unsigned long ChunkSize>
inline Item * TMPtrList<Item, ChunkSize>::operator[](long index) {
  return m_list(index);
}

template <class Item, unsigned long ChunkSize>
inline Item * TMPtrList<Item, ChunkSize>::operator()(long index) 
{
  return m_list(index);
}

template <class Item, unsigned long ChunkSize>
inline TMPtrList<Item> & TMPtrList<Item, ChunkSize>::operator = (const TMPtrList<Item> & toCopy)
{
  removeAll();
  for (long i=0; i<toCopy.length(); i++) {
    Item * pNew = new Item;
    *pNew = *(toCopy(i));
    add(pNew);
  }
  return *this;
}
    

template <class Item, unsigned long ChunkSize>
inline long TMPtrList<Item, ChunkSize>::length() const
{
  return m_realLength;
}

//References the pointer
template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::add(Item * pItem)
{
  if (m_realLength == m_list.length())
	m_list.reshape(m_realLength + ChunkSize);

  m_list(m_realLength) = pItem;
  m_realLength++;
}

//Makes a copy
template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::add(const Item & item)
{
  Item * pItem = new Item(item);
  add(pItem);
}

//Creates an element
template <class Item, unsigned long ChunkSize>
inline long TMPtrList<Item, ChunkSize>::add()
{
  Item * pItem = new Item;
  add(pItem);
  return m_realLength - 1;
}

//References the pointer
template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::insert(Item * pItem, long pos)
{
  if (m_realLength == 0) {
    add(pItem);
	return;
  }

  add(m_list(m_realLength - 1));

  for (long i=m_realLength - 1; i>pos; i--) {
	m_list(i) = m_list(i-1);
  }
  m_list(pos) = pItem;
}

//Makes a copy
template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::insert(const Item & item, long pos)
{
  Item * pItem = new Item(item);
  insert(pItem, pos);
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::remove(const Item * pItem)
{
  for (long i=0; i<m_realLength; i++) {
	if (m_list(i) == pItem) {
	  remove(i);
	  break;
	}
  }
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::remove(long index)
{
  delete m_list(index);
  for (long i=index + 1; i<m_realLength; i++) {
  	m_list(i - 1) = m_list(i);
  }
  m_realLength--;
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::removeNoDelete(const Item * pItem)
{
  for (long i=0; i<m_realLength; i++) {
	if (m_list(i) == pItem) {
	  removeNoDelete(i);
	  break;
	}
  }
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::removeNoDelete(long index)
{
  for (long i=index + 1; i<m_realLength; i++) {
  	m_list(i - 1) = m_list(i);
  }
  m_realLength--;
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::removeNoDeleteAll()
{
  m_realLength = 0;
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::nullify(long i)
{
  delete m_list(i);
  m_list(i) = NULL;
}

template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::replaceWithLast(long i)
{
  m_list(i) = m_list(m_realLength - 1);
  m_realLength--;
}


template <class Item, unsigned long ChunkSize>
inline void TMPtrList<Item, ChunkSize>::removeAll()
{
  for (long i=0; i<m_realLength; i++) {
    if (m_list(i) != NULL) {
	  delete m_list(i);
    }
  }
  
  m_realLength = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item, unsigned long ChunkSize = 64>
class TMNoOwnPtrList : public TMPtrList<Item, ChunkSize>
{
public:
    ~TMNoOwnPtrList() { this->removeNoDeleteAll(); }
     inline TMNoOwnPtrList<Item> & operator = (const TMNoOwnPtrList<Item> & toCopy);
};

template <class Item, unsigned long ChunkSize>
inline TMNoOwnPtrList<Item> & TMNoOwnPtrList<Item, ChunkSize>::operator = (const TMNoOwnPtrList<Item> & toCopy)
{
  this->removeNoDeleteAll();
  //Just "borrows" the objects...
  for (long i=0; i<toCopy.length(); i++) {
    add((Item*)toCopy(i));
  }
  return *this;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class Item>
class TMIOPtrList : public TMPtrList<Item>
{
public:
    bool Read(IMReadStream & stream) 
    {

      long len;
      stream.Read(len);

      this->removeAll();

      for (int i=0; i<len; i++) {
        Item * pItem = new Item;
        add(pItem);
        pItem->Read(stream);
      }

      return true;
    }

    bool Write(IMWriteStream & stream) 
    {

      stream.Write(this->m_realLength);
      
      for (int i = 0; i < this->length(); i++)
        (this->operator())(i)->Write(stream);

      return true;
    }
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef TMPtrList<CMString> CMPtrStringList;
typedef TMValueVector<CMString> CMStringList;


class CMStringMap
{
public:
	CMStringMap(const CMString & string1, const CMString & string2) : m_string1(string1), m_string2(string2)
	{
	}

	CMStringMap(const CMString & string1) : m_string1(string1), m_string2(string1)
	{
	}

	const CMString & GetString() {return m_string1;}
	const CMString & GetMap() {return m_string2;}
	
	void SetString(const CMString & s) {m_string1 = s;}
	void SetMap(const CMString & s) {m_string2 = s;}

private:
	CMString m_string1;
	CMString m_string2;
};


typedef TMPtrList<CMStringMap> CMPtrStringMapList;

MDLLEXPORT bool Tokenize(CMPtrStringList & result, const CMString & source, char delimiter = ' ', long limit = 0x7FFFFFFF);
MDLLEXPORT bool GetNextToken(CMString & result, CMString & source, char delimiter = ' ');


class CMTokenizer
{
public:
	MDLLEXPORT CMTokenizer();
	MDLLEXPORT ~CMTokenizer();
  
	MDLLEXPORT bool Tokenize(CMPtrStringList & result, const CMString & source);
    MDLLEXPORT void AddDelimiter(const CMString & delimiter, const CMString & replacement = "");

private:
	bool IsDelimiter(const char * pBuffer, long & inc, long & mapIndex);

	CMPtrStringMapList m_map;
};




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class CMAsciiReadFileStream : public IMReadStream
{
public:
    MDLLEXPORT CMAsciiReadFileStream();
    MDLLEXPORT virtual ~CMAsciiReadFileStream();


    MDLLEXPORT virtual bool Open(const CMString &);
    MDLLEXPORT virtual bool Close();
    MDLLEXPORT virtual bool IsOpen();
    MDLLEXPORT virtual bool IsEnd();

    virtual IMReadStream * CloneAndOpen(const CMString &name);

    virtual long BytesProcessed() {return m_bytesProcessed;}

protected:
    virtual bool ReadSimpleType(void * pData, long lenInBytes);
    virtual bool ReadBlob(void * pData, long lenInElements, long elSize);
    virtual bool ReadString(CMString & string);
    virtual bool ReadStringLine(CMString & string);


private:

    FILE * m_pFile;
    CMString m_fileName;
    bool m_bIsOpen;
    bool m_bIsEof;
    long m_bytesProcessed;

};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class CMAsciiWriteFileStream : public IMWriteStream
{
public:
    MDLLEXPORT CMAsciiWriteFileStream();
    MDLLEXPORT virtual ~CMAsciiWriteFileStream();


    MDLLEXPORT virtual bool Open(const CMString &);
    MDLLEXPORT virtual bool Close();
    MDLLEXPORT virtual bool IsOpen();
    MDLLEXPORT virtual bool IsEnd();

    MDLLEXPORT virtual bool WriteLine(const CMString & line);
    virtual IMWriteStream * CloneAndOpen(const CMString &) {return NULL;}

    virtual long BytesProcessed() {return m_bytesProcessed;}

protected:
    virtual bool WriteSimpleType(const void * pData, long lenInBytes);
    virtual bool WriteBlob(const void * pData, long lenInElements, long elSize);
    virtual bool WriteString(const CMString & string);
    virtual bool WriteStringLine(const CMString & string);

private:

    FILE * m_pFile;
    CMString m_fileName;
    bool m_bIsOpen;
    bool m_bIsEof;
    long m_bytesProcessed;

};



///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class CMReadFileStream : public IMReadStream
{
public:
    MDLLEXPORT CMReadFileStream();
    MDLLEXPORT virtual ~CMReadFileStream();


    MDLLEXPORT virtual bool Open(const CMString & name);
    MDLLEXPORT virtual bool Close();
    MDLLEXPORT virtual bool IsOpen();
    MDLLEXPORT virtual bool IsEnd();
    virtual IMReadStream * CloneAndOpen(const CMString &) {return NULL;}
    virtual long BytesProcessed() {return m_bytesProcessed;}

protected:
    virtual bool ReadSimpleType(void * pData, long lenInBytes);
    virtual bool ReadBlob(void * pData, long lenInElements, long elSize);
    virtual bool ReadString(CMString & string);
    virtual bool ReadStringLine(CMString & string);


private:

    FILE * m_pFile;
    CMString m_fileName;
    bool m_bIsOpen;
    bool m_bIsEof;
    long m_bytesProcessed;

};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class CMWriteFileStream : public IMWriteStream
{
public:
    MDLLEXPORT CMWriteFileStream();
    MDLLEXPORT virtual ~CMWriteFileStream();


    MDLLEXPORT virtual bool Open(const CMString & name);
    MDLLEXPORT virtual bool Close();
    MDLLEXPORT virtual bool IsOpen();
    MDLLEXPORT virtual bool IsEnd();
    virtual IMWriteStream * CloneAndOpen(const CMString &) {return NULL;}
    virtual long BytesProcessed() {return m_bytesProcessed;}

protected:
    virtual bool WriteSimpleType(const void * pData, long lenInBytes);
    virtual bool WriteBlob(const void * pData, long lenInElements, long elSize);
    virtual bool WriteString(const CMString & string);
    virtual bool WriteStringLine(const CMString & string) {return WriteString(string);}

private:

    FILE * m_pFile;
    CMString m_fileName;
    bool m_bIsOpen;
    bool m_bIsEof;
    long m_bytesProcessed;

};

#ifdef MCL_NO_EXCEPTIONS
#define MCL_TRY
#else //MCL_NO_EXCEPTIONS
#define MCL_TRY try
#endif//MCL_NO_EXCEPTIONS

MDLLEXPORT void ThrowException();
MDLLEXPORT void ThrowException(const CMString & reason);
MDLLEXPORT void ThrowException(const CMString & comment, const CMString & reason);


class CMException
{
public:
	MDLLEXPORT CMException(const CMString & error);
	MDLLEXPORT virtual ~CMException();

	MDLLEXPORT const CMString & GetErrorText();

	MDLLEXPORT void Print();
private:
	CMString m_text;
};





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <class Item>
class TMBinaryArraySearcher
{
public:
	TMBinaryArraySearcher() : m_pSearchPtr(NULL), m_length(0), m_virtLength(0) {}
	~TMBinaryArraySearcher() {}

	inline void SetTo(const Item * pSearchList, long length);
	inline bool Search(long & index, const Item & toFind);

private:
	Item * m_pSearchPtr;

	long m_length;
	long m_virtLength;

};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <class Item, class Index>
class TMBinaryIndexSearcher
{
public:
	TMBinaryIndexSearcher() : m_pSearchPtr(NULL), m_index(), m_length(0), m_virtLength(0) {}
	~TMBinaryIndexSearcher() {}

	inline void Add(const Item & , Index);

    inline void SetDataPtr(const TMValueVector<Item> * pSearchList);
	inline bool Search(Index & index, const Item & toFind);

	TMValueVector<Index> & GetIndexTable() {return m_index;}
	long & GetLength() {return m_length;}
	long & GetVirtualLength() {return m_virtLength;}

	inline const Index IndexToIndex(const long index) const { return m_index(index); }
	inline bool Peek(const long & index, Item & toFind);
	inline bool InternalSearch(long & i, const Item & toFind);

	bool Write(IMWriteStream &);
	bool Read(IMReadStream &);

	TMBinaryIndexSearcher & operator=(const TMBinaryIndexSearcher & idx)
	{
	  m_pSearchPtr = NULL;
	  m_index = idx.m_index;
      m_length = idx.m_length;
	  m_virtLength = idx.m_virtLength;
	  return *this;
	}



protected:
	const TMValueVector<Item> * m_pSearchPtr;
	TMSTValueVector<Index> m_index;

	long m_length;
	long m_virtLength;
};




//--------------------------------------------------------------
//--------------------------------------------------------------

template <class Item>
inline void TMBinaryArraySearcher<Item>::SetTo(const Item * pSearchList, long length)
{
  m_pSearchPtr = (Item *)pSearchList;

  m_length = (long)length;
  m_virtLength = 1;
  while (m_virtLength < m_length)
	m_virtLength = m_virtLength << 1;

}

template <class Item>
inline bool TMBinaryArraySearcher<Item>::Search(long & index, const Item & toFind)
{

  long i = 0;  
  long interval = (long)m_virtLength / 2;
  index = 0;

  if (m_pSearchPtr[0] > toFind) {
    index = 0;
	return false;
  }

  if (m_pSearchPtr[m_length - 1] < toFind) {
    index = m_length;
	return false;
  }

  if (m_pSearchPtr[0] == toFind) {
    index = 0;
	return true;
  }

  while (interval > 0) {

	if (i < m_length && m_pSearchPtr[i] == toFind) {
	  index = i;
	  return true;
	}

	if (i >= m_length || m_pSearchPtr[i] > toFind) {
	  interval = interval >> 1;
	  i -= interval;
	} else {
	  i += interval;
	}
  
  }

  if (i >= m_length) {
	index = m_length;
  } else {
    index = (long)i;
  }

  return false;

} 


//-----------------------------------------------------
//-----------------------------------------------------
template <class Item, class Index>
bool TMBinaryIndexSearcher<Item, Index>::Write(IMWriteStream & stream)
{
  stream.Write(m_length);
  stream.Write(m_virtLength);
  m_index.Write(stream);

  return true;
}

template <class Item, class Index>
bool TMBinaryIndexSearcher<Item, Index>::Read(IMReadStream & stream)
{
  stream.Read(m_length);
  stream.Read(m_virtLength);
  m_index.Read(stream);
  return true;
}

template <class Item, class Index>
inline void TMBinaryIndexSearcher<Item, Index>::SetDataPtr(const TMValueVector<Item> * pSearchList)
{
  m_pSearchPtr = pSearchList;

}


template <class Item, class Index>
inline void TMBinaryIndexSearcher<Item, Index>::Add(const Item & item, Index wordIndex)
{
  if (m_length == m_virtLength) {
    m_virtLength = 1;
    while (m_virtLength < m_length + 1)
	  m_virtLength = m_virtLength << 1;
  }


  long index;
  if (m_length > 0) {
	  if (InternalSearch(index, item)) {
	  }
  } else {
	index = 0;
  }

  if (m_length + 1 >= m_index.length())
	m_index.reshape(m_index.length() + 256);

  for (long i=m_length; i>(long)index; i--) {
  	m_index(i) = m_index(i - 1);
  }
  m_index(index) = wordIndex;
  m_length++;

}


template <class Item, class Index>
inline bool TMBinaryIndexSearcher<Item, Index>::InternalSearch(long & i, const Item & toFind)
{
  if(m_length<=0)
	  return false;

  i=0;
  long interval = (long)m_virtLength / 2;

  if ((*m_pSearchPtr)[m_index(0)] > toFind) {
	i=0;
	return false;
  }

  if ((*m_pSearchPtr)[m_index(m_length - 1)] < toFind) {
	i=m_length;
	return false;
  }

  if (m_length > 0 && (*m_pSearchPtr)[m_index(0)] == toFind) {
	i=0;
	return true;
  }

  while (interval > 0) {

	if (i < m_length && (*m_pSearchPtr)[m_index(i)] == toFind) 
	  return true;
	

	if (i >= m_length || (*m_pSearchPtr)[m_index(i)] > toFind) {
	  interval = interval >> 1;
	  i -= interval;
	} else {
	  i += interval;
	}
  
  }

  if (i >= m_length) 
	i = (long)m_length;


  return false;
	
}

template <class Item, class Index>
inline bool TMBinaryIndexSearcher<Item, Index>::Peek(const long & index, Item & toFind)
{
		if(index>=m_length)
			return false;
		toFind=(*m_pSearchPtr)[m_index(index)];
		return true;
}

template <class Item, class Index>
inline bool TMBinaryIndexSearcher<Item, Index>::Search(Index & index, const Item & toFind)
{

	long i = 0;  
	bool res=InternalSearch(i,toFind);
	if(res)
		index=m_index(i);
	else
		index=(Index)i;
	return res;
} 


class CMFileHeader
{
public:
	MDLLEXPORT CMFileHeader();
	MDLLEXPORT ~CMFileHeader();

	MDLLEXPORT long GetVersion();
	MDLLEXPORT long GetRevision();

	MDLLEXPORT void SetVersion(long v);
	MDLLEXPORT void SetRevision(long r);

	MDLLEXPORT bool Read(IMReadStream & stream);
	MDLLEXPORT bool Write(IMWriteStream & stream);

private:
	long m_version;
	long m_revision;
};


#define INVALID_STRING_DICT_ID -1

class CMStringDictionary
{
public:
	MDLLEXPORT CMStringDictionary();
	MDLLEXPORT ~CMStringDictionary();

	MDLLEXPORT long GetDictID(const CMString & word);
	MDLLEXPORT long GetDictIDList(CMInt32List & result, const CMString & word);
    MDLLEXPORT const CMString & GetWord(long id);

    MDLLEXPORT long GetWordCount();
	MDLLEXPORT long GetDictIDByIndex(long index);
	MDLLEXPORT const CMString & GetWordByIndex(long index);


    //Add a word...
	MDLLEXPORT long AddWord(const CMString & word);
	MDLLEXPORT long AddWordDontCheck(const CMString & word);

    MDLLEXPORT bool Read(IMReadStream & stream);
    MDLLEXPORT bool Write(IMWriteStream & stream);
	MDLLEXPORT long IndexToIndex(long index) {return m_searcher.IndexToIndex(index);}
    MDLLEXPORT CMStringDictionary & operator = (const CMStringDictionary & dict);

private:
	CMStringList m_dict;
	long m_wordCount;
	TMBinaryIndexSearcher<CMString, long> m_searcher;

};




MDLLEXPORT const char * GetUTF8Sig();

MDLLEXPORT bool IsUTF8(const CMString & string);

MDLLEXPORT bool RemoveUTF8Sig(CMString & string);

MDLLEXPORT bool AddUTF8Sig(CMString & string);

#endif


