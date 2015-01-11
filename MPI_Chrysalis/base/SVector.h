#ifndef SVECTOR_H_
#define SVECTOR_H_

using namespace std;

#include <vector>
#include <iterator>
#include <assert.h>
#include <algorithm>
#include "base/ErrorHandling.h"

template <class T> 
class svec : public vector<T> 
{
 public:

  int isize() const {return (int)vector<T>::size();}
  long long lsize() const {return vector<T>::size();}

  const T & operator[] (long long i) const {
#ifndef NDEBUG
    if (i>= lsize()) {
      cout << "ERROR: i=" << i << " lsize=" << lsize() << endl;
      ThrowError("i>=lsize()", "svec");
    }
    if (i<0) {
      cout << "ERROR: i=" << i << " lsize=" << lsize() << endl;
      ThrowError("i<0", "svec");
    }
    //assert(i<lsize());
    //assert(i>=0);
#endif
    return vector<T>::operator[](i);
  }

  T & operator[] (long long i) {
#ifndef NDEBUG
    if (i>= lsize()) {
      cout << "ERROR: i=" << i << " lsize=" << lsize() << endl;
      ThrowError("i>=lsize()", "svec");
    }
    if (i<0) {
      cout << "ERROR: i=" << i << " lsize=" << lsize() << endl;
      ThrowError("i<0", "svec");
    }
    //assert(i<lsize());
    //assert(i>=0);
#endif
    return vector<T>::operator[](i);
  }



};


template <class T> 
class vec : public svec<T> 
{
};

template<class T> 
void Sort(svec<T>& v)
{
  sort(v.begin(), v.end());
}


template<class T> 
void UniqueSort(svec<T>& v)
{
  sort(v.begin(), v.end());

  long long i;
  long long k = 0;
  for (i=0; i<v.lsize(); i++) {
    v[k] = v[i];
    while (i+1<v.lsize() &&  !(v[k] < v[i+1]))
      i++;
  
    k++;    
  }
  v.resize(k);
}

template<class T> 
long long BinSearch(svec<T> & v, const T & item) 
{
  typename svec<T>::iterator iter;    
  typename svec<T>::iterator begin = v.begin();
  typename svec<T>::iterator end = v.end();
  iter = lower_bound(begin, end, item);
  long long pos = distance(begin, iter);
  
  if (pos < 0)
    return -1;

  if (pos == v.isize())
    return -1;
  
  if (v[pos] < item || item < v[pos])
    return -1;

  return pos;

}

template<class T> 
long long BinSearchFuzzy(svec<T> & v, const T & item) 
{
  typename svec<T>::iterator iter;    
  typename svec<T>::iterator begin = v.begin();
  typename svec<T>::iterator end = v.end();
  iter = lower_bound(begin, end, item);
  long pos = distance(begin, iter);
  
  if (pos < 0)
    return -1;

  return pos;

}

 


#endif //SVECTOR_H_


