#ifndef ERRORHANDLING_H_
#define ERRORHANDLING_H_

using namespace std;

#include <string>
#include <iostream>
#include <stdio.h>

class SException
{
 public:
  SException(const string & info, const string & type = "") {
    m_type = type;
    m_info = info;
    cout << "EXCEPTION: " << type << " " << info << endl; 
  }


  void Throw();

 private:
  string m_type;
  string m_info;

};


inline void ThrowError(const string & info, const string & type) {
  SException up(info, type);
  up.Throw();
}

 
#endif //ERRORHANDLING_H_





