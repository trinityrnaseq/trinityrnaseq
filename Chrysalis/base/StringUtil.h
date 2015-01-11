
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

string After(string &s, string &t);
string Before(string &s, string &t);
bool Contains(string &s, string &t);

string After(string &s, char *t);
string Before(string &s, char *t);
bool Contains(string &s, char *t);

bool ContainsAt(string &s, string &t, int at);

int PositionAfter(string &in, string& s, int startSearchAt);

inline string Stringify(int x)
{
  ostringstream out;
  out << x;
  return (out.str());
}


int Tokenize( const string &a_string,
	      vector<char> &separators,
	      vector<string> &tokens );


int Tokenize( const string &a_string,
	      vector<string> &tokens );




#endif
