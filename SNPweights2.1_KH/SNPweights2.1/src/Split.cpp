#include "Split.h"
#include "StringTokenizer.h"

#include <algorithm>     // need transform

// private - remove white space
void trim(std::string& s)  {
  if ( s.empty() )
    return;
  std::string::iterator p;

  for(p=s.end(),p--; p!= s.begin() && (*p == ' ' || *p == '\t');p--)
    ;
  if ( *p != ' ' && *p != '\t' && *p != '\r' )
    p++;
  s.erase(p,s.end());

  for(p=s.begin(); p!= s.end() && (*p == ' ' || *p == '\t');p++)
    ;
  s.erase(s.begin(),p);

}

Split::Split(string l, vector<string> *e, const char* delim)  {

  string line = l;
  vector<string> *elements;

  // if the vector is not empty, empty it
  while (!e->empty())  {
    e->pop_back();
  }

  // if the line is empty, return empty vector
  if ( line.empty() )  {
    return;
  }

  // remove leading and trailing whitespace
  trim(line);

  // split the line and put in the vector
  StringTokenizer st(line,delim);
  int ntokens = (int) st.countTokens();

  for(int i=0;i<ntokens;i++)  {
    string s;
    st.nextToken(s);
    trim(s);
    e->push_back(s);
  }

  // remove trailing '\r'
  string se = e->back();
  if ( se == "\r" )  {
    e->pop_back();
  }

  return;

}

// default - split on whitespace
Split::Split(string l, vector<string> *e)  {
  Split::Split(l, e, " \t");
  return;
}
