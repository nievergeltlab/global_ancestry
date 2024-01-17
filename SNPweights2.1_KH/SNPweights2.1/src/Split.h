#ifndef _SPLIT_
#define _SPLIT_

#include <string>
#include <vector>

using namespace std;

// utility class to split a line of text  
//   constructed with a string of text and empty vector
//   remove leading white-space
//   split on token passed or (default) white-space 
//   return a vector of strings

// Note: StringTokenizer separates on 1+ tabs - will not work for
// phenofiles with empty strings signifying missing data

class Split  {

  public:
    Split(string l, vector<string> *e); 
    Split(string l, vector<string> *e, const char *delim); 

};

#endif

