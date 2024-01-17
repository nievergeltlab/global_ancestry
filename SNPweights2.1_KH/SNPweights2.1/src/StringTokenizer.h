#ifndef _STRINGTOKENIZER_
#define _STRINGTOKENIZER_

#include <string>
#include <iostream>

using namespace std;

// String tokenizer class.
class StringTokenizer {

public:

   StringTokenizer(const string& s, const char* delim = NULL);
   size_t countTokens( );
   bool hasMoreTokens( ) {return(begin_ != end_);}
   void nextToken(string& s);

private:
   StringTokenizer( ) {};
   string delim_;
   string str_;
   int count_;
   int begin_;
   int end_;
};

#endif
