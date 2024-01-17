#include "StringTokenizer.h"

StringTokenizer::StringTokenizer(const string& s, const char* delim) {
  str_   = s;
  count_ = -1;
  begin_ = 0;
  end_   = 0;

  if (!delim)
     delim_ = " \f\n\r\t\v";  //default to whitespace
  else
     delim_ = delim;

  // Point to the first token
  begin_ = (int) str_.find_first_not_of(delim_);
  end_ = (int) str_.find_first_of(delim_, begin_);
}


size_t StringTokenizer::countTokens( ) {
  if (count_ >= 0) // return if we've already counted
    return(count_);

  string::size_type n = 0;
  string::size_type i = 0;

  for (;;) {
     // advance to the first token
     if ((i = str_.find_first_not_of(delim_, i)) == string::npos)
        break;
     // advance to the next delimiter
     i = str_.find_first_of(delim_, i+1);
     n++;
     if (i == string::npos)
       break;
  }
  return (count_ = (int) n);
}

void StringTokenizer::nextToken(string& s) {
  if (begin_ != (int) string::npos && end_ != (int) string::npos) {
     s = str_.substr(begin_, end_-begin_);
     begin_ = (int) str_.find_first_not_of(delim_, end_);
     end_ = (int) str_.find_first_of(delim_, begin_);
  }
  else if (begin_ != (int) string::npos && end_ == (int) string::npos)  {
     s = str_.substr(begin_, str_.length( )-begin_);
     begin_ = (int) str_.find_first_not_of(delim_, end_);
  }

}
