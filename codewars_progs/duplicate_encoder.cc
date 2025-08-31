#include <string>
#include <cctype>
#include <map>
#include "all_codewars.h"

std::string codewars_mysol::duplicate_encoder(const std::string& word) {
  std::map<char, int> mydict;
  for (const char& c : word)
  {
      const char c1 = std::tolower(c);
      if (mydict.find(c1) == mydict.end())
      {
        mydict[c1] = 1;
      } else
      {
        mydict[c1] += 1;
      }
   }
  std::string out;
  for (const char& c : word)
  {
    const char c1 = std::tolower(c);
    if (mydict[c1] > 1)
    {
      out += ')';
    } else
    {
      out += '(';
    }
  }
  return out;
}
