#include <string>
#include <map>
#include <cctype>
#include <iostream>
#include "all_codewars.h"

std::size_t codewars_mysol::duplicateCount(const std::string& in)
{
  std::map<char, int> mydict;
  for (const char &c : in)
  {
    char c1 = c;
    if (std::isalpha(c))
      c1 = std::tolower(c);
    if (mydict.find(c1) == mydict.end())
    {
      mydict[c1] = 1;
    } else
    {
      mydict[c1] += 1;
    }
  }
  size_t count_duplicates = 0;
  for (const auto& pair : mydict)
  {
    if (pair.second > 1)
      ++count_duplicates;
  }
  return count_duplicates;
}
