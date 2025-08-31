#include <string>
#include <vector>
#include <iostream>
#include "all_codewars.h"

std::string codewars_mysol::longestConsec(const std::vector<std::string> &strarr, int k)
{
  if ((k > 0) && ((int)(strarr.size() - k + 1) >= 0))
  {
    std::string first_longest_str = "";
    std::vector<size_t> arr_lengths = {};
    std::vector<std::string> arr_concat = {};
    size_t max_len = 0;
    for (int j = 0; j != (int)(strarr.size() - k + 1); ++j)
    {
      std::string word1 = "";
      for (int i = 0; i != k; ++i)
      {
        word1.append(strarr[j + i]);
      }
      const size_t word1_len = word1.length();
      arr_lengths.push_back(word1_len);
      arr_concat.push_back(word1);
      if (word1_len > max_len)
      {
        max_len = word1_len;
      }
    }
    for (const std::string& word1 : arr_concat)
    {
      if (word1.length() == max_len)
      {
        first_longest_str = word1;
        break;
      }
    }
    return first_longest_str;
  } else
  {
    return "";
  }
}
