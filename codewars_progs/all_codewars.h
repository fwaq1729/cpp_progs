#pragma once
#include <vector>
#include <string>

namespace codewars_mysol
{
  std::vector<int> digits(int n);
  std::string duplicate_encoder(const std::string& word);
  std::string to_camel_case(std::string text);
  std::string longestConsec(const std::vector<std::string> &strarr, int k);
  bool is_pangram(const std::string& s);
  int bouncingBall(double h, double bounce, double window);
  std::size_t duplicateCount(const std::string& in);
  int solution(int number);
}
