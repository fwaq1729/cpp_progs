#include "all_codewars.h"
#include <limits.h>
#include "gtest/gtest.h"

// CODEWARS problems solved:
// Problem1   : Every possible sum of two digits
// Description: Given a long number, return all possible sum of two digits of it.
//   For example, 12345 : all possible sum of two digits from that number are:
//   [ 1 + 2, 1 + 3, 1 + 4, 1 + 5, 2 + 3, 2 + 4, 2 + 5, 3 + 4, 3 + 5, 4 + 5 ]
//   Therefore the result must be:
//   [ 3, 4, 5, 6, 5, 6, 7, 7, 8, 9 ]
// 
// Problem2   : Convert string to camel case
// Description: Complete the method/function so that it converts dash/underscore 
//   delimited words into camel casing.  The first word within the output should be
//   capitalized only if the original word was capitalized (known as Upper Camel Case,
//   also often referred to as Pascal case).  The next words should be always capitalized.
//   Examples:
//   "the-stealth-warrior" gets converted to "theStealthWarrior"
//   "The_Stealth_Warrior" gets converted to "TheStealthWarrior"
//   "The_Stealth-Warrior" gets converted to "TheStealthWarrior"
//   
// Problem3   : Consecutive strings
// Description: You are given an array(list) strarr of strings and an integer k.
//   Your task is to return the first longest string consisting of k consecutive strings
//   taken in the array.
//   Examples:
//   strarr = ["tree", "foling", "trashy", "blue", "abcdef", "uvwxyz"], k = 2
// 
//   Concatenate the consecutive strings of strarr by 2, we get:
// 
//   treefoling   (length 10)  concatenation of strarr[0] and strarr[1]
//   folingtrashy (       12)  concatenation of strarr[1] and strarr[2]
//   trashyblue   (       10)  concatenation of strarr[2] and strarr[3]
//   blueabcdef   (       10)  concatenation of strarr[3] and strarr[4]
//   abcdefuvwxyz (       12)  concatenation of strarr[4] and strarr[5]
// 
//   Two strings are the longest: "folingtrashy" and "abcdefuvwxyz".
//   The first that came is "folingtrashy" so 
//   longest_consec(strarr, 2) should return "folingtrashy".
// 
// Problem4    : Bouncing balls
// Description : A child is playing with a ball on the nth floor of a tall building.
//   The height of this floor above ground level, h, is known.
//   He drops the ball out of the window. The ball bounces (for example),
//   to two-thirds of its height (a bounce of 0.66).
//   His mother looks out of a window 1.5 meters from the ground.
//   How many times will the mother see the ball pass in front of her window
//   (including when it's falling and bouncing)?
//   Three conditions must be met for a valid experiment:
// 
//   Float parameter "h" in meters must be greater than 0
//   Float parameter "bounce" must be greater than 0 and less than 1
//   Float parameter "window" must be less than h.
// 
//   If all three conditions above are fulfilled, return a positive integer, otherwise return -1.
//   Note.- The ball can only be seen if the height of the rebounding ball is 
//          strictly greater than the window parameter.
//   Examples:
//   Case 1: h = 3, bounce = 0.66, window = 1.5, result is 3
//   Case 2: h = 3, bounce = 1, window = 1.5, result is -1
// 
// Problem5   : Duplicate Encoder
// Description: The goal of this exercise is to convert a string to a new string where
//   each character in the new string is "(" if that character appears only once in the
//   original string, or ")" if that character appears more than once in the original string.
//   Ignore capitalization when determining if a character is a duplicate.
//   Examples:
//   "din"      =>  "((("
//   "recede"   =>  "()()()"
//   "Success"  =>  ")())())"
//   "(( @"     =>  "))((" 
// 
// Problem6   : Counting Duplicates
// Description: Write a function that will return the count of distinct case-insensitive
//   alphabetic characters and numeric digits that occur more than once in the input string.
//   The input string can be assumed to contain only alphabets (both uppercase and lowercase)
//   and numeric digits.
//   Examples:
//   "abcde" -> 0 # no characters repeats more than once
//   "aabbcde" -> 2 # 'a' and 'b'
//   "aabBcde" -> 2 # 'a' occurs twice and 'b' twice (`b` and `B`)
//   "indivisibility" -> 1 # 'i' occurs six times
//   "Indivisibilities" -> 2 # 'i' occurs seven times and 's' occurs twice
//   "aA11" -> 2 # 'a' and '1'
//   "ABBA" -> 2 # 'A' and 'B' each occur twice
// 
// Problem7: Detect Pangram
// Description: A pangram is a sentence that contains every single letter of the alphabet
//   at least once. For example, the sentence "The quick brown fox jumps over the lazy dog"
//   is a pangram, because it uses the letters A-Z at least once (case is irrelevant).
//   Given a string, detect whether or not it is a pangram. Return True if it is,
//   False if not. Ignore numbers and punctuation.
// 
// Problem8   : Multiple of 3 or 5
// Description: If we list all the natural numbers below 10 that are multiples of 3 or 5,
//   we get 3, 5, 6 and 9. The sum of these multiples is 23.
//   Finish the solution so that it returns the sum of all the multiples of 3 or 5
//   below the number passed in.
//   Additionally, if the number is negative, return 0.
//   Note: If a number is a multiple of both 3 and 5, only count it once.


namespace {
TEST(Every_sum, test1)
{
  const std::vector<int> v1 = { 6, 7, 11 };
  EXPECT_EQ(codewars_mysol::digits(156), v1);
  const std::vector<int> v2 = { 9, 13, 17, 14, 6, 10, 7, 14, 11, 15 };
  EXPECT_EQ(codewars_mysol::digits(81596), v2);
  const std::vector<int> v3 = { 11, 8, 5, 13, 10, 7 };
  EXPECT_EQ(codewars_mysol::digits(3852), v3);
  const std::vector<int> v4 = { 5, 9, 7, 4, 5, 11, 8, 6, 3, 4, 10, 10, 7, 8, 14, 5, 6, 12, 3, 9, 10};
  EXPECT_EQ(codewars_mysol::digits(3264128), v4);
  const std::vector<int> v5 = { 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 };
  EXPECT_EQ(codewars_mysol::digits(999999), v5);
}

TEST(Strings_to_camelcase, test3)
{
  EXPECT_EQ(codewars_mysol::to_camel_case(""), "");
  EXPECT_EQ(codewars_mysol::to_camel_case("the_stealth_warrior"), "theStealthWarrior");
  EXPECT_EQ(codewars_mysol::to_camel_case("The-Stealth-Warrior"), "TheStealthWarrior");
  EXPECT_EQ(codewars_mysol::to_camel_case("A-B-C"), "ABC");
}

TEST(Consecutive_strings, test4)
{
  std::vector<std::string> arr = {"zone", "abigail", "theta", "form", "libe", "zas", "theta", "abigail"};
  EXPECT_EQ(codewars_mysol::longestConsec(arr, 2), "abigailtheta");
  arr = {"ejjjjmmtthh", "zxxuueeg", "aanlljrrrxx", "dqqqaaabbb", "oocccffuucccjjjkkkjyyyeehh"};
  EXPECT_EQ(codewars_mysol::longestConsec(arr, 1), "oocccffuucccjjjkkkjyyyeehh");
}

TEST(Bouncing_balls, test6)
{
  EXPECT_EQ(codewars_mysol::bouncingBall(3, 0.66, 1.5), 3);
  EXPECT_EQ(codewars_mysol::bouncingBall(30, 0.66, 1.5), 15);
}

TEST(Duplicate_encoder, test2)
{
   EXPECT_EQ(codewars_mysol::duplicate_encoder("din"), "(((");
   EXPECT_EQ(codewars_mysol::duplicate_encoder("recede"), "()()()");
   EXPECT_EQ(codewars_mysol::duplicate_encoder("Success"), ")())())");
   EXPECT_EQ(codewars_mysol::duplicate_encoder("CodeWarrior"), "()(((())())");
   EXPECT_EQ(codewars_mysol::duplicate_encoder("Supralapsarian"), ")()))()))))()(");
   EXPECT_EQ(codewars_mysol::duplicate_encoder("(( @"), "))((");
   EXPECT_EQ(codewars_mysol::duplicate_encoder(" ( ( )"), ")))))(");
}

TEST(Counting_duplicates, test7)
{
  EXPECT_EQ(codewars_mysol::duplicateCount("asdfghjkl54"), 0);
  EXPECT_EQ(codewars_mysol::duplicateCount("abcdeaa"), 1);
  EXPECT_EQ(codewars_mysol::duplicateCount("93917949902"), 1);
  EXPECT_EQ(codewars_mysol::duplicateCount("hhhhhhHHhhHHHHhhhhhHhH"), 1);
  EXPECT_EQ(codewars_mysol::duplicateCount("asdfghjkl55"), 1);
  EXPECT_EQ(codewars_mysol::duplicateCount("aabbcde"), 2);
  EXPECT_EQ(codewars_mysol::duplicateCount("aabBcde"), 2);
  EXPECT_EQ(codewars_mysol::duplicateCount("abcdeaB"), 2);
  EXPECT_EQ(codewars_mysol::duplicateCount("0"), 0);
  EXPECT_EQ(codewars_mysol::duplicateCount("000000000112"), 2);
  EXPECT_EQ(codewars_mysol::duplicateCount("Indivisibility"), 1);
  EXPECT_EQ(codewars_mysol::duplicateCount("Indivisibilities"), 2);
}

TEST(Is_pangram, test5)
{
 EXPECT_EQ(codewars_mysol::is_pangram("The quick, brown fox jumps over the lazy dog!"), true);
 EXPECT_EQ(codewars_mysol::is_pangram("1bcdefghijklmnopqrstuvwxyz"), false);
}

TEST(Multiples_of_3_or_5, test8)
{
  EXPECT_EQ(codewars_mysol::solution(10), 23);
}
}
