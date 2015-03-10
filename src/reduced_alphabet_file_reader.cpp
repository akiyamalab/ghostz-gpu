/*
 * reduced_alphabet_file_reader.cpp
 *
 *  Created on: 2012/11/28
 *      Author: shu
 */

#include <limits.h>
#include <vector>
#include <string>
#include "reduced_alphabet_file_reader.h"

using namespace std;

bool ReducedAlphabetFileReader::Read(istream &in, vector<string> &alphabet_sets) {
  alphabet_sets.resize(0);
  string::size_type index;
  string line;
  std::getline(in, line);
  string delim = " ";
  string alphabet_set = "";
  index = line.find_first_of(delim);
  while (index != string::npos) {
    if (index > 0) {
      alphabet_set = line.substr(0, index);
      alphabet_sets.push_back(alphabet_set);
    }
    line = line.substr(index + 1);
    index = line.find_first_of(delim);
  }
  if (line.length() > 0) {
    alphabet_set = line;
    alphabet_sets.push_back(alphabet_set);
  }
  return true;
}
