/*
 * fasta_sequence_reader.cpp
 *
 *  Created on: 2009/06/19
 *      Author: shu
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <algorithm>
#include <functional>
#include "fasta_sequence_reader.h"

using namespace std;

FastaSequenceReader::FastaSequenceReader() :
    in_(NULL) {
}

FastaSequenceReader::FastaSequenceReader(std::istream &in) :
    in_(&in) {
}

bool FastaSequenceReader::Read(std::string &name, std::string &sequence) {
  name.clear();
  sequence.clear();
  std::string line;

  if (!in_->eof()) {
    // then jump to the next header line
    while (!in_->eof() && (line.length() == 0 || line.at(0) != '>')) {
      std::getline(*in_, line);
    }
    if (in_->eof()) {
      return false;
    }

    // set name
    if (line.at(line.length() - 1) == '\r') {
      line = line.substr(0, line.length() - 1);
    }
    size_t position = line.find_first_not_of("> ");
    if (position != std::string::npos) {
      name = line.substr(position);
    }

    // set sequence
    while (!in_->eof()) {
      std::getline(*in_, line);
      if (line.length() != 0) {
        if (line.at(0) != '>') {
          std::string::iterator left = std::find_if(line.begin(), line.end(),
              std::not1(std::ptr_fun<int, int>(isspace)));

          std::string::reverse_iterator right = std::find_if(line.rbegin(),
              line.rend(), std::not1(std::ptr_fun<int, int>(isspace)));
          line = std::string(left, right.base());
          sequence += line;
        } else {
          break;
        }
      }
    }

    if (in_->eof()) {
      in_->clear();
    }
    in_->seekg(-(line.size() + 1), std::ios::cur);
    return true;
  }
  return false;
}
