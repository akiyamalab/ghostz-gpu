/*
 * reduced_alphabet_file_reader.h
 *
 *  Created on: 2012/11/28
 *      Author: shu
 */

#ifndef REDUCED_ALPHABET_FILE_READER_H_
#define REDUCED_ALPHABET_FILE_READER_H_

#include <fstream>
#include <string>
#include <vector>
#include "alphabet_coder.h"

class ReducedAlphabetFileReader {
public:
  bool Read(std::istream &in, std::vector<std::string> &alphabet_sets);
};

#endif /* REDUCED_ALPHABET_FILE_READER_H_ */
