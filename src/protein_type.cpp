/*
 * protein.cpp
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#include "protein_type.h"

using namespace std;

const string ProteinType::kRegularLetters = "ARNDCQEGHILKMFPSTWYV";
const string ProteinType::kAmbiguousLetters = "BJZ*";
const char ProteinType::kUnknownLetter = 'X';
