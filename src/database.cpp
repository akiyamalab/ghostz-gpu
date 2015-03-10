/*
 * database.cpp
 *
 *  Created on: 2012/10/31
 *      Author: shu
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <string>
#include "database_chunk.h"
#include "alphabet_coder.h"
#include "reduced_alphabet_coder.h"
#include "fasta_sequence_reader.h"
#include "sequence.h"
#include "database.h"

using namespace std;

