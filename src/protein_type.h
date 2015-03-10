/*
 * protein.h
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#ifndef PROTEIN_TYPE_H_
#define PROTEIN_TYPE_H_

#include "sequence_type.h"
#include <string>

class ProteinType : public SequenceType {
public:
	std::string GetRegularLetters() const {
    return kRegularLetters;
  }

	std::string GetAmbiguousLetters() const {
    return kAmbiguousLetters;
  }

  char GetUnknownLetter() const {
    return kUnknownLetter;
  }

private:
  static const char kUnknownLetter;
  static const std::string kRegularLetters;
  static const std::string kAmbiguousLetters;
};

#endif /* PROTEIN_TYPE_H_ */
