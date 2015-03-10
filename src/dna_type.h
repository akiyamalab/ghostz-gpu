/*
 * dna.h
 *
 *  Created on: 2010/09/28
 *      Author: shu
 */

#ifndef DNA_TYPE_H_
#define DNA_TYPE_H_

#include "sequence_type.h"
#include <string>

class DnaType : public SequenceType {
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

#endif /* DNA_TYPE_H_ */
