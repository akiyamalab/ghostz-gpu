/*
 * protein_sequence.h
 *
 *  Created on: 2011/04/11
 *      Author: shu
 */

#ifndef PROTEIN_SEQUENCE_H_
#define PROTEIN_SEQUENCE_H_

#include "sequence.h"
#include <string>

class ProteinSequence : public Sequence{
public:
  ProteinSequence(const std::string &name, const std::string &sequence_data)
  : Sequence(name, sequence_data)
  {}

  virtual ~ProteinSequence()
  {}
};

#endif /* PROTEIN_SEQUENCE_H_ */
