/*
 * sequence.h
 *
 *  Created on: 2009/06/19
 *      Author: shu
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>

class Sequence {
public:
  Sequence(const std::string &name, const std::string &sequence_data)
  :name_(name), sequence_data_(sequence_data)
  {}

  virtual ~Sequence()
  {}

  std::string GetName() const {
    return name_;
  }

  std::string GetSequenceData() const {
    return sequence_data_;
  }




private:
  std::string name_;
  std::string sequence_data_;

};

#endif /* SEQUENCE_H_ */
