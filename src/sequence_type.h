/*
 * sequence_type.h
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#ifndef SEQUENCE_TYPE_H_
#define SEQUENCE_TYPE_H_

#ifdef __GNUC__
#include <tr1/memory>
#else
#include <string>
// for nvcc 
#define __aligned__ ignored
#include <boost/tr1/memory.hpp>
#undef __aligned__
#endif


class SequenceType {
public:
  virtual ~SequenceType();
  virtual std::string GetRegularLetters() const = 0;
  virtual std::string GetAmbiguousLetters() const = 0;
  virtual char GetUnknownLetter() const = 0;

};

typedef std::tr1::shared_ptr<SequenceType> SequenceTypePtr;

#endif /* SEQUENCE_TYPE_H_ */
