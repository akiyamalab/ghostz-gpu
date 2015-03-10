/*
 * fasta_sequence_reader.h
 *
 *  Created on: 2009/06/19
 *      Author: shu
 */

#ifndef FASTA_SEQUENCE_READER_H_
#define FASTA_SEQUENCE_READER_H_

#include <iostream>
#include <fstream>
#include <string>

class FastaSequenceReader {
public:
  FastaSequenceReader();
  FastaSequenceReader(std::istream &in);
  void SetFileStream(std::istream &in) {
    in_ = &in;
  }
  bool IsRead() {
    if(in_ != NULL && *in_) {
      return true;
    } else {
      return false;
    }
  }

  void Seek(std::istream::pos_type position) {
    in_->clear();
    in_->seekg(position);
  }

  std::istream::pos_type Tell(){
    return in_->tellg();
  }

  bool Read(std::string &name, std::string &sequence);
private:
  std::istream *in_;

};

#endif /* FASTA_SEQUENCE_READER_H_ */
