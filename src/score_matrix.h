/*
 * score_matrix.h
 *
 *  Created on: 2010/04/22
 *      Author: shu
 */

#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include <stdint.h>
#include <string>
#include <vector>

class AlphabetCoder;

class ScoreMatrix {
public:
  ScoreMatrix();
  ScoreMatrix(const std::string &name, uint32_t number_letters, int match, int mismatch);
  ScoreMatrix(const std::string &name, int *matrix, uint32_t number_letters);
  virtual ~ScoreMatrix() {}

  std::string GetName() const {
    return name_;
  }

  const int* GetMatrix() const {
    return &matrix_[0];
  }

  int GetHighestValue() const {
    return highest_value_;
  }

  int GetLowestValue() const {
    return lowest_value_;
  }

  uint32_t GetNumberLetters() const {
    return number_letters_;
  }

private:
  std::string name_;
  std::vector<int> matrix_;
  int highest_value_;
  int lowest_value_;
  uint32_t number_letters_;
};

#endif /* SCORE_MATRIX_H_ */
