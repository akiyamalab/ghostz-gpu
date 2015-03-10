/*
 * score_matrix.cpp
 *
 *  Created on: 2010/04/22
 *      Author: shu
 */

#include "score_matrix.h"
#include <iostream>
#include <stdint.h>

using namespace std;

ScoreMatrix::ScoreMatrix()
:name_(""), matrix_(0), highest_value_(0), lowest_value_(0), number_letters_(0)
{}

ScoreMatrix::ScoreMatrix(const string &name, uint32_t number_letters, int match, int mismatch)
:name_(name), matrix_(0), highest_value_(match), lowest_value_(mismatch), number_letters_(number_letters)
{
  matrix_.resize(number_letters_*number_letters_);
  for (uint32_t i = 0; i < number_letters_; ++i) {
    for (uint32_t j = 0; j < number_letters_; ++j) {
      if (i == j) {
        matrix_[i*number_letters_ + j] = match;
      } else {
        matrix_[i*number_letters_ + j] = mismatch;
      }
    }
  }
}

ScoreMatrix::ScoreMatrix(const string &name, int *matrix, uint32_t number_letters)
:name_(name), matrix_(0), highest_value_(0), lowest_value_(0), number_letters_(number_letters)
{
  matrix_.insert(matrix_.begin(), matrix, matrix + number_letters*number_letters);
  for (uint32_t i = 0; i < number_letters; ++i) {
    for (uint32_t j = 0; j < number_letters; ++j) {
      int s = matrix_[number_letters_*i + j];
      if (highest_value_ < s) {
        highest_value_ = s;
      }
      if (lowest_value_ > s) {
        lowest_value_ = s;
      }
    }
  }
}
