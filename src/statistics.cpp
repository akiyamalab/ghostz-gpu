/*
 * statistics.cpp
 *
 *  Created on: 2010/10/18
 *      Author: shu
 */


#include "sequence_type.h"
#include "alphabet_coder.h"
#include "statistics.h"
#include "score_matrix.h"
#include "../ext/karlin/src/karlin.h"
#include "protein_type.h"
#include "dna_type.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include <stdexcept>

using namespace std;


float Statistics::Nominal2Normalized(int normal_score, const KarlinParameters &paramters) {
  float normalized_score;
  normalized_score = (((static_cast<float>(normal_score)*paramters.lambda) - log(paramters.K))/static_cast<float>(log(2.0)));
  return normalized_score;
}
int Statistics::Normalized2Nominal(float normalized_score, const KarlinParameters &paramters) {
  float normal_score;
  normal_score = (normalized_score*static_cast<float>(log(2.0)) + log(paramters.K))/paramters.lambda;
  return static_cast<int>(floor(normal_score));
}

int Statistics::NormalizedCutoff2NominalCutoff(float normalized_cutoff, const KarlinParameters &paramters) {
  float normal_cutoff;
  normal_cutoff = (normalized_cutoff*static_cast<float>(log(2.0)))/paramters.lambda;
  return static_cast<int>(floor(normal_cutoff));
}

double Statistics::Nominal2EValue(int normal_score, uint64_t search_space, const KarlinParameters &paramters) {
  double e_value;
  e_value = search_space*paramters.K*exp(static_cast<double>(-1.0*normal_score*paramters.lambda));
  return e_value;
}

uint64_t Statistics::CalculateSearchSpace(uint64_t query_length, uint64_t db_length, uint64_t db_number_sequences, const KarlinParameters &parameters, double alpha, double beta) {
  uint64_t length_adjustment = CalculateLengthAdjustment(query_length, db_length, db_number_sequences, parameters, alpha, beta);
  uint64_t effective_db_length = db_length -(db_number_sequences*length_adjustment);
  return effective_db_length*(query_length-length_adjustment);
}

uint64_t Statistics::CalculateLengthAdjustment(uint64_t query_length, uint64_t db_length, uint64_t db_number_sequences, const KarlinParameters &paramters, double alpha, double beta) {
  uint64_t length_adjustment = 0;
  if (typeid(*type_) == typeid(ProteinType)) {
    const int kMaxIteration = 20;
    double m = static_cast<double> (query_length);
    double n = static_cast<double> (db_length);
    double N = static_cast<double> (db_number_sequences);
    double log_k = log(paramters.K);
    double alpha_d_lambda = static_cast<double> (alpha / paramters.lambda);

    double e, e_min = 0, e_max;
    double search_space;

    double a = N;
    double mb = m * N + n;
    double c = n * m - max(m, n) / paramters.K;
    if (c < 0) {
      return 0;
    } else {
      e_max = 2 * c / (mb + sqrt((mb * mb - 4 * a * c)));
    }

    bool converged = false;
    double e_next = 0.0;
    for (int i = 0; i < kMaxIteration; ++i) {
      double e_bar;
      e = e_next;
      search_space = (m - e) * (n - N * e);
      e_bar = alpha_d_lambda * (log_k + log(search_space)) + beta;
      if (e_bar >= e) {
        e_min = e;
        if (e_bar - e_min <= 1.0) {
          converged = true;
          break;
        }
        if (e_min >= e_max) {
          break;
        }
      } else {
        e_max = e;
      }
      if (e_min <= e_bar && e_bar <= e_max) {
        e_next = e_bar;
      } else {
        e_next = (i == 0) ? e_max : (e_min + e_max) / 2;
      }
    }
    if (converged) {
      length_adjustment = static_cast<int> (e_min);
      e = ceil(e_min);

      if (e <= e_max) {
        search_space = (m - e) * (n - N * e);
        if (alpha_d_lambda * (log_k + log(search_space)) + beta > e) {
          length_adjustment = static_cast<int> (e);
        }
      }
    } else {
      length_adjustment = static_cast<int> (e_min);
    }
  } else {
    float length;
    float min_query_length = 1.0/paramters.K;

    length = 0;
    for (uint32_t i = 0; i < 5; ++i) {
      length = (paramters.K +
          log((query_length - length) *
              (db_length - db_number_sequences*length)))/paramters.H;
      if (length > query_length - min_query_length) {
        length = query_length - min_query_length;
      }
    }
    length_adjustment = (uint64_t)length;
  }
  return length_adjustment;
}

Statistics::Statistics(SequenceType &type) :
    min_regular_letter_code_(0), max_regular_letter_code_(0), min_code_(0), max_code_(0) {
 if (typeid(type) == typeid(ProteinType)) {
    type_ = SequenceTypePtr(new ProteinType);
  } else {
    type_ = SequenceTypePtr(new DnaType);
  }

  AlphabetCoder coder(type);

  min_regular_letter_code_ = coder.GetMinRegularLetterCode();
  max_regular_letter_code_ = coder.GetMaxRegularLetterCode();
  min_code_ = coder.GetMinCode();
  max_code_ = coder.GetMaxCode();
  letter_prob_.resize(max_code_ + 1);
  for (uint32_t i = 0; i <= max_code_; ++i) {
    letter_prob_[i] = 0.0;
  }

  if (typeid(type) == typeid(ProteinType)) {
    // robinson prob
    letter_prob_[coder.Encode('L')] = 90.19;
    letter_prob_[coder.Encode('A')] = 78.05;
    letter_prob_[coder.Encode('G')] = 73.77;
    letter_prob_[coder.Encode('S')] = 71.20;
    letter_prob_[coder.Encode('V')] = 64.41;
    letter_prob_[coder.Encode('E')] = 62.95;
    letter_prob_[coder.Encode('T')] = 58.41;
    letter_prob_[coder.Encode('K')] = 57.44;
    letter_prob_[coder.Encode('D')] = 53.64;
    letter_prob_[coder.Encode('P')] = 52.03;
    letter_prob_[coder.Encode('I')] = 51.42;
    letter_prob_[coder.Encode('R')] = 51.29;
    letter_prob_[coder.Encode('N')] = 44.87;
    letter_prob_[coder.Encode('Q')] = 42.64;
    letter_prob_[coder.Encode('F')] = 38.56;
    letter_prob_[coder.Encode('Y')] = 32.16;
    letter_prob_[coder.Encode('M')] = 22.43;
    letter_prob_[coder.Encode('H')] = 21.99;
    letter_prob_[coder.Encode('C')] = 19.25;
    letter_prob_[coder.Encode('W')] = 13.30;
  } else {
    // dna
    letter_prob_[coder.Encode('A')] = 255;
    letter_prob_[coder.Encode('C')] = 255;
    letter_prob_[coder.Encode('G')] = 255;
    letter_prob_[coder.Encode('T')] = 255;
  }

  Normalize(&letter_prob_[0], 1.0);

}

Statistics::~Statistics() {
}

bool Statistics::CalculateUngappedIdealKarlinParameters(const ScoreMatrix &score_matrix, KarlinParameters *paramters)  {
  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();
  vector<double> score_probabilities0(highest - lowest + 1);
  double *score_probabilities = &score_probabilities0[-lowest];
  bool ret;
  ret = CalculateScoreProbabilities(score_matrix, &letter_prob_[0], &letter_prob_[0], score_probabilities);
  if (!ret) {
    return false;
  }
  BlastKarlinBlkCalc(score_probabilities, lowest, highest, &(paramters->lambda), &(paramters->K), &(paramters->H));
  return true;
}

bool Statistics::CalculateUngappedKarlinParameters(const AlphabetCoder::Code *sequence,
    uint32_t length, const ScoreMatrix &score_matrix, KarlinParameters *paramters) {
  vector<double> query_frequency(max_code_ + 1);
  bool ret;
  ret = ComposeSequenceFrequency(sequence, length, &query_frequency[0]);
  if (!ret) {
    return false;
  }
  ret = Normalize(&query_frequency[0], 1.0);
  if (!ret) {
    return false;
  }

  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();

  vector<double> score_probabilities0(highest - lowest + 1);
  double *score_probabilities = &score_probabilities0[-lowest];
  ret = CalculateScoreProbabilities(score_matrix, &query_frequency[0], &letter_prob_[0], score_probabilities);
  if (!ret) {
    return false;
  }
  BlastKarlinBlkCalc(score_probabilities, lowest, highest, &(paramters->lambda), &(paramters->K), &(paramters->H));
  return true;
}


bool Statistics::CalculateGappedKarlinParameters(const ScoreMatrix &score_matrix, int gap_open, int gap_extension, KarlinParameters *parameters) {
  if (typeid(*type_) == typeid(DnaType)) {
    return CalculateUngappedIdealKarlinParameters(score_matrix, parameters);
  } else {
    if (score_matrix.GetName() == "BLOSUM62" && gap_open == 11 && gap_extension == 1) {
      parameters->lambda = 0.267;
      parameters->K = 0.041;
      parameters->H = 0.14;
      return true;
    } else if (score_matrix.GetName() == "PAM30" && gap_open == 9 && gap_extension == 1) {
      parameters->lambda = 0.294;
      parameters->K = 0.11;
      parameters->H = 0.61;
      return true;
    }
  }
  return false;
}

bool Statistics::CalculateAlphaBeta(const ScoreMatrix &score_matrix, int gap_open, int gap_extension, double *alpha, double *beta) {
  if (typeid(*type_) == typeid(DnaType)) {
    *alpha = 0;
    *beta = 0;
    return true;
  } else {
    if (score_matrix.GetName() == "BLOSUM62" && gap_open == 11
        && gap_extension == 1) {
      *alpha = 1.9;
      *beta = -30;
      return true;
    } else if (score_matrix.GetName() == "PAM30" && gap_open == 9 && gap_extension == 1) {
      *alpha = 0.48;
      *beta = -6;
      return true;
    }
  }
  return false;
}


bool Statistics::CalculateScoreProbabilities(const ScoreMatrix &score_matrix, double *score_frequency1, double *score_frequency2, double *score_probabilities)  {
  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();

  uint32_t number_letters = score_matrix.GetNumberLetters();

  for (int i = lowest; i <= highest; ++i) {
    score_probabilities[i] = 0.0;
  }

  const int *matrix = score_matrix.GetMatrix();
  for (AlphabetCoder::Code c1 = min_code_; c1 <= max_code_; ++c1) {
    uint32_t offset = c1*number_letters;
    for (AlphabetCoder::Code c2 = min_code_; c2 <= max_code_; ++c2) {
      int score = matrix[offset + c2];
      score_probabilities[score] += score_frequency1[c1]*score_frequency2[c2];
    }
  }

  double sum_score = 0.0;
  for (int i = lowest; i <= highest; ++i) {
    if (score_probabilities[i] > 0.0) {
      sum_score += score_probabilities[i];
    }
  }
  if (sum_score <= 0.0) {
    return false;
  }
  for (int i = lowest; i <= highest; ++i) {
    if (score_probabilities[i] > 0.0) {
      score_probabilities[i] /= sum_score;
    }
  }
  return true;
}

bool Statistics::Normalize(double *frequency, double norm)  {
  assert((frequency != NULL) && (norm > 0.0) && "invalid arguments");
  double sum = 0.0;
  for (AlphabetCoder::Code c = min_code_; c <= max_code_; ++c) {
    double p = frequency[c];
    assert((p >= 0.0) && "invalid frequencies");
    sum += p;
  }

  if (sum <= 0.0) {
    return false;
  }

  for (AlphabetCoder::Code c = min_code_; c <= max_code_; ++c) {
    frequency[c] /= sum;
    frequency[c] *= norm;
  }
  return true;
}

bool Statistics::ComposeSequenceFrequency(const AlphabetCoder::Code *sequence, uint32_t length, double *frequency)  {
  assert(sequence != NULL && frequency != NULL && "invalid arguments");
  for (AlphabetCoder::Code c = min_code_; c <= max_code_; ++c) {
    frequency[c] = 0.0;
  }

  for (uint32_t i = 0; i < length; ++i) {
    AlphabetCoder::Code c = sequence[i];
    if (min_regular_letter_code_ <= c && c <= max_regular_letter_code_) {
      ++frequency[c];
    }
  }

  for (AlphabetCoder::Code c = min_code_; c <= max_code_; ++c) {
    if (frequency[c] < 1.0) {
      ++frequency[c];
    }
  }
  return true;
}







