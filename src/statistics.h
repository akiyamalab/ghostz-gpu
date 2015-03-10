/*
 * statistics.h
 *
 *  Created on: 2010/10/18
 *      Author: shu
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "score_matrix.h"
#include "alphabet_coder.h"
#include "sequence_type.h"
#include <vector>
#include <stdint.h>

class Statistics {
public:
  typedef struct  {
    float lambda;
    float K;
    float H;
  } KarlinParameters;

  static float Nominal2Normalized(int normal_score, const KarlinParameters &paramters);
  static int Normalized2Nominal(float normalized_score, const KarlinParameters &paramters);
  static int NormalizedCutoff2NominalCutoff(float normalized_cutoff, const KarlinParameters &paramters);

  static double Nominal2EValue(int normal_score, uint64_t search_space, const KarlinParameters &paramters);

  Statistics(SequenceType &type);
  virtual ~Statistics();

  uint64_t CalculateSearchSpace(uint64_t query_length, uint64_t db_length, uint64_t db_number_sequences, const KarlinParameters &parameters, double alpha, double beta);
  uint64_t CalculateLengthAdjustment(uint64_t query_length, uint64_t db_length, uint64_t db_number_sequences, const KarlinParameters &parameters, double alpha, double beta);


  bool CalculateUngappedIdealKarlinParameters(
      const ScoreMatrix &score_matrix, KarlinParameters *paramters) ;
  bool CalculateUngappedKarlinParameters(const AlphabetCoder::Code *query, uint32_t length,
      const ScoreMatrix &score_matrix, KarlinParameters *paramters) ;
  bool CalculateGappedKarlinParameters(const ScoreMatrix &score_matrix,
      int gap_open, int gap_extension, KarlinParameters *parameters) ;
  bool CalculateAlphaBeta(const ScoreMatrix &score_matrix, int gap_open, int gap_extension, double *alpha, double *beta) ;
private:
  bool CalculateScoreProbabilities(const ScoreMatrix &score_matrix,
      double *score_frequency1, double *score_frequency2,
      double *score_probabilities);
  bool Normalize(double *frequency, double norm);
  bool ComposeSequenceFrequency(const AlphabetCoder::Code *sequence, uint32_t length,
      double *frequency);

  SequenceTypePtr type_;
  std::vector<double> letter_prob_;
  AlphabetCoder::Code min_regular_letter_code_;
  AlphabetCoder::Code max_regular_letter_code_;
  AlphabetCoder::Code min_code_;
  AlphabetCoder::Code max_code_;

};

#endif /* STATISTICS_H_ */
