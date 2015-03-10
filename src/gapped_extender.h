/*
 * gapped_extender.h
 *
 *  Created on: 2012/11/12
 *      Author: shu
 */

#ifndef GAPPED_EXTENDER_H_
#define GAPPED_EXTENDER_H_

#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "sequence_type.h"
#include "score_matrix.h"
#include <stdint.h>
#include <vector>

class GappedExtender {
public:
	GappedExtender();
	virtual ~GappedExtender();

	bool ExtendOneSide(const AlphabetCoder::Code *sequence0, uint32_t sequence0_length,
	    const AlphabetCoder::Code *sequence1, const AlphabetCoder::Code sequence_delimiter,
	    bool reversed, const ScoreMatrix &score_matrix, int gap_open, int gap_extention, int cutoff,
	    int *best_score, int *best_sequence0_position,
	    int *best_sequence1_position, EditBlocks *edit_blocks);

private:
  typedef struct {
    int best;
    int best_gap;
  }DpCell;

  typedef struct {
    std::vector<uint8_t> state_array;
    uint32_t used;
  }GapState;

  static const int kInitScore = -(1 << 10);
  static const int kMaxStateArrayLength = (1 << 21);
  static const int kEditOpMask = 3;

  enum EditOpType {
    kSubstitution = 0,
    kGapInSeq0 = 1,
    kGapInSeq1 = 2,
    kExtendGapInSeq0 = 4,
    kExtendGapInSeq1 = 8
  };

  void ClearGapState();
  uint32_t GetGapStateIndex(uint32_t length);
  bool ExtendOneSideScoreOnly(const AlphabetCoder::Code *sequence0, uint32_t sequence0_length,
        const AlphabetCoder::Code *sequence1, const AlphabetCoder::Code sequence_delimiter,
        bool reversed, const ScoreMatrix &score_matrix, int gap_open, int gap_extention, int cutoff,
        int *best_score, int *best_sequence0_position,
        int *best_sequence1_position);
  void PickOptimalPath(const std::vector<uint8_t*> &edit_script, const std::vector<uint32_t> &edit_script_offset, int query_position, int db_position, int increment, EditBlocks *edit_blocks);

  std::vector<DpCell> dp_cells_;
  std::vector<GapState> gap_states_;
};

#endif /* GAPPED_EXTENDER_H_ */
