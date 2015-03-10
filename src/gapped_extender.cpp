/*
 * gapped_extender.cpp
 *
 *  Created on: 2012/11/12
 *      Author: shu
 */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <algorithm>
#include <assert.h>
#include <limits.h>
#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "gapped_extender.h"

using namespace std;

GappedExtender::GappedExtender() :
		dp_cells_(0), gap_states_(0) {
}

GappedExtender::~GappedExtender() {
}

void GappedExtender::ClearGapState() {
	for (uint32_t i = 0; i < gap_states_.size(); ++i) {
		gap_states_[i].used = 0;
	}
}

uint32_t GappedExtender::GetGapStateIndex(uint32_t length) {
	uint32_t l = length + length / 3; // Add on about 30% so the end will get reused.
	uint32_t max_length = kMaxStateArrayLength;
	uint32_t array_length = max(max_length, l);
	if (gap_states_.empty()) {
		GapState s;
		s.used = 0;
		gap_states_.push_back(s);
		gap_states_[0].state_array.resize(array_length);
		return 0;
	} else {
		uint32_t return_i = gap_states_.size();
		for (uint32_t i = 0; i < gap_states_.size(); ++i) {
			if (l < (gap_states_[i].state_array.size() - gap_states_[i].used)) {
				return_i = i;
				break;
			} else if (gap_states_[i].used == 0) {
				gap_states_[i].state_array.resize(array_length);
				return_i = i;
				break;
			}
		}

		if (return_i == gap_states_.size()) {
			GapState s;
			s.used = 0;
			gap_states_.push_back(s);
			return_i = gap_states_.size() - 1;

			gap_states_[return_i].state_array.resize(array_length);
		}
		return return_i;
	}
}

bool GappedExtender::ExtendOneSide(const AlphabetCoder::Code *sequence0,
		uint32_t sequence0_length, const AlphabetCoder::Code *sequence1,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const ScoreMatrix &score_matrix, int gap_open, int gap_extention,
		int cutoff, int *best_score, int *best_sequence0_position,
		int *best_sequence1_position, EditBlocks *edit_blocks) {
	*best_score = 0;
	if (edit_blocks == NULL) {
		return ExtendOneSideScoreOnly(sequence0, sequence0_length, sequence1,
				sequence_delimiter, reversed, score_matrix, gap_open,
				gap_extention, cutoff, best_score, best_sequence0_position,
				best_sequence1_position);
	}
	int increment = 0;
	if (reversed) {
		increment = -1;
	} else {
		increment = 1;
	}
	*best_sequence0_position = -increment;
	*best_sequence1_position = -increment;
	edit_blocks->Clear();
	if (sequence0[0] == sequence_delimiter
			|| sequence1[0] == sequence_delimiter) {
		return 0;
	}

	uint32_t number_extra_cells = 0;
	if (gap_extention > 0 && cutoff != INT_MAX) {
		number_extra_cells = cutoff / gap_extention + 3;
	} else {
		number_extra_cells = sequence0_length + 3;
	}

	uint32_t dp_cells_length = sequence0_length + 1;
	if (dp_cells_.size() < dp_cells_length) {
		dp_cells_.resize(dp_cells_length);
	}

	ClearGapState();
	vector<uint8_t*> edit_script(100);
	vector<uint32_t> edit_script_offset(edit_script.size());

	uint32_t state_index = GetGapStateIndex(number_extra_cells);
	edit_script[0] = &(gap_states_[state_index].state_array[0]);
	edit_script_offset[0] = 0;
	uint8_t *edit_script_array = &(gap_states_[state_index].state_array[0]);

	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *score_matrix_ptr = &score_matrix.GetMatrix()[0];
	int sequence0_position = 0;
	int sequence1_position = 0;
	DpCell *score_array = NULL;
	score_array = &dp_cells_[0];

	int array_start = 0;
	int array_end = 0;
	int array_index = 0;
	int array_last_index = 0;
	int gap_init = gap_open + gap_extention;

	int score = -gap_init;
	score_array[0].best = 0;
	score_array[0].best_gap = -gap_init;
	edit_script_array[0] = kSubstitution;
	sequence0_position = 0;
	for (array_index = 1; sequence0[sequence0_position] != sequence_delimiter;
			++array_index, sequence0_position += increment) {
		if (score < -cutoff) {
			break;
		}
		score_array[array_index].best = score;
		score_array[array_index].best_gap = score - gap_init;
		score -= gap_extention;
	}
	if (sequence0[sequence0_position] != sequence_delimiter) {
		score_array[array_index].best = kInitScore;
		score_array[array_index].best_gap = kInitScore;
		array_end = array_index + 1;
	} else {
		array_end = array_index;
	}

	int max_score = 0;
	int max_score_sequence0_position = -increment;
	int max_score_sequence1_position = -increment;
	int prev_score = 0;
	int score_gap_row = 0;
	int score_gap_column = 0;

#if 0
	fprintf(stderr, "\n");
	fprintf(stderr, "     ");
	fprintf(stderr, "     ");
	for (int x = 0; sequence0[x] != sequence_delimiter; x += increment) {
		fprintf(stderr, "%5d", sequence0[x]);
	}
	fprintf(stderr, "\n");
#endif
	for (sequence1_position = 0;
			sequence1[sequence1_position] != sequence_delimiter;
			sequence1_position += increment) {
#if 0
		fprintf(stderr, "%5d", sequence1[sequence1_position - increment]);
		for (int x = 0; x < array_start; ++x) {
			fprintf(stderr, "     ");
		}
		for (int x = array_start; x < array_end; ++x) {
			fprintf(stderr, "%5d", score_array[x].best);
			//fprintf(stderr, "%3d", insertion_sequence1_row[x]);
		}
		fprintf(stderr, "\n");
#endif

		if (gap_extention > 0) {
			state_index = GetGapStateIndex(
					array_end - array_start + number_extra_cells);
		} else {
			state_index = GetGapStateIndex(number_extra_cells - array_start);
		}
		uint32_t edit_script_index = sequence1_position * increment;

		if (edit_script_index >= edit_script.size()) {
			uint32_t new_size = edit_script.size() * 2;
			edit_script.resize(new_size);
			edit_script_offset.resize(new_size);
		}
		edit_script[edit_script_index] =
				&(gap_states_[state_index].state_array[gap_states_[state_index].used
						+ 1]);
		;
		edit_script_offset[edit_script_index] = array_start;
		edit_script_array = edit_script[edit_script_index]
				- edit_script_offset[edit_script_index];

		const int *matrix_row = score_matrix_ptr
				+ sequence1[sequence1_position] * number_letters;
		sequence0_position = array_start * increment;
		prev_score = score_array[array_start].best;
		if (max_score - score_array[array_start].best_gap > cutoff) {
			score_array[array_start].best = kInitScore;
		} else {
			score_array[array_start].best = score_array[array_start].best_gap;
		}
		score_gap_row = score_array[array_start].best - gap_init;
		score_array[array_start].best_gap -= gap_extention;
		edit_script_array[array_start] = kGapInSeq0 + kExtendGapInSeq0;
		array_last_index = array_start;
		score = kInitScore;
		for (array_index = array_start + 1; array_index < array_end;
				++array_index, sequence0_position += increment) {
			score = prev_score + matrix_row[sequence0[sequence0_position]];
			prev_score = score_array[array_index].best;
			score_gap_column = score_array[array_index].best_gap;

			uint8_t op = kSubstitution;
			uint8_t op_col = kExtendGapInSeq0;
			uint8_t op_row = kExtendGapInSeq1;

			if (score < score_gap_column) {
				op = kGapInSeq0;
				score = score_gap_column;
			}
			if (score < score_gap_row) {
				op = kGapInSeq1;
				score = score_gap_row;
			}
			if (max_score - score > cutoff) {
				if (array_start + 1 == array_index) {
					++array_start;
				}
				score_array[array_index].best = kInitScore;
				score_array[array_index].best_gap = kInitScore;
				score_gap_row = kInitScore;
			} else {
				array_last_index = array_index;
				if (score > max_score) {
					max_score = score;
					max_score_sequence0_position = sequence0_position;
					max_score_sequence1_position = sequence1_position;
				}

				score_gap_column -= gap_extention;
				score_gap_row -= gap_extention;
				if (score_gap_column < score - gap_init) {
					score_array[array_index].best_gap = score - gap_init;
				} else {
					score_array[array_index].best_gap = score_gap_column;
					op += op_col;
				}

				if (score_gap_row < score - gap_init) {
					score_gap_row = score - gap_init;
				} else {
					op += op_row;
				}
				edit_script_array[array_index] = op;
				score_array[array_index].best = score;
			}
		}

		if (array_start + 1 == array_end) {
			break;
		}

		if (array_last_index < array_end - 1) {
			array_end = array_last_index + 1;
		} else {
			while (score_gap_row >= (max_score - cutoff)
					&& sequence0[sequence0_position] != sequence_delimiter) {
				score_array[array_end].best = score_gap_row;
				score_array[array_end].best_gap = score_gap_row - gap_init;
				edit_script_array[array_end] = kGapInSeq0 + kExtendGapInSeq0;
				score_gap_row -= gap_extention;
				++array_end;
				sequence0_position += increment;
			}
		}
		if (sequence0[sequence0_position] != sequence_delimiter) {
			score_array[array_end].best = kInitScore;
			score_array[array_end].best_gap = kInitScore;
			++array_end;
		}

		gap_states_[state_index].used += array_end
				- edit_script_offset[edit_script_index];
	}

	*best_score = max_score;
	*best_sequence0_position = max_score_sequence0_position;
	*best_sequence1_position = max_score_sequence1_position;
	PickOptimalPath(edit_script, edit_script_offset, *best_sequence0_position,
			*best_sequence1_position, increment, edit_blocks);
	return true;
}

bool GappedExtender::ExtendOneSideScoreOnly(
		const AlphabetCoder::Code *sequence0, uint32_t sequence0_length,
		const AlphabetCoder::Code *sequence1,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const ScoreMatrix &score_matrix, int gap_open, int gap_extention,
		int cutoff, int *best_score, int *best_sequence0_position,
		int *best_sequence1_position) {
	*best_score = 0;
	int increment = 0;
	if (reversed) {
		increment = -1;
	} else {
		increment = 1;
	}
	*best_sequence0_position = -increment;
	*best_sequence1_position = -increment;
	if (sequence0[0] == sequence_delimiter
			|| sequence1[0] == sequence_delimiter) {
		return 0;
	}

	uint32_t dp_cells_length = sequence0_length + 1;
	if (dp_cells_.size() < dp_cells_length) {
		dp_cells_.resize(dp_cells_length);
	}

	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *score_matrix_ptr = &score_matrix.GetMatrix()[0];
	int sequence0_position = 0;
	int sequence1_position = 0;
	DpCell *score_array = NULL;
	score_array = &dp_cells_[0];

	int array_start = 0;
	int array_end = 0;
	int array_index = 0;
	int array_last_index = 0;
	int gap_init = gap_open + gap_extention;

	int score = -gap_init;
	score_array[0].best = 0;
	score_array[0].best_gap = -gap_init;
	sequence0_position = 0;
	for (array_index = 1; sequence0[sequence0_position] != sequence_delimiter;
			++array_index, sequence0_position += increment) {
		if (score < -cutoff) {
			break;
		}
		score_array[array_index].best = score;
		score_array[array_index].best_gap = score - gap_init;
		score -= gap_extention;
	}
	if (sequence0[sequence0_position] != sequence_delimiter) {
		score_array[array_index].best = kInitScore;
		score_array[array_index].best_gap = kInitScore;
		array_end = array_index + 1;
	} else {
		array_end = array_index;
	}
	int max_score = 0;
	int max_score_sequence0_position = -increment;
	int max_score_sequence1_position = -increment;
	int prev_score = 0;
	int score_gap_row = 0;
	int score_gap_column = 0;

#if 0
	fprintf(stderr, "\n");
	fprintf(stderr, "     ");
	fprintf(stderr, "     ");
	for (int x = 0; sequence0[x] != sequence_delimiter; x += increment) {
		fprintf(stderr, "%5d", sequence0[x]);
	}
	fprintf(stderr, "\n");
#endif

	for (sequence1_position = 0;
			sequence1[sequence1_position] != sequence_delimiter;
			sequence1_position += increment) {
#if 0
		fprintf(stderr, "%5d", sequence1[sequence1_position - increment]);
		for (int x = 0; x < array_start; ++x) {
			fprintf(stderr, "     ");
		}
		for (int x = array_start; x < array_end; ++x) {
			fprintf(stderr, "%5d", score_array[x].best);
			//fprintf(stderr, "%3d", insertion_sequence1_row[x]);
		}
		fprintf(stderr, "\n");
#endif

		const int *matrix_row = score_matrix_ptr
				+ sequence1[sequence1_position] * number_letters;
		sequence0_position = array_start * increment;
		prev_score = score_array[array_start].best;
		if (max_score - score_array[array_start].best_gap > cutoff) {
			score_array[array_start].best = kInitScore;
		} else {
			score_array[array_start].best = score_array[array_start].best_gap;
		}
		score_gap_row = score_array[array_start].best - gap_init;
		score_array[array_start].best_gap -= gap_extention;
		array_last_index = array_start;
		score = kInitScore;
		for (array_index = array_start + 1; array_index < array_end;
				++array_index, sequence0_position += increment) {
			score = prev_score + matrix_row[sequence0[sequence0_position]];
			prev_score = score_array[array_index].best;
			score_gap_column = score_array[array_index].best_gap;

			if (score < score_gap_column) {
				score = score_gap_column;
			}
			if (score < score_gap_row) {
				score = score_gap_row;
			}
			if (max_score - score > cutoff) {
				if (array_start + 1 == array_index) {
					++array_start;
				}
				score_array[array_index].best = kInitScore;
				score_array[array_index].best_gap = kInitScore;
				score_gap_row = kInitScore;
			} else {
				array_last_index = array_index;
				if (score > max_score) {
					max_score = score;
					max_score_sequence0_position = sequence0_position;
					max_score_sequence1_position = sequence1_position;
				}

				score_array[array_index].best_gap = max(score - gap_init,
						score_gap_column - gap_extention);
				score_gap_row = max(score - gap_init,
						score_gap_row - gap_extention);
				score_array[array_index].best = score;
			}
		}

		if (array_start + 1 == array_end) {
			break;
		}

		if (array_last_index < array_end - 1) {
			array_end = array_last_index + 1;
		} else {
			while (score_gap_row >= (max_score - cutoff)
					&& sequence0[sequence0_position] != sequence_delimiter) {
				score_array[array_end].best = score_gap_row;
				score_array[array_end].best_gap = score_gap_row - gap_init;
				score_gap_row -= gap_extention;
				++array_end;
				sequence0_position += increment;
			}
		}
		if (sequence0[sequence0_position] != sequence_delimiter) {
			score_array[array_end].best = kInitScore;
			score_array[array_end].best_gap = kInitScore;
			++array_end;
		}
	}
	*best_score = max_score;
	*best_sequence0_position = max_score_sequence0_position;
	*best_sequence1_position = max_score_sequence1_position;
	return true;
}

void GappedExtender::PickOptimalPath(const vector<uint8_t*> &edit_script,
		const vector<uint32_t> &edit_script_offset, int sequence0_position,
		int sequence1_position, int increment, EditBlocks *edit_blocks) {
	uint8_t op = kSubstitution;
	uint8_t next_op = kSubstitution;

	while (sequence1_position * increment >= 0
			&& sequence0_position * increment >= 0) {
		uint32_t edit_script_index = sequence1_position * increment;
		next_op =
				edit_script[edit_script_index][(sequence0_position * increment)
						- edit_script_offset[edit_script_index] + 1];
		switch (op) {
		case kGapInSeq0:
			op = next_op & kEditOpMask;
			if (next_op & kExtendGapInSeq0) {
				op = kGapInSeq0;
			}
			break;

		case kGapInSeq1:
			op = next_op & kEditOpMask;
			if (next_op & kExtendGapInSeq1) {
				op = kGapInSeq1;
			}
			break;

		default:
			op = next_op & kEditOpMask;
			break;
		}

		switch (op) {
		case kGapInSeq0:
			edit_blocks->Add(EditBlocks::kGapInSeq0, 1);
			sequence1_position = sequence1_position - increment;
			break;
		case kGapInSeq1:
			edit_blocks->Add(EditBlocks::kGapInSeq1, 1);
			sequence0_position = sequence0_position - increment;
			break;
		default:
			edit_blocks->Add(EditBlocks::kSubstitution, 1);
			sequence1_position = sequence1_position - increment;
			sequence0_position = sequence0_position - increment;
			break;
		}
	}
	if (sequence1_position * increment >= 0) {
		edit_blocks->Add(EditBlocks::kGapInSeq0,
				sequence1_position * increment + 1);
	} else if (sequence0_position * increment >= 0) {
		edit_blocks->Add(EditBlocks::kGapInSeq1,
				sequence0_position * increment + 1);
	}
	if (increment > 0) {
		edit_blocks->Reverse();
	}
}

