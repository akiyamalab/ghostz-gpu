/*
 * ungapped_extender.h
 *
 *  Created on: 2012/11/12
 *      Author: shu
 */

#ifndef UNGAPPED_EXTENDER_H_
#define UNGAPPED_EXTENDER_H_

#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "score_matrix.h"
#include <stdint.h>

class UngappedExtender {
public:
	static bool ExtendOneSideScoreOnly(const AlphabetCoder::Code *sequence0,
			const AlphabetCoder::Code *sequence1,
			const AlphabetCoder::Code sequence_delimiter, bool reversed,
			const ScoreMatrix &score_matrix, int cutoff, int *best_score); //__attribute__((noinline));

	static bool ExtendOneSide(const AlphabetCoder::Code *sequence0,
			const AlphabetCoder::Code *sequence1,
			const AlphabetCoder::Code sequence_delimiter, bool reversed,
			const ScoreMatrix &score_matrix, int cutoff, int *best_score,
			int *best_sequence0_position, int *best_sequence1_position,
			EditBlocks *edit_blocks); // __attribute__((noinline));

	//UngappedExtender();
	//virtual ~UngappedExtender();
};

inline bool UngappedExtender::ExtendOneSideScoreOnly(
		const AlphabetCoder::Code *sequence0,
		const AlphabetCoder::Code *sequence1,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const ScoreMatrix &score_matrix, int cutoff, int *best_score) {

	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *sm = score_matrix.GetMatrix();
	int increment = 0;
	if (reversed) {
		increment = -1;
	} else {
		increment = 1;
	}
	int threshold = -cutoff;
	int position = 0;
	int b_score = 0;
	if (sequence1[position] != sequence_delimiter
			&& sequence0[position] != sequence_delimiter) {
		int score = 0;
		do {
			 //std::cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score << std::endl;
			if (sequence1[position] == sequence_delimiter
					|| sequence0[position] == sequence_delimiter) {
				break;
			}
			score += sm[sequence1[position] * number_letters
					+ sequence0[position]];
			position += increment;
			if (score > 0) {
				do {
					b_score += score;
					 //std::cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score << std::endl;
					if (sequence1[position] == sequence_delimiter
							|| sequence0[position] == sequence_delimiter) {
						break;
					}
					score = sm[sequence1[position] * number_letters
							+ sequence0[position]];
					position += increment;
				} while (score > 0);
			}
		} while (score > threshold);
	}

	*best_score = b_score;

	return true;
}

inline bool UngappedExtender::ExtendOneSide(
		const AlphabetCoder::Code *sequence0,
		const AlphabetCoder::Code *sequence1,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const ScoreMatrix &score_matrix, int cutoff, int *best_score,
		int *best_sequence0_position, int *best_sequence1_position,
		EditBlocks *edit_blocks) {
	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *sm = score_matrix.GetMatrix();
	int increment = 0;
	if (reversed) {
		increment = -1;
	} else {
		increment = 1;
	}
	*best_sequence0_position = -increment;
	*best_sequence1_position = -increment;
	int threshold = -cutoff;
	int position = 0;
	int b_position = 0;
	int b_score = 0;
	if (sequence1[position] != sequence_delimiter
			&& sequence0[position] != sequence_delimiter) {
		int score = 0;
		do {
			//std::cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score << std::endl;
			if (sequence1[position] == sequence_delimiter
					|| sequence0[position] == sequence_delimiter) {
				break;
			}
			score += sm[sequence1[position] * number_letters
					+ sequence0[position]];
			position += increment;
			if (score > 0) {
				do {
					b_score += score;
					b_position = position;
					//std::cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score << std::endl;
					if (sequence1[position] == sequence_delimiter
							|| sequence0[position] == sequence_delimiter) {
						break;
					}
					score = sm[sequence1[position] * number_letters
							+ sequence0[position]];
					position += increment;
				} while (score > 0);
			}
		} while (score > threshold);
	}
	b_position = b_position - increment;
	*best_sequence0_position = b_position;
	*best_sequence1_position = b_position;
	*best_score = b_score;

	if (edit_blocks != NULL) {
		int edit_length = (b_position + increment) * increment;
		edit_blocks->Add(EditBlocks::kSubstitution, edit_length);
	}

	return true;
}

#endif /* UNGAPPED_EXTENDER_H_ */
