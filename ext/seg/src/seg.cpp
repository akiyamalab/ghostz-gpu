/*
 * seg.cpp
 *
 *   Copyright (c) 2013, Shuji Suzuki
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "seg.h"
#include "ln_factrial_table.h"
#include <iostream>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <limits.h>
#include <string.h>
#include <math.h>

using namespace std;

namespace seg {

const std::string Seg::kAminoAcidCharacters = "ACDEFGHIKLMNPQRSTVWY";
const double Seg::kLn2 = 0.69314718055994530;
const double Seg::kLn20 = 2.9957322735539909;

Seg::Seg() :
		non_aa_char_(kAminoAcidCharacters.size()), aa_code_map_(CHAR_MAX,
				non_aa_char_) {
	Initialize(12);
}
Seg::Seg(int window_inp) :
		non_aa_char_(kAminoAcidCharacters.size()), aa_code_map_(CHAR_MAX,
				non_aa_char_) {
	Initialize(window_inp);
}

void Seg::Initialize(int window_inp) {
	window_ = window_inp;
	low_cut_ = 2.2;
	high_cut_ = 2.5;

	down_set_ = 0;
	up_set_ = 1;

	high_len_min_ = 0;
	max_trim_ = 100;

	InitAaCodeMap();
	InitEntropies(window_);

}

Seg::~Seg(void) {
}

void Seg::InitAaCodeMap(void) {
	int next_code = 0;
	for (string::const_iterator it = kAminoAcidCharacters.begin();
			it != kAminoAcidCharacters.end(); ++it, ++next_code) {
		uint8_t c = (uint8_t) *it;
		aa_code_map_[c] = next_code;
	}
	return;
}

int Seg::MaskSequence(const std::string &seq, std::string *masked_seq) {
  if (seq.size() < (size_t)window_) {
		return 1;
	}
	vector<char> copyed_seq(seq.length());
	memcpy(&copyed_seq[0], &seq[0], seq.length());
	Window masked_sequence_basic_window;
	BuildSequence(&copyed_seq[0], seq.length(), &masked_sequence_basic_window);
	SegmentPtr segs;
	MakeSeqments(masked_sequence_basic_window, 0, &segs);
	MergeSegments(masked_sequence_basic_window, segs);
	MakeMaskedSequence(masked_sequence_basic_window, segs, masked_seq);
	return 0;
}

void Seg::MakeSeqments(Window &seq, int offset, SegmentPtr *segs) {
	vector<double> H(0);
	CalculateSequenceEntropies(seq, &H);
	if (H.empty()) {
		return;
	}
	int first = down_set_;
	int last = seq.length - up_set_;
	int left_limit = first;

	for (int i = first; i <= last; i++) {
		if (H[i] <= low_cut_ && H[i] != -1) {
			int left_i = FindLeft(i, left_limit, &H[0]);
			int right_i = FindRight(i, last, &H[0]);

			int left_end = left_i - down_set_;
			int right_end = right_i + up_set_ - 1;
			Window window;
			BuildWindow(seq, left_end, right_end - left_end + 1, &window);
			Trim(window, &left_end, &right_end);
			if (i + up_set_ - 1 < left_end) {
				int lend = left_i - down_set_;
				int rend = left_end - 1;
				Window left_seq;
				BuildWindow(seq, lend, rend - lend + 1, &left_seq);
				SegmentPtr left_segs;
				MakeSeqments(left_seq, offset + lend, &left_segs);
				if (left_segs) {
					if (!(*segs)) {
						*segs = left_segs;
					} else {
						AppendSequment(*segs, left_segs);
					}
				}
			}
			SegmentPtr seg = SegmentPtr(new Segment());
			seg->begin = left_end + offset;
			seg->end = right_end + offset;
			seg->next = SegmentPtr();

			if (!(*segs)) {
				*segs = seg;
			} else {
				AppendSequment(*segs, seg);
			}

			i = std::min(right_i, right_end + down_set_);
			left_limit = i + 1;
		}
	}
	return;
}

bool Seg::BuildWindow(Window &parent, int start, int length,
		Seg::Window* window) {
	if (start
			< 0|| length < 0 || start + length > parent.length || window == NULL) {
		return false;
	}

	window->start = start;
	window->length = length;
	window->seq = parent.seq + start;
	window->parent_length = parent.length;

	window->entropy = -2.0;
	window->composition.clear();
	window->state.clear();

	SetState(window);

	return window;
}

void Seg::SetState(Window *window) {
	if (window->composition.empty()) {
		SetComposition(window);
	}
	window->state.resize(21, 0);
	int nel = 0;
	for (int aa = 0; aa < 20; ++aa) {
		int c = window->composition[aa];
		if (c != 0) {
			window->state[nel] = c;
			++nel;
		}
	}
	window->state[nel] = 0;
	std::sort(window->state.begin(), window->state.begin() + nel,
			std::greater<int>());
	window->entropy = CalculateEntropy(window->state);

	return;
}

void Seg::SetComposition(Window *window) {
	window->composition.resize(20, 0);
	int *comp = &(window->composition[0]);
	char *seq = &(window->seq[0]);
	int length = window->length;

	for (int i = 0; i < length; ++i) {
		int aa = seq[i];
		if (aa_code_map_[aa] != non_aa_char_) {
			++comp[aa_code_map_[aa]];
		}
	}
	return;
}

void Seg::SetEntropy(Window *window) {
	if (window->state.empty()) {
		SetState(window);
	}
	window->entropy = CalculateEntropy(window->state);
	return;
}

void Seg::CalculateSequenceEntropies(Window &seq, std::vector<double> *H) {
	if (window_ > seq.length) {
		H->clear();
		return;
	}

	H->resize(seq.length, -1);
	Window window;
	BuildWindow(seq, 0, window_, &window);
	SetEntropy(&window);

	int first = down_set_;
	int last = seq.length - up_set_;
	for (int i = first; i <= last; i++) {
		(*H)[i] = window.entropy;
		SiftWindow1(&window);
	}
	return;
}

bool Seg::SiftWindow1(Window *window) {
	int length = window->length;
	if ((++window->start + length) > window->parent_length) {
		--window->start;
		return false;
	}
	int *comp = &(window->composition[0]);
	int sequence_first = window->seq[0];
	if (aa_code_map_[sequence_first] != non_aa_char_) {
		DecrementState(&(window->state[0]), comp[aa_code_map_[sequence_first]]);
		--comp[aa_code_map_[sequence_first]];
	}

	int sequence_last = window->seq[length];
	++window->seq;
	if (aa_code_map_[sequence_last] != non_aa_char_) {
		IncrementState(&(window->state[0]), comp[aa_code_map_[sequence_last]]);
		++comp[aa_code_map_[sequence_last]];
	}
	if (window->entropy > -2.0) {
		window->entropy = CalculateEntropy(window->state);
	}
	return true;
}

int Seg::FindLeft(int start, int limit, double *H) {
	int position = 0;
	for (position = start; position >= limit; --position) {
		if (H[position] == -1)
			break;
		if (H[position] > high_cut_)
			break;
	}

	return position + 1;
}

int Seg::FindRight(int start, int limit, double *H) {
	int positoin = 0;
	for (positoin = start; positoin <= limit; ++positoin) {
		if (H[positoin] == -1)
			break;
		if (H[positoin] > high_cut_)
			break;
	}

	return positoin - 1;
}

void Seg::Trim(Window &seq, int *left_end, int *right_end) {
	int current_left_end = 0;
	int current_rend = seq.length - 1;
	int min_length = 1;
	if ((seq.length - max_trim_) > min_length) {
		min_length = seq.length - max_trim_;
	}
	double min_prob = 1.0;
	for (int len = seq.length; len > min_length; len--) {
		Window window;
		BuildWindow(seq, 0, len, &window);
		int i = 0;

		bool shift = true;
		while (shift) {
			double prob = GetProbability(&(window.state[0]), len);
			if (prob < min_prob) {
				min_prob = prob;
				current_left_end = i;
				current_rend = len + i - 1;
			}
			shift = SiftWindow1(&window);
			++i;
		}
	}

	*left_end = *left_end + current_left_end;
	*right_end = *right_end - (seq.length - current_rend - 1);
	return;
}

double Seg::GetProbability(int *state, int total) {
	double totseq = ((double) total) * kLn20;
	double ans = LnAss(state) + LnPerm(state, total) - totseq;
	return ans;
}

double Seg::LnPerm(int *state, int total) {
	double ans = kLnFactrialTable[total];
	for (int i = 0; state[i] != 0; ++i) {
		ans -= kLnFactrialTable[state[i]];
	}
	return ans;
}

double Seg::LnAss(int *state) {
	double ans = kLnFactrialTable[20];
	if (state[0] == 0) {
		return ans;
	}
	int total = 20;
	int current_class = 1;
	int previous_state_value = state[0];
	for (int i = 1; i < 20; ++i) {
		int state_value = state[i];
		if (state_value == previous_state_value) {
			++current_class;
		} else {
			total -= current_class;
			ans -= kLnFactrialTable[current_class];
			if (state_value == 0) {
				current_class = total;
				break;
			} else {
				current_class = 1;
			}
		}
		previous_state_value = state_value;
	}
	ans -= kLnFactrialTable[current_class];

	return ans;
}

void Seg::MergeSegments(Window &seq, SegmentPtr &segs) {
	if (!segs) {
		return;
	}
	if (segs->begin < high_len_min_) {
		segs->begin = 0;
	}
	SegmentPtr seg = segs;
	SegmentPtr nextseg = seg->next;

	while (nextseg) {
		int len = nextseg->begin - seg->end - 1;
		if (len < high_len_min_) {
			seg->end = nextseg->end;
			seg->next = nextseg->next;
			nextseg = seg->next;
		} else {
			seg = nextseg;
			nextseg = seg->next;
		}
	}

	int len = seq.length - seg->end - 1;
	if (len < high_len_min_) {
		seg->end = seq.length - 1;
	}
	return;
}

void Seg::MakeMaskedSequence(Window &seq, SegmentPtr &segs,
		std::string *masked_sequence) {
	stringstream masked_sequence_stream("");
	char *sequence = &(seq.seq[0]);
	int sequence_length = seq.length;
	Upper(sequence, seq.length);

	for (SegmentPtr seg = segs; seg; seg = seg->next) {
		int begin = seg->begin;
		int end = seg->end;
		memset(sequence + begin, 'x', end - begin + 1);
	}

	for (int i = 0; i < sequence_length; ++i) {
		masked_sequence_stream << sequence[i];
	}
	*masked_sequence = masked_sequence_stream.str();
	return;
}

void Seg::AppendSequment(SegmentPtr &segs, SegmentPtr &seg) {
	SegmentPtr temp = segs;
	while (1) {
		if (temp->next == NULL) {
			temp->next = seg;
			break;
		} else {
			temp = temp->next;
		}
	}
	return;
}

void Seg::InitEntropies(int window) {
	entropies_.resize(window + 1);
	double double_window = window;
	for (int i = 1; i <= window; ++i) {
		double x = i / double_window;
		entropies_[i] = -x * log(x) / kLn2;
	}
	the_window_ = window;
}

double Seg::CalculateEntropy(const std::vector<int> &state) {
	int total = 0;
	int state_length = 0;
	for (state_length = 0; state[state_length] != 0; ++state_length) {
		total += state[state_length];
	}

	double ent = 0.0;
	if (total == the_window_) {
		for (int i = 0; i < state_length; ++i) {
			ent += entropies_[state[i]];
		}
		return ent;
	}
	if (total == 0) {
		return 0.0;
	}
	double total_reciprocal = 1. / (double) total;
	for (int i = 0; i < state_length; ++i) {
		double v = state[i];
		ent += v * log(v * total_reciprocal);
	}
	return -ent * total_reciprocal / kLn2;
}

void Seg::BuildSequence(char *seq, int length, Window *sequence) {
	sequence->length = length;
	sequence->seq = seq;
	sequence->parent_length = length;
	sequence->state.clear();
	sequence->composition.clear();

	SetState(sequence);
	return;
}

void Seg::DecrementState(int *state, int current_state_value) {
	for (int i = 0; state[i] != 0; ++i) {
		if (state[i] == current_state_value
				&& state[i + 1] < current_state_value) {
			--state[i];
			break;
		}
	}
}

void Seg::IncrementState(int *state, int current_state_value) {
	for (int i = 0;; ++i) {
		if (state[i] == current_state_value) {
			++state[i];
			break;
		}
	}
}

void Seg::Upper(char *string, size_t len) {
	for (uint32_t i = 0; i < len; ++i) {
		string[i] = toupper(string[i]);
	}
}

}
