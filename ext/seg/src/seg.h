/*
 * seg.h
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

#ifndef _SEG_H_
#define _SEG_H_

#include <vector>
#include <string>
#include <tr1/memory>

namespace seg {

class Seg {
public:
  Seg();
  Seg(int window);
  ~Seg();
  int MaskSequence(const std::string &seq, std::string *masked_seq);
private:
  typedef struct WindowType {
    char *seq;
    int start;
    int length;
    int parent_length;
    std::vector<int> state;
    std::vector<int> composition;
    double entropy;

  } Window;

  typedef struct SeqmentType {
    int begin;
    int end;
    std::tr1::shared_ptr<SeqmentType> next;
  } Segment;
  typedef std::tr1::shared_ptr<Segment> SegmentPtr;

  static const std::string kAminoAcidCharacters;
  static const double kLn2;
  static const double kLn20;

  void Initialize(int window_inp);
  void InitAaCodeMap(void);
  void MakeSeqments(Window &seq, int offset, SegmentPtr *segs);
  void CalculateSequenceEntropies(Window &seq, std::vector<double> *H);
  int FindLeft(int i, int limit, double *H);
  int FindRight(int i, int limit, double *H);
  void Trim(Window &seq, int *leftend, int *rightend);
  double GetProbability(int *sv, int total);
  double LnPerm(int *sv, int tot);
  double LnAss(int *sv);
  void MergeSegments(Window &seq, SegmentPtr &segs);
  void MakeMaskedSequence(Window &seq, SegmentPtr &segs,
      std::string *masked_sequence);
  void AppendSequment(SegmentPtr &segs, SegmentPtr &seg);
  void InitEntropies(int window);
  double CalculateEntropy(const std::vector<int> &state);
  void BuildSequence(char *seq, int length, Window *sequence);
  bool SiftWindow1(Window *win);
  void SetState(Window *seq);
  void SetComposition(Window *win);
  void SetEntropy(Window *win);
  bool BuildWindow(Window &parent, int start, int length, Window* window);
  void DecrementState(int *state, int thisclass);
  void IncrementState(int *state, int thisclass);
  void Upper(char *string, size_t len);

  int window_;
  int down_set_;
  int up_set_;
  double low_cut_;
  double high_cut_;
  int high_len_min_;
  int max_trim_;
  int the_window_;
  std::vector<double> entropies_;
  char non_aa_char_;
  std::vector<int> aa_code_map_;
};

}

#endif
