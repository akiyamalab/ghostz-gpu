/*
 * seg_test.cpp
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

#include <gtest/gtest.h>
#include <stdint.h>
#include <iostream>
#include "../src/seg.h"

using namespace std;

TEST(SegTest, MaskSequence) {
  uint32_t window = 12;
  seg::Seg seg(window);

  string seq;
  string masked_seq;
  seq = "SESEWWAKKPILQNFMKGAYCPLND";
  seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("SESEWWAKKPILQNFMKGAYCPLND", masked_seq);

  seq = "ARASGGQKNLYCKIL*KVLIVL*M";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("ARASGGQKNLYCKIL*KVLIVL*M", masked_seq);

  seq = "RERVVGKKTYTAKFYERCLLSSE*";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("RERVVGKKTYTAKFYERCLLSSE*", masked_seq);

  seq = "IIQRTISTFHKILQYRFFCPPLALA";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("IIQRTISTFHKILQYRFFCPPLALA", masked_seq);

  seq = "SFRGQ*APFIKFCSIGFFAHHSLS";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("SFRGQ*APFIKFCSIGFFAHHSLS", masked_seq);

  seq = "HSEDNKHLS*NFAV*VFLPTTRSR";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("HSEDNKHLS*NFAV*VFLPTTRSR", masked_seq);

  seq = "AAAAAAAAAAAAAAAAAAAAAAAAA";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "LLLLLLLLLLLLLLLLLLLLLLLL";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "CCCCCCCCCCCCCCCCCCCCCCCC";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "SSSSSSSSSSSSSSSSSSSSSSSSS";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "AAAAAAAAAAAAAAAAAAAAAAAA";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "QQQQQQQQQQQQQQQQQQQQQQQQ";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxxxxxxxxx", masked_seq);

  seq = "AAAAAAAAAAAAAAAMKGAYCPLND";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxAAAAMKGAYCPLND", masked_seq);

  seq = "LLLLLLLLLLLLLLL*KVLIVL*M";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("xxxxxxxxxxxxxxxxKxxxxxxx", masked_seq);

  seq = "CCCCCCCCCCCCCCYERCLLSSE*";
 seg.MaskSequence(seq, &masked_seq);

  EXPECT_EQ("xxxxxxxxxxxxxxYERCLLSSE*", masked_seq);

  seq = "IIQRTISTFHSSSSSSSSSSSSSSS";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("IIQRTISTFHxxxxxxxxxxxxxxx", masked_seq);

  seq = "SFRGQ*APFIAAAAAAAAAAAAAA";
 seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("SFRGQ*APFIxxxxxxxxxxxxxx", masked_seq);

  seq = "HSEDNKHLS*QQQQQQQQQQQQQQ";
  seg.MaskSequence(seq, &masked_seq);
  EXPECT_EQ("HSEDNKHLSxxxxxxxxxxxxxxx", masked_seq);

}

