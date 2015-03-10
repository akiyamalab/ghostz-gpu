#include <assert.h>
#include "packed_alphabet_code.h"

namespace packed_alphabet_code {
size_t GetPackedAlphabetCodeSequenceLength(size_t code_sequence_length) {
	return (code_sequence_length + kNumberOfCodesInBlock - 1)
			/ kNumberOfCodesInBlock;
}

int PackAlphabetCodeSequence(const AlphabetCoder::Code *code_sequence,
		const size_t code_sequence_length,
		PackedAlphabetCode* packed_code_sequence) {
	size_t packed_code_sequence_i = 0;
	size_t sequence_i = 0;
	if (code_sequence_length > kNumberOfCodesInBlock) {
		for (; sequence_i < (code_sequence_length - kNumberOfCodesInBlock);
				sequence_i += kNumberOfCodesInBlock, ++packed_code_sequence_i) {
			PackedAlphabetCode new_packed_block = 0;
			for (size_t i = 0; i < kNumberOfCodesInBlock; ++i) {
				const PackedAlphabetCode c =
						(PackedAlphabetCode) code_sequence[sequence_i + i];

				new_packed_block |= c << i * kCodeBitSize;
			}
#if 0
			// debug //////////////////////////
			std::cout << packed_code_sequence_i << " : " << new_packed_block << std::endl;
			////////////////////////////////////
#endif
			packed_code_sequence[packed_code_sequence_i] = new_packed_block;
		}
	}
	PackedAlphabetCode new_packed_block = 0;
	size_t remain_length = code_sequence_length - sequence_i;
	for (size_t i = 0; i < remain_length; ++i) {
		const PackedAlphabetCode c =
				(PackedAlphabetCode) code_sequence[sequence_i + i];

		new_packed_block |= c << i * kCodeBitSize;
	}
#if 0
	// debug //////////////////////////
	std::cout << packed_code_sequence_i << " : " << new_packed_block << std::endl;
	////////////////////////////////////
#endif
	packed_code_sequence[packed_code_sequence_i] = new_packed_block;
	return 0;
}

}
