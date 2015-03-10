/*
 * database_chunk.h
 *
 *  Created on: 2012/11/08
 *      Author: shu
 */

#ifndef DATABASE_CHUNK_H_
#define DATABASE_CHUNK_H_

#include <fstream>
#include <string>
#include <stdint.h>
#include <tr1/memory>
#include <boost/thread/mutex.hpp>
#include "alphabet_coder.h"
#include "seed_searcher.h"
#include "seed_searcher_gpu.h"

template<typename TSeedSearcher>
class DatabaseChunk {
public:
	typedef TSeedSearcher SeedSearcher;
	typedef typename SeedSearcher::DatabaseParameters SeedSearcherParameters;
	typedef typename SeedSearcherParameters::BuildParameters SeedSearcherParametersBuildParameters;
	DatabaseChunk();
	DatabaseChunk(std::vector<Sequence> &sequences, AlphabetCoder &coder,
			AlphabetCoder::Code sequence_delimiter,
			SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters);
	DatabaseChunk(std::string filename_prefix);

	bool Build(std::vector<Sequence> &sequences, AlphabetCoder &coder,
			AlphabetCoder::Code sequence_delimiter,
			SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters);
	bool Clear();
	uint32_t GetNumberSequences() {
		return number_sequences_;
	}

	uint32_t GetConcatenatedSequenceLength() {
		return concatenated_sequences_length_;
	}

	std::string GetName(uint32_t id) {
		boost::mutex::scoped_lock look(names_loading_mutex_);
		{
			if (names_.size() == 0) {
				LoadNames(filename_prefix_);
			}
		}
		return names_[id];
	}

	uint32_t GetOffsets(uint32_t id) {
		boost::mutex::scoped_lock look(offsets_loading_mutex_);
		{
			if (offsets_.size() == 0) {
				LoadOffsets(filename_prefix_);
			}
		}
		return offsets_[id];
	}

	AlphabetCoder::Code *GetConcatenatedSequence() {
		boost::mutex::scoped_lock look(concatenated_sequence_loading_mutex_);
		{
			if (concatenated_sequence_.size() == 0) {
				LoadConcatenatedSequence(filename_prefix_);
			}
		}
		return &concatenated_sequence_[0];
	}

	typename SeedSearcherParameters::CommonParameters &GetSeedSearcherCommonParameters() {
		boost::mutex::scoped_lock look(
				seed_searcher_common_parameters_loading_mutex_);
		{
			if (!setted_seed_searcher_common_parameters_) {
				LoadSeedSearcherCommonParameters(filename_prefix_);
			}
		}
		return seed_searcher_common_parameters_;
	}

	SeedSearcherParameters &GetSeedSearcherParameters() {
		boost::mutex::scoped_lock look(seed_searcher_parameters_loading_mutex_);
		{
			if (!setted_seed_searcher_parameters_) {
				LoadSeedSearcherParameters(filename_prefix_);
			}
		}
		return seed_searcher_parameters_;
	}

	uint32_t GetId(uint32_t position);

	bool Load(std::string filename_prefix);
	bool Save(std::string filename_prefix);
private:
	DatabaseChunk& operator =(const DatabaseChunk& rhs);
	DatabaseChunk(const DatabaseChunk& rhs);

	static std::string GetInformationFileName(std::string filename_prefix) {
		return filename_prefix + ".inf";
	}
	static std::string GetOffsetsFileName(std::string filename_prefix) {
		return filename_prefix + ".off";
	}
	static std::string GetNamesFileName(std::string filename_prefix) {
		return filename_prefix + ".nam";
	}
	static std::string GetSequencesFileName(std::string filename_prefix) {
		return filename_prefix + ".seq";
	}

	static std::string GetSeedSearcherCommonParameters(
			std::string filename_prefix) {
		return filename_prefix + ".scp";
	}

	static std::string GetSeedSearcherParameters(std::string filename_prefix) {
		return filename_prefix + ".sdp";
	}

	void EncodeSequences(AlphabetCoder &coder, std::vector<Sequence> &sequences,
			AlphabetCoder::Code sequence_delimiter);

	bool LoadInfomation(std::string filename_prefix);
	bool LoadOffsets(std::string filename_prefix);
	bool LoadNames(std::string filename_prefix);
	bool LoadConcatenatedSequence(std::string filename_prefix);
	bool LoadSeedSearcherCommonParameters(std::string filename_prefix);
	bool LoadSeedSearcherParameters(std::string filename_prefix);

	bool SaveInfomation(std::string filename_prefix);
	bool SaveOffsets(std::string filename_prefix);
	bool SaveNames(std::string filename_prefix);
	bool SaveConcatenatedSequence(std::string filename_prefix);
	bool SaveSeedSearcherCommonParameters(std::string filename_prefix);
	bool SaveSeedSearcherParameters(std::string filename_prefix);

	boost::mutex names_loading_mutex_;
	boost::mutex offsets_loading_mutex_;
	boost::mutex concatenated_sequence_loading_mutex_;
	boost::mutex seed_searcher_common_parameters_loading_mutex_;
	boost::mutex seed_searcher_parameters_loading_mutex_;

	bool building_;
	bool setted_seed_searcher_common_parameters_;
	bool setted_seed_searcher_parameters_;
	std::string filename_prefix_;
	uint32_t number_sequences_;
	uint32_t concatenated_sequences_length_;
	std::vector<std::string> names_;
	std::vector<uint32_t> offsets_;
	std::vector<AlphabetCoder::Code> concatenated_sequence_;
	typename SeedSearcherParameters::CommonParameters seed_searcher_common_parameters_;
	SeedSearcherParameters seed_searcher_parameters_;
};

template<typename TSeedSearcher>
DatabaseChunk<TSeedSearcher>::DatabaseChunk() :
		building_(false), setted_seed_searcher_common_parameters_(false), setted_seed_searcher_parameters_(
				false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(
				0), names_(0), offsets_(0), concatenated_sequence_(0), seed_searcher_parameters_() {
}

template<typename TSeedSearcher>
DatabaseChunk<TSeedSearcher>::DatabaseChunk(std::vector<Sequence> &sequences,
		AlphabetCoder &coder, AlphabetCoder::Code sequence_delimiter,
		SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters) :
		building_(false), setted_seed_searcher_common_parameters_(false), setted_seed_searcher_parameters_(
				false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(
				0), names_(0), offsets_(0), concatenated_sequence_(0), seed_searcher_parameters_() {
	Build(sequences, coder, sequence_delimiter,
			seed_search_parameters_build_parameters);
}

template<typename TSeedSearcher>
DatabaseChunk<TSeedSearcher>::DatabaseChunk(std::string filename_prefix) :
		building_(false), setted_seed_searcher_common_parameters_(false), setted_seed_searcher_parameters_(
				false), filename_prefix_(""), number_sequences_(0), concatenated_sequences_length_(
				0), names_(0), offsets_(0), concatenated_sequence_(0), seed_searcher_parameters_() {
	Load(filename_prefix);
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::Build(std::vector<Sequence> &sequences,
		AlphabetCoder &coder, AlphabetCoder::Code sequence_delimiter,
		SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters) {
	number_sequences_ = sequences.size();
	names_.resize(sequences.size());
	offsets_.resize(sequences.size() + 1);
	concatenated_sequences_length_ = 1;
	for (unsigned int i = 0; i < sequences.size(); ++i) {
		concatenated_sequences_length_ += sequences[i].GetSequenceData().size()
				+ 1;
		names_[i] = sequences[i].GetName();
	}
	EncodeSequences(coder, sequences, sequence_delimiter);
	seed_searcher_parameters_.Build(&concatenated_sequence_[0],
			concatenated_sequence_.size(),
			seed_search_parameters_build_parameters);
	setted_seed_searcher_parameters_ = true;
	seed_searcher_parameters_.BuildCommonParameters(
			&seed_searcher_common_parameters_);
	setted_seed_searcher_common_parameters_ = true;
	building_ = true;
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::Clear() {
	number_sequences_ = 0;
	concatenated_sequences_length_ = 0;
	names_.clear();
	offsets_.clear();
	concatenated_sequence_.clear();
	setted_seed_searcher_common_parameters_ = false;
	setted_seed_searcher_parameters_ = false;
	seed_searcher_parameters_ = SeedSearcherParameters();
	building_ = false;
	return true;
}

template<typename TSeedSearcher>
uint32_t DatabaseChunk<TSeedSearcher>::GetId(uint32_t position) {
	boost::mutex::scoped_lock look(offsets_loading_mutex_);
	{
		if (offsets_.size() == 0) {
			LoadOffsets(filename_prefix_);
		}
	}
	if (offsets_[number_sequences_ - 1] <= position
			&& position < concatenated_sequences_length_) {
		return number_sequences_ - 1;
	}

	uint32_t left = 0;
	uint32_t right = number_sequences_ - 2;
	uint32_t mid;
	while (left <= right) {
		mid = (left + right) / 2;
		if (offsets_[mid] <= position && position < offsets_[mid + 1]) {
			return mid;
		} else if (offsets_[mid] < position) {
			left = mid + 1;
		} else {
			right = mid - 1;
		}
	}
	return UINT_MAX;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::Load(std::string filename_prefix) {
	filename_prefix_ = filename_prefix;
	names_.clear();
	offsets_.clear();
	concatenated_sequence_.clear();
	setted_seed_searcher_common_parameters_ = false;
	setted_seed_searcher_parameters_ = false;
	return LoadInfomation(filename_prefix);
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::Save(std::string filename_prefix) {
	bool ret = false;
	if (building_ == true) {
		ret = true;
		ret &= SaveInfomation(filename_prefix);
		ret &= SaveNames(filename_prefix);
		ret &= SaveOffsets(filename_prefix);
		ret &= SaveConcatenatedSequence(filename_prefix);
		ret &= SaveSeedSearcherCommonParameters(filename_prefix);
		ret &= SaveSeedSearcherParameters(filename_prefix);
	}
	return ret;
}

template<typename TSeedSearcher>
void DatabaseChunk<TSeedSearcher>::EncodeSequences(AlphabetCoder &coder,
		std::vector<Sequence> &sequences,
		AlphabetCoder::Code sequence_delimiter) {
	concatenated_sequence_.resize(concatenated_sequences_length_);
	offsets_.resize(sequences.size() + 1);

	concatenated_sequence_[0] = sequence_delimiter;
	uint32_t offset = 1;
	for (unsigned int i = 0; i < sequences.size(); ++i) {
		offsets_[i] = offset;
		std::string sequence = sequences[i].GetSequenceData();
		coder.Encode(&sequence[0], sequence.size(),
				&concatenated_sequence_[offset]);
		offset += sequence.size();
		concatenated_sequence_[offset] = sequence_delimiter;
		++offset;
	}
	offsets_[sequences.size()] = offset;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadInfomation(std::string filename_prefix) {
	std::string filename = GetInformationFileName(filename_prefix);
	std::ifstream in;
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		in.read((char *) &number_sequences_, sizeof(number_sequences_));
		in.read((char *) &concatenated_sequences_length_,
				sizeof(concatenated_sequences_length_));

		in.close();
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadOffsets(std::string filename_prefix) {
	std::string filename = GetOffsetsFileName(filename_prefix);
	std::ifstream in;
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		offsets_.resize(number_sequences_ + 1);
		in.read((char *) &offsets_[0],
				sizeof(offsets_[0]) * (number_sequences_ + 1));
		in.close();
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadNames(std::string filename_prefix) {
	std::ifstream in;
	std::string line;
	std::string filename = GetNamesFileName(filename_prefix);
	in.open(filename.c_str());
	if (in) {
		names_.resize(number_sequences_);
		uint32_t i;
		for (i = 0; i < number_sequences_ && !in.eof(); ++i) {
			std::getline(in, line);
			names_[i] = line;
		}
		assert(i == number_sequences_);
		in.close();
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadConcatenatedSequence(
		std::string filename_prefix) {
	std::string filename = GetSequencesFileName(filename_prefix);
	std::ifstream in;
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		concatenated_sequence_.resize(concatenated_sequences_length_);
		in.read((char *) &concatenated_sequence_[0],
				sizeof(concatenated_sequence_[0])
						* concatenated_sequences_length_);
		in.close();
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadSeedSearcherCommonParameters(
		std::string filename_prefix) {
	std::string filename = GetSeedSearcherCommonParameters(filename_prefix);
	std::ifstream in;
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		seed_searcher_parameters_.BuildCommonParameters(in,
				&seed_searcher_common_parameters_);
		in.close();
		setted_seed_searcher_common_parameters_ = true;
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::LoadSeedSearcherParameters(
		std::string filename_prefix) {
	std::string filename = GetSeedSearcherParameters(filename_prefix);
	std::ifstream in;
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		seed_searcher_parameters_.Build(GetConcatenatedSequence(),
				GetConcatenatedSequenceLength(), in);
		in.close();
		setted_seed_searcher_parameters_ = true;
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveInfomation(std::string filename_prefix) {
	std::string filename = GetInformationFileName(filename_prefix);
	std::ofstream out;
	out.open(filename.c_str(), std::ios::binary);
	out.write((char *) &number_sequences_, sizeof(number_sequences_));
	out.write((char *) &concatenated_sequences_length_,
			sizeof(concatenated_sequences_length_));
	out.close();
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveOffsets(std::string filename_prefix) {
	std::ofstream out;
	std::string filename = GetOffsetsFileName(filename_prefix);
	out.open(filename.c_str(), std::ios::binary);
	out.write((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_.size());
	out.close();
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveNames(std::string filename_prefix) {
	std::ofstream out;
	std::string filename = GetNamesFileName(filename_prefix);
	out.open(filename.c_str());
	for (std::vector<std::string>::iterator i = names_.begin();
			i != names_.end(); ++i) {
		out << *i << std::endl;
	}
	out.close();
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveConcatenatedSequence(
		std::string filename_prefix) {
	std::ofstream out;
	std::string filename = GetSequencesFileName(filename_prefix);
	out.open(filename.c_str(), std::ios::binary);
	out.write((char *) &concatenated_sequence_[0],
			sizeof(concatenated_sequence_[0]) * concatenated_sequence_.size());
	out.close();
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveSeedSearcherCommonParameters(
		std::string filename_prefix) {
	std::ofstream out;
	std::string filename = GetSeedSearcherCommonParameters(filename_prefix);
	out.open(filename.c_str(), std::ios::binary);
	seed_searcher_parameters_.SaveCommonParameters(out);
	out.close();
	return true;
}

template<typename TSeedSearcher>
bool DatabaseChunk<TSeedSearcher>::SaveSeedSearcherParameters(
		std::string filename_prefix) {
	std::ofstream out;
	std::string filename = GetSeedSearcherParameters(filename_prefix);
	out.open(filename.c_str(), std::ios::binary);
	seed_searcher_parameters_.Save(out);
	out.close();
	return true;
}

#endif /* DATABASE_CHUNK_H_ */
