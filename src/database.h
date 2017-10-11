/*
 * database.h
 *
 *  Created on: 2012/10/31
 *      Author: shu
 */

#ifndef DATABASE_H_
#define DATABASE_H_

#include <vector>
#include <string>
#include <stdint.h>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "sequence_type.h"
#include "fasta_sequence_reader.h"
#include "database_chunk.h"
#include "sequence.h"
#include "seed_searcher.h"
#include "seed_searcher_gpu.h"

#include <boost/thread.hpp>
#include <list>
template<typename TSeedSearcher>
class Database {
public:
	typedef int ChunkBuildOption;
	typedef TSeedSearcher SeedSearcher;
	typedef typename SeedSearcher::DatabaseParameters SeedSearcherParameters;
	typedef typename SeedSearcherParameters::BuildParameters SeedSearcherParametersBuildParameters;
	typedef typename SeedSearcherParameters::CommonParameters SeedSearcherParametersCommonParameters;

	typedef int PreloadTarget;
	typedef struct {
		unsigned int chunk_size;
		ChunkBuildOption chunk_build_option;
		std::tr1::shared_ptr<SequenceType> sequence_type_ptr;
		AlphabetCoder::Code sequence_delimiter;
		SeedSearcherParametersBuildParameters seed_search_parameters_build_parameters;
	} Parameters;

	static const ChunkBuildOption ChunkBuildOptionNotBuild = -2;
	static const ChunkBuildOption ChunkBuildOptionAllBuild = -1;

	enum PreloadTargetFlag {
		kSequence = 1 << 0,
		kName = 1 << 1,
		kOffset = 1 << 2,
		kSeedSearcherCommonParameters = 1 << 3,
		kSeedSearcherParameters = 1 << 4,
	};

	Database();
	Database(std::string filename_prefix);
	Database(std::istream &is, Parameters &parameters);
	Database(std::istream &is, Parameters &parameters,bool parallel);

	bool Load(std::string filename_prefix);
	bool Save(std::string filename_prefix);
	
	bool BuildParallel(std::string filename_prefix,int threads);
	
	
	bool SetChunk(int id);
	void ResetChunk();
	bool NextChunk();

	uint32_t GetNumberChunks() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return number_chunks_;
	}

	int GetChunkId() {
		return chunk_id_;
	}

	int GetPreloadedChunkId() {
		return preload_chunk_id_;
	}

	uint32_t GetNumberSequencesInChunk() {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetNumberSequences();
	}

	uint32_t GetConcatenatedSequenceLength() {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetConcatenatedSequenceLength();
	}

	std::string GetName(uint32_t id) {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetName(id);
	}

	uint32_t GetOffset(uint32_t id) {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetOffsets(id);
	}

	AlphabetCoder::Code *GetConcatenatedSequence() {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetConcatenatedSequence();
	}

	uint32_t GetId(uint32_t position) {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetId(position);
	}

	SeedSearcherParametersCommonParameters &GetSeedSearcherCommonParameters() {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetSeedSearcherCommonParameters();
	}

	SeedSearcherParameters &GetSeedSearcherParameters() {
		return chunk_list_[chunk_ids_[kUsedChunkId]]->GetSeedSearcherParameters();
	}

	uint64_t GetDatabaseTotalLenght() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return database_length_;
	}

	uint64_t GetNumberTotalSequences() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return number_sequences_;
	}

	uint32_t GetMaxSequenceLength() {
		if (!setting_informations_) {
			SetInfomations();
		}
		return max_sequence_length_;
	}

	AlphabetCoder::Code GetSequenceDelimiter() {
		return sequence_delimiter_;
	}

	int Preload(int chunk_id, PreloadTarget target);

private:
	typedef DatabaseChunk<SeedSearcher> DatabaseChunkType;
	typedef std::tr1::shared_ptr<DatabaseChunkType> DatabaseChunkTypePtr;
	static const size_t kNumberOfChunkIds = 2;
	enum {
		kUsedChunkId = 0, kPreloadedChunkId = 1,
	};
	static std::string GetInformationFileName(std::string filename_prefix) {
		return filename_prefix + ".inf";
	}

	void ExchangeChunkList() {
		std::swap(chunk_id_, preload_chunk_id_);
		std::swap(chunk_ids_[kUsedChunkId], chunk_ids_[kPreloadedChunkId]);
	}

	unsigned int ReadSequences(std::vector<Sequence> &sequences);
	bool SetInfomations();
	bool BuildNextChunk();
	bool LoadInfomation(std::string filename_prefix);
	bool SaveInformation(std::string filename_prefix);

	bool BuildParallelThread(int chunk_id,std::vector<Sequence> &sequences,AlphabetCoder &coder,
				 AlphabetCoder::Code sequence_delimiter,
				 SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters);

	bool all_on_memory_flag_;
	bool saving_;
	bool setting_informations_;
	Parameters parameters_;
	FastaSequenceReader reader_;
	std::string filename_prefix_;
	int number_chunks_;
	uint32_t max_sequence_length_;
	uint64_t database_length_;
	uint64_t number_sequences_;
	AlphabetCoder::Code sequence_delimiter_;
	int chunk_id_;
	int preload_chunk_id_;
	int chunk_ids_[kNumberOfChunkIds];
	std::vector<DatabaseChunkTypePtr> chunk_list_;
	Sequence next_sequence_;
};

template<typename TSeedSearcher>
Database<TSeedSearcher>::Database(std::string filename_prefix) :
		all_on_memory_flag_(false), saving_(false), setting_informations_(
				false), filename_prefix_(filename_prefix), number_chunks_(0), max_sequence_length_(
				0), database_length_(0), number_sequences_(0), sequence_delimiter_(
				0), chunk_id_(-1), preload_chunk_id_(-1), next_sequence_("", "") {
	for (size_t i = 0; i < kNumberOfChunkIds; ++i) {
		chunk_ids_[i] = i;
	}
	Load(filename_prefix);
}

template<typename TSeedSearcher>
Database<TSeedSearcher>::Database(std::istream &is, Parameters &parameters) :
		all_on_memory_flag_(false), saving_(false), setting_informations_(
				false), parameters_(parameters), reader_(is), filename_prefix_(
				""), number_chunks_(0), max_sequence_length_(0), database_length_(
				0), number_sequences_(0), sequence_delimiter_(
				parameters.sequence_delimiter), chunk_id_(-1), preload_chunk_id_(
				-1), next_sequence_("", "") {
	for (size_t i = 0; i < kNumberOfChunkIds; ++i) {
		chunk_ids_[i] = i;
	}
	BuildNextChunk();
}

template<typename TSeedSearcher>
Database<TSeedSearcher>::Database(std::istream &is, Parameters &parameters, bool parallel) :
		all_on_memory_flag_(false), saving_(false), setting_informations_(
				false), parameters_(parameters), reader_(is), filename_prefix_(
				""), number_chunks_(0), max_sequence_length_(0), database_length_(
				0), number_sequences_(0), sequence_delimiter_(
				parameters.sequence_delimiter), chunk_id_(-1), preload_chunk_id_(
				-1), next_sequence_("", "") {
	for (size_t i = 0; i < kNumberOfChunkIds; ++i) {
		chunk_ids_[i] = i;
	}
	if(!parallel){
	BuildNextChunk();
	}
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::Load(std::string filename_prefix) {
	all_on_memory_flag_ = false;
	saving_ = true;
	setting_informations_ = true;
	LoadInfomation(filename_prefix);
	SetChunk(0);
	return true;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::Save(std::string filename_prefix) {
	filename_prefix_ = filename_prefix;
	int current_id = GetChunkId();
	ResetChunk();
	do {
		std::stringstream ss;
		ss << filename_prefix_ << "_" << GetChunkId();
		chunk_list_[chunk_ids_[kUsedChunkId]]->Save(ss.str());
	} while (BuildNextChunk());
	SaveInformation(filename_prefix_);
	saving_ = true;
	setting_informations_ = true;
	SetChunk(current_id);
	return true;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::BuildParallel(std::string filename_prefix,int n_threads){
  filename_prefix_ = filename_prefix;
  int current_id = 0;
  bool build =true;
  boost::thread_group threads;
    std::list<boost::thread*> thread_list;
  
  /////////
  chunk_id_=0;
   
  while(build){
    //std::cout<<"chunk_id:"<<chunk_id_<<std::endl;
    std::vector<Sequence> sequences;
    AlphabetCoder coder(*(parameters_.sequence_type_ptr));
    chunk_list_.resize(chunk_id_ + 1);
    chunk_list_[chunk_id_] = DatabaseChunkTypePtr( new DatabaseChunkType);
    
    std::vector<AlphabetCoder::Code> encoded_sequences;
    std::vector<unsigned int> encoded_sequence_offsets;
    database_length_ += ReadSequences(sequences);
    if (sequences.empty()) {
      
      build= false;
    }
    
    if(build){
      
      for (size_t i = 0; i < sequences.size(); ++i) {
	max_sequence_length_ = std::max(max_sequence_length_,
					(uint32_t) (sequences[i].GetSequenceData().size()));
      }
      number_sequences_ += sequences.size();
      
      
      
      if (parameters_.chunk_build_option
	  == Database<TSeedSearcher>::ChunkBuildOptionAllBuild
	  || parameters_.chunk_build_option == (chunk_id_ + 1)) {
	if(thread_list.size()>=n_threads){
	  thread_list.front()->join();
	  thread_list.pop_front();
	  
	}
	
	
	
	std::cout<<"build chunk:"<<chunk_id_<<std::endl; 
	boost::thread *t = threads.create_thread(boost::bind(&Database::BuildParallelThread,this,
							     chunk_id_,sequences,coder,
							     parameters_.sequence_delimiter,
							     parameters_.seed_search_parameters_build_parameters));
	thread_list.push_back(t);
	// t->join();
	/*
	  chunk_list_[chunk_id_]->Build(sequences, coder,
	  parameters_.sequence_delimiter,
	  parameters_.seed_search_parameters_build_parameters);
	  
	  
	  
	  std::cout<<"save chunk:"<<chunk_id_<<std::endl;
	  std::stringstream ss;
	  ss << filename_prefix_ << "_" << GetChunkId();
	  std::cout<<ss.str()<<std::endl;
	  chunk_list_[chunk_id_]->Save(ss.str());
	  
	*/
	
      } else {
	std::cout<<"skip chunk:"<<chunk_id_<<std::endl;
	chunk_list_[chunk_id_]->Clear();
      }
      
      ++chunk_id_;
      build=true;
    }
  }
  threads.join_all();
  
  number_chunks_=chunk_id_;
  SaveInformation(filename_prefix_);
  saving_ = true;
  setting_informations_ = true;
  SetChunk(current_id);
  return true;
  
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::BuildParallelThread(int chunk_id,
						  std::vector<Sequence> &sequences,AlphabetCoder &coder,
						  AlphabetCoder::Code sequence_delimiter,
						  SeedSearcherParametersBuildParameters &seed_search_parameters_build_parameters)
{
  DatabaseChunkTypePtr chunk = DatabaseChunkTypePtr( new DatabaseChunkType);
  chunk->Build(sequences, coder,
	       parameters_.sequence_delimiter,
	       parameters_.seed_search_parameters_build_parameters);
  
    std::stringstream ss;
  ss << filename_prefix_ << "_" << chunk_id;
  chunk->Save(ss.str());
  chunk->Clear();
  sequences.clear();
   return true;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::SetChunk(int id) {
	if (id == chunk_id_) {
		return true;
	}
	if (!saving_) {
		reader_.Seek(0);
		for (int i = 0; i < id; ++i) {
			BuildNextChunk();
		}
	} else {
		if (id >= number_chunks_) {
			return false;
		}
		if (preload_chunk_id_ == id) {
			//std::cout << "set chunk id " << id << std::endl;
			//std::cout << "preload chunk id " << preload_chunk_id_ << std::endl;
			ExchangeChunkList();
			//std::cout << "exchange preload chunk" << std::endl;
		} else {
			std::stringstream prefix;
			prefix.str("");
			prefix << filename_prefix_ << "_" << id;
			if (all_on_memory_flag_) {
				chunk_ids_[kUsedChunkId] = id;
				if (!chunk_list_[chunk_ids_[kUsedChunkId]]) {
					chunk_list_[chunk_ids_[kUsedChunkId]] =
							DatabaseChunkTypePtr(new DatabaseChunkType);
					chunk_list_[chunk_ids_[kUsedChunkId]]->Load(prefix.str());
				}
			} else {
				if (!chunk_list_[chunk_ids_[kUsedChunkId]]) {
					chunk_list_[chunk_ids_[kUsedChunkId]] =
							DatabaseChunkTypePtr(new DatabaseChunkType);
				}
				chunk_list_[chunk_ids_[kUsedChunkId]]->Load(prefix.str());
			}
			chunk_id_ = id;
		}
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
void Database<TSeedSearcher>::ResetChunk() {
	SetChunk(0);
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::NextChunk() {
	if (!saving_) {
		return BuildNextChunk();
	} else if (chunk_id_ + 1 < number_chunks_) {
		SetChunk(chunk_id_ + 1);
		return true;
	} else {
		chunk_id_ = number_chunks_;
		return false;
	}
}

template<typename TSeedSearcher>
unsigned int Database<TSeedSearcher>::ReadSequences(
		std::vector<Sequence> &sequences) {
	unsigned int sum_length = 0;
	std::string name;
	std::string sequence;
	sequences.clear();
	bool reader_ret = true;
	if (next_sequence_.GetName().size() == 0
			&& next_sequence_.GetSequenceData().size() == 0) {
		reader_ret = reader_.Read(name, sequence);
		next_sequence_ = Sequence(name, sequence);
	}

	while (reader_ret) {
		unsigned int new_sum_length = sum_length
				+ next_sequence_.GetSequenceData().size();
		if (new_sum_length >= parameters_.chunk_size) {
			break;
		}
		sequences.push_back(next_sequence_);
		sum_length = new_sum_length;
		reader_ret = reader_.Read(name, sequence);
		if (reader_ret) {
			next_sequence_ = Sequence(name, sequence);
		} else {
			next_sequence_ = Sequence("", "");
			break;
		}
	}
	return sum_length;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::SetInfomations() {
	if (setting_informations_) {
		return true;
	} else if (reader_.IsRead()) {
		uint32_t id = GetChunkId();
		while (BuildNextChunk()) {
		}
		setting_informations_ = true;
		SetChunk(id);
		return true;
	} else {
		return false;
	}
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::BuildNextChunk() {
	std::vector<Sequence> sequences;
	AlphabetCoder coder(*(parameters_.sequence_type_ptr));
	if (all_on_memory_flag_) {
		chunk_ids_[kUsedChunkId] = chunk_id_;
		chunk_list_.resize(chunk_id_ + 1);
	} else if (chunk_list_.empty()) {
		chunk_list_.resize(1);
		assert(chunk_ids_[kUsedChunkId] == 0);
	}
	if (!chunk_ids_[kUsedChunkId]) {
		chunk_list_[chunk_ids_[kUsedChunkId]] = DatabaseChunkTypePtr(
				new DatabaseChunkType);
	}

	std::vector<AlphabetCoder::Code> encoded_sequences;
	std::vector<unsigned int> encoded_sequence_offsets;
	database_length_ += ReadSequences(sequences);
	if (sequences.empty()) {
		number_chunks_ = chunk_id_ + 1;
		return false;
	}
	for (size_t i = 0; i < sequences.size(); ++i) {
		max_sequence_length_ = std::max(max_sequence_length_,
				(uint32_t) (sequences[i].GetSequenceData().size()));
	}
	number_sequences_ += sequences.size();
	if (parameters_.chunk_build_option
			== Database<TSeedSearcher>::ChunkBuildOptionAllBuild
			|| parameters_.chunk_build_option == (chunk_id_ + 1)) {
		chunk_list_[chunk_ids_[kUsedChunkId]]->Build(sequences, coder,
				parameters_.sequence_delimiter,
				parameters_.seed_search_parameters_build_parameters);
	} else {
		chunk_list_[chunk_ids_[kUsedChunkId]]->Clear();
	}
	++chunk_id_;
	return true;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::LoadInfomation(std::string filename_prefix) {
	filename_prefix_ = filename_prefix;
	std::ifstream in;
	std::string filename = GetInformationFileName(filename_prefix_);
	in.open(filename.c_str(), std::ios::binary);
	if (in) {
		in.read((char *) &number_chunks_, sizeof(number_chunks_));
		in.read((char *) &max_sequence_length_, sizeof(max_sequence_length_));
		in.read((char *) &database_length_, sizeof(database_length_));
		in.read((char *) &number_sequences_, sizeof(number_sequences_));
		in.read((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
		in.close();
		if (all_on_memory_flag_) {
			chunk_list_.resize(number_chunks_);
		} else {
			chunk_list_.resize(kNumberOfChunkIds);
		}
		return true;
	}
	return false;
}

template<typename TSeedSearcher>
bool Database<TSeedSearcher>::SaveInformation(std::string filename_prefix) {
	filename_prefix_ = filename_prefix;
	std::string filename = GetInformationFileName(filename_prefix);
	std::ofstream out;
	out.open(filename.c_str(), std::ios::binary);
	out.write((char *) &number_chunks_, sizeof(number_chunks_));
	out.write((char *) &max_sequence_length_, sizeof(max_sequence_length_));
	out.write((char *) &database_length_, sizeof(database_length_));
	out.write((char *) &number_sequences_, sizeof(number_sequences_));
	out.write((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
	out.close();

	return true;
}

template<typename TSeedSearcher>
int Database<TSeedSearcher>::Preload(int id, PreloadTarget target) {
	//std::cout << "preload chunk " << id << std::endl;
	if (!saving_) {
		return 1;
	}
	if (id >= number_chunks_) {
		return 1;
	}
	size_t chunk_list_id = 0;
	if (all_on_memory_flag_) {
		chunk_list_id = id;
		chunk_ids_[kPreloadedChunkId] = chunk_list_id;
		if (!chunk_list_[chunk_list_id]) {
			chunk_list_[chunk_list_id] = DatabaseChunkTypePtr(
					new DatabaseChunkType);
			std::stringstream prefix;
			prefix.str("");
			prefix << filename_prefix_ << "_" << id;
			chunk_list_[chunk_list_id]->Load(prefix.str());
		}
		preload_chunk_id_ = id;
	} else if (id == chunk_id_) {
		chunk_list_id = chunk_ids_[kUsedChunkId];
	} else {
		chunk_list_id = chunk_ids_[kPreloadedChunkId];
		if (!chunk_list_[chunk_list_id]) {
			chunk_list_[chunk_list_id] = DatabaseChunkTypePtr(
					new DatabaseChunkType);
		}
		std::stringstream prefix;
		prefix.str("");
		prefix << filename_prefix_ << "_" << id;
		chunk_list_[chunk_list_id]->Load(prefix.str());
		preload_chunk_id_ = id;
	}
	if (target & kSequence) {
		chunk_list_[chunk_list_id]->GetConcatenatedSequence();
	}
	if (target & kName) {
		chunk_list_[chunk_list_id]->GetName(0);
	}
	if (target & kOffset) {
		chunk_list_[chunk_list_id]->GetOffsets(0);
	}
	if (target & kSeedSearcherCommonParameters) {
		chunk_list_[chunk_list_id]->GetSeedSearcherCommonParameters();
	}
	if (target & kSeedSearcherParameters) {
		chunk_list_[chunk_list_id]->GetSeedSearcherParameters();
	}
	//std::cout << "finished preloading chunk " << id << std::endl;
	return 0;
}

#endif /* DATABASE_H_ */
