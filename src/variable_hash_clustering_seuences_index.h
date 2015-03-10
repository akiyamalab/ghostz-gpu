/*
 * variable_hash_clustering_seuences_index.h
 *
 *  Created on: 2013/05/30
 *      Author: shu
 */

#ifndef VARIABLE_HASH_CLUSTERING_SEUENCES_INDEX_H_
#define VARIABLE_HASH_CLUSTERING_SEUENCES_INDEX_H_

#include "seed_searcher_common.h"
#include "dense_hash_inverted_index.h"
#include "perfect_hash_inverted_index.h"

class VariableHashClusteringSeuencesIndex {
public:
	typedef SeedSearcherCommon::HashFunction HashFunction;
	typedef HashFunction::Hash HashKey;
	typedef uint32_t ClusterId;
	typedef uint32_t Position;
	typedef std::pair<ClusterId, Position> RepresentativeData;

	VariableHashClusteringSeuencesIndex();
	VariableHashClusteringSeuencesIndex(std::istream &is);
	virtual ~VariableHashClusteringSeuencesIndex();

	int Build(
			std::vector<std::pair<HashKey, RepresentativeData> >::iterator hash_representation_pairs_first,
			std::vector<std::pair<HashKey, RepresentativeData> >::iterator hash_representation_pairs_last,
			std::vector<std::pair<ClusterId, Position> >::iterator cluster_id_member_pairs_first,
			std::vector<std::pair<ClusterId, Position> >::iterator cluster_id_pairs_last,
			const HashFunction &hash_function);

	int GetHashKey(const AlphabetCoder::Code *sequence,
			HashKey *hash_key, uint32_t *hash_length) const;

	int GetHashKey(const AlphabetCoder::Code *sequence,
			HashKey *hash_key) const;
	int GetRepresentatives(const HashKey &hash_key,
			RepresentativeData const **representations,
			size_t *representations_length) const;
	int GetMembers(const ClusterId &cluster_id, Position const **members,
			size_t *members_length) const;
	int Load(std::istream &is);
	int Save(std::ostream &os);
private:
	//typedef InvertedIndex<HashFunction::Hash, Value> Index;
	typedef DenseHashInvertedIndex<HashFunction::Hash, RepresentativeData> RepresentationsIndex;
	typedef PerfectHashInvertedIndex<ClusterId, Position> MembersIndex;
	HashFunction hash_function_;
	RepresentationsIndex representations_index_;
	MembersIndex members_index_;
};

inline int VariableHashClusteringSeuencesIndex::GetRepresentatives(
		const HashKey &hash_key, RepresentativeData const **representations,
		size_t *representations_length) const {
	return representations_index_.GetValues(hash_key, representations,
			representations_length);
}

inline int VariableHashClusteringSeuencesIndex::GetMembers(
		const ClusterId &cluster_id, Position const **members,
		size_t *members_length) const {
	return members_index_.GetValues(cluster_id, members,
			members_length);
}

#endif /* VARIABLE_HASH_CLUSTERING_SEUENCES_INDEX_H_ */
