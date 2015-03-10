/*
 * variable_hash_clustering_seuences_index.cpp
 *
 *  Created on: 2013/05/30
 *      Author: shu
 */

#include <assert.h>
#include <algorithm>
#include "alphabet_coder.h"
#include "variable_hash_clustering_seuences_index.h"

using namespace std;

VariableHashClusteringSeuencesIndex::VariableHashClusteringSeuencesIndex() {
}

VariableHashClusteringSeuencesIndex::VariableHashClusteringSeuencesIndex(
		std::istream &is) {
	Load(is);
}

VariableHashClusteringSeuencesIndex::~VariableHashClusteringSeuencesIndex() {

}

int VariableHashClusteringSeuencesIndex::Build(
		std::vector<std::pair<HashKey, RepresentativeData> >::iterator hash_representation_pairs_first,
		std::vector<std::pair<HashKey, RepresentativeData> >::iterator hash_representation_pairs_last,
		std::vector<std::pair<ClusterId, Position> >::iterator cluster_id_member_pairs_first,
		std::vector<std::pair<ClusterId, Position> >::iterator cluster_id_pairs_last,
		const HashFunction &hash_function) {
	hash_function_ = hash_function;
	representations_index_.Build(hash_representation_pairs_first,
			hash_representation_pairs_last);
	members_index_.Build(cluster_id_member_pairs_first, cluster_id_pairs_last);
	return 0;
}

int VariableHashClusteringSeuencesIndex::GetHashKey(
		const AlphabetCoder::Code *sequence,
		VariableHashClusteringSeuencesIndex::HashKey *hash_key,
		uint32_t *hash_length) const {
	return hash_function_.CalculateHash(sequence, hash_key, hash_length);
}

int VariableHashClusteringSeuencesIndex::GetHashKey(
		const AlphabetCoder::Code *sequence,
		VariableHashClusteringSeuencesIndex::HashKey *hash_key) const {
	return hash_function_.CalculateHash(sequence, hash_key);
}

int VariableHashClusteringSeuencesIndex::Load(std::istream &is) {
	hash_function_.Load(is);
	representations_index_.Load(is);
	members_index_.Load(is);
	return 0;
}
int VariableHashClusteringSeuencesIndex::Save(std::ostream &os) {
	hash_function_.Save(os);
	representations_index_.Save(os);
	members_index_.Save(os);
	return 0;
}
