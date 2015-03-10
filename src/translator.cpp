/*
 * translator.cpp
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#include "translator.h"
#include "dna_sequence.h"
#include "protein_sequence.h"
#include "protein_type.h"
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <stdint.h>

using namespace std;

Translator::Translator() {
  start_codons_.push_back("ATG");

  codon_table_.insert(map<string, string>::value_type("AAA", "K"));
  codon_table_.insert(map<string, string>::value_type("AAC", "N"));
  codon_table_.insert(map<string, string>::value_type("AAG", "K"));
  codon_table_.insert(map<string, string>::value_type("AAT", "N"));

  codon_table_.insert(map<string, string>::value_type("ACA", "T"));
  codon_table_.insert(map<string, string>::value_type("ACC", "T"));
  codon_table_.insert(map<string, string>::value_type("ACG", "T"));
  codon_table_.insert(map<string, string>::value_type("ACT", "T"));

  codon_table_.insert(map<string, string>::value_type("AGA", "R"));
  codon_table_.insert(map<string, string>::value_type("AGC", "S"));
  codon_table_.insert(map<string, string>::value_type("AGG", "R"));
  codon_table_.insert(map<string, string>::value_type("AGT", "S"));

  codon_table_.insert(map<string, string>::value_type("ATA", "I"));
  codon_table_.insert(map<string, string>::value_type("ATC", "I"));
  codon_table_.insert(map<string, string>::value_type("ATG", "M"));
  codon_table_.insert(map<string, string>::value_type("ATT", "I"));

  codon_table_.insert(map<string, string>::value_type("CAA", "Q"));
  codon_table_.insert(map<string, string>::value_type("CAC", "H"));
  codon_table_.insert(map<string, string>::value_type("CAG", "Q"));
  codon_table_.insert(map<string, string>::value_type("CAT", "H"));

  codon_table_.insert(map<string, string>::value_type("CCA", "P"));
  codon_table_.insert(map<string, string>::value_type("CCC", "P"));
  codon_table_.insert(map<string, string>::value_type("CCG", "P"));
  codon_table_.insert(map<string, string>::value_type("CCT", "P"));

  codon_table_.insert(map<string, string>::value_type("CGA", "R"));
  codon_table_.insert(map<string, string>::value_type("CGC", "R"));
  codon_table_.insert(map<string, string>::value_type("CGG", "R"));
  codon_table_.insert(map<string, string>::value_type("CGT", "R"));

  codon_table_.insert(map<string, string>::value_type("CTA", "L"));
  codon_table_.insert(map<string, string>::value_type("CTC", "L"));
  codon_table_.insert(map<string, string>::value_type("CTG", "L"));
  codon_table_.insert(map<string, string>::value_type("CTT", "L"));

  codon_table_.insert(map<string, string>::value_type("GAA", "E"));
  codon_table_.insert(map<string, string>::value_type("GAC", "D"));
  codon_table_.insert(map<string, string>::value_type("GAG", "E"));
  codon_table_.insert(map<string, string>::value_type("GAT", "D"));

  codon_table_.insert(map<string, string>::value_type("GCA", "A"));
  codon_table_.insert(map<string, string>::value_type("GCC", "A"));
  codon_table_.insert(map<string, string>::value_type("GCG", "A"));
  codon_table_.insert(map<string, string>::value_type("GCT", "A"));

  codon_table_.insert(map<string, string>::value_type("GGA", "G"));
  codon_table_.insert(map<string, string>::value_type("GGC", "G"));
  codon_table_.insert(map<string, string>::value_type("GGG", "G"));
  codon_table_.insert(map<string, string>::value_type("GGT", "G"));

  codon_table_.insert(map<string, string>::value_type("GTA", "V"));
  codon_table_.insert(map<string, string>::value_type("GTC", "V"));
  codon_table_.insert(map<string, string>::value_type("GTG", "V"));
  codon_table_.insert(map<string, string>::value_type("GTT", "V"));

  codon_table_.insert(map<string, string>::value_type("TAA", "*"));
  codon_table_.insert(map<string, string>::value_type("TAC", "Y"));
  codon_table_.insert(map<string, string>::value_type("TAG", "*"));
  codon_table_.insert(map<string, string>::value_type("TAT", "Y"));

  codon_table_.insert(map<string, string>::value_type("TCA", "S"));
  codon_table_.insert(map<string, string>::value_type("TCC", "S"));
  codon_table_.insert(map<string, string>::value_type("TCG", "S"));
  codon_table_.insert(map<string, string>::value_type("TCT", "S"));

  codon_table_.insert(map<string, string>::value_type("TGA", "*"));
  codon_table_.insert(map<string, string>::value_type("TGC", "C"));
  codon_table_.insert(map<string, string>::value_type("TGG", "W"));
  codon_table_.insert(map<string, string>::value_type("TGT", "C"));

  codon_table_.insert(map<string, string>::value_type("TTA", "L"));
  codon_table_.insert(map<string, string>::value_type("TTC", "F"));
  codon_table_.insert(map<string, string>::value_type("TTG", "L"));
  codon_table_.insert(map<string, string>::value_type("TTT", "F"));
}

void Translator::Translate(const DnaSequence &dna, vector<ProteinSequence> &proteins) {
  string sequences[2];
  sequences[0] = dna.GetSequenceData();
  sequences[1] = dna.GetComplementaryStrand();
  stringstream new_seq;
  char unknown_letter = ProteinType().GetUnknownLetter();
  uint32_t dna_length = sequences[0].length();
  bool stoped = false;

  proteins.clear();

  if (dna_length < 3) {
    return;
  }

  for (int i = 0; i < 2; ++i) {
    string seq = sequences[i];
    for (int offset = 0; offset < 3; ++offset) {
      new_seq.str("");
      bool translated = true;
      for (uint32_t j = offset; j < dna_length - 2; j += 3) {
        string codon = seq.substr(j, 3);
        if (find(start_codons_.begin(), start_codons_.end(), codon) != start_codons_.end()) {
          translated = true;
        }
        if ((translated || !stoped) && codon_table_.find(codon) != codon_table_.end()) {
          new_seq << codon_table_[codon];
          if (codon_table_[codon] == "*") {
            translated = false;
          }
        } else {
          new_seq << unknown_letter;
        }
      }
      string protein_sequece = new_seq.str();
      proteins.push_back(ProteinSequence(dna.GetName(), protein_sequece));
    }
  }
}
