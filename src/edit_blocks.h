/*
 * edit_array.h
 *
 *  Created on: 2010/11/28
 *      Author: shu
 */

#ifndef EDIT_BLOCKS_H_
#define EDIT_BLOCKS_H_

#include <vector>

class EditBlocks {
public:
  enum EditOpType{
    kGapInSeq0,
    kSubstitution,
    kGapInSeq1
  };

  typedef struct {
    EditOpType op;
    int length;
  }Block;

  void Clear();
  void Add(EditBlocks::EditOpType op_type , int length);
  void Add(EditBlocks other);
  void Reverse();
  std::vector<EditBlocks::EditOpType> ToVector();

private:
  std::vector<Block> blocks_;
};

#endif /* EDIT_ARRAY_H_ */
