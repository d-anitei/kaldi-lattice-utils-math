// kwsbin2/utils_math.h

// Copyright (c) 2017 Joan Puigcerver <joapuipe@upv.es>
// Copyright (c) 2026 Dan Anitei <danitei at prhlt dot upv dot es>

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#ifndef KALDI_KWSBIN2_UTILS_MATH_H_
#define KALDI_KWSBIN2_UTILS_MATH_H_

#include "kwsbin2/utils.h"

namespace fst {

// Convert a CompactLattice to a FST where the cost of each arc is the
// total cost of the CompactLatticeArc (acoustic + graph costs), the
// input label is the word/char in the CompactLatticeArc, and the output label
// is an integer that maps to the initial and end frames segmentation.
template <typename Arc, typename Int>
void CompactLatticeToWordCountSegmentFst(
    const kaldi::CompactLattice& clat,
    fst::MutableFst<Arc>* ofst,
    std::vector< std::tuple<Int, Int> >* olabels) {
  using namespace fst;
  using namespace kaldi;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;

  ofst->DeleteStates();
  olabels->clear();
  if (clat.Start() == kNoStateId) return;

  std::vector<int32> state_times;
  CompactLatticeStateTimes(clat, &state_times);

  for (StateId s = 0; s < clat.NumStates(); ++s) {
    ofst->SetFinal(ofst->AddState(), ConvertToCost(clat.Final(s)));
  }
  ofst->SetStart(clat.Start());

  std::map<std::tuple<Int, Int>, Label> olabels_map;
  olabels_map[std::make_tuple(0, 0)] = 0;

  for (StateIterator<CompactLattice> sit(clat); !sit.Done(); sit.Next()) {
    const StateId s = sit.Value();
    for (ArcIterator<CompactLattice> ait(clat, s); !ait.Done(); ait.Next()) {
      const auto& arc = ait.Value();
      const auto new_weight = ConvertToCost(arc.weight);
      const auto otuple =
          std::make_tuple(state_times[s], state_times[arc.nextstate]);
      const auto olabel =
          olabels_map.emplace(otuple, olabels_map.size()).first->second;
      ofst->AddArc(s, Arc(arc.ilabel, olabel, new_weight, arc.nextstate));
    }
  }

  olabels->resize(olabels_map.size());
  for (const auto& kv : olabels_map) {
    (*olabels)[kv.second] = kv.first;
  }

  ofst->SetProperties(ofst->Properties(kFstProperties, true), kFstProperties);
}

}  // namespace fst

#endif  // KALDI_KWSBIN2_UTILS_MATH_H_
