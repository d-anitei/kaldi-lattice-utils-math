// fstext/fstext-utils2_math.h

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

#ifndef KALDI_FSTEXT_FSTEXT_UTILS2_MATH_H_
#define KALDI_FSTEXT_FSTEXT_UTILS2_MATH_H_

#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "fstext/fstext-utils2.h"

namespace fst {

template <typename Arc>
void SegmentWFST(
    std::unordered_map<typename Arc::StateId,
                       std::set<typename Arc::StateId>>& segments,
    std::unordered_map<typename Arc::StateId,
                       std::set<typename Arc::StateId>>& segments_bin,
    std::unordered_map<typename Arc::StateId,
                       std::set<std::tuple<typename Arc::StateId,
                                           typename Arc::Label,
                                           typename Arc::StateId>>>& scope_closes,
    MutableFst<Arc>& fst,
    int32& final_state,
    std::set<typename Arc::StateId>& final_states,
    const std::unordered_set<typename Arc::Label>& open_scope_labels,
    const std::unordered_map<typename Arc::Label, typename Arc::Label>&
        close_to_open_scope_label,
    const std::unordered_set<typename Arc::Label>&
        operator_segmentation_labels) {
  typedef MutableFst<Arc> Fst;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  if (fst.Start() == kNoStateId) return;

  const StateId super_final = fst.AddState();
  std::unordered_map<StateId, std::set<std::tuple<StateId, Label, StateId>>>
      segmentation_points;

  for (StateIterator<Fst> sit(fst); !sit.Done(); sit.Next()) {
    const StateId state = sit.Value();

    if (fst.Final(state) != Weight::Zero()) {
      final_states.insert(state);
    }
    if (state == super_final) continue;

    if (state == fst.Start()) {
      segmentation_points[state] = {std::make_tuple(state, 0, state)};
    }

    for (ArcIterator<Fst> ait(fst, state); !ait.Done(); ait.Next()) {
      const Arc& arc = ait.Value();
      const StateId next_state = arc.nextstate;
      const auto label = arc.ilabel;
      const auto current_points = segmentation_points[state];

      if (open_scope_labels.count(label) != 0) {
        segmentation_points[next_state].insert(
            std::make_tuple(state, label, next_state));
      } else if (close_to_open_scope_label.count(label) != 0) {
        const auto expected_open_label =
            close_to_open_scope_label.find(label)->second;
        for (const auto& point : current_points) {
          const auto open_label = std::get<1>(point);
          if (open_label != expected_open_label) continue;

          const auto previous = segmentation_points[std::get<0>(point)];
          segmentation_points[next_state].insert(previous.begin(),
                                                 previous.end());
          const auto scope_start = std::get<2>(point);
          segments[scope_start].insert(next_state);
          scope_closes[scope_start].emplace(state, label, next_state);
        }
      } else if (operator_segmentation_labels.count(label) != 0) {
        for (const auto& point : current_points) {
          segments_bin[std::get<2>(point)].insert(next_state);
          segmentation_points[next_state].insert(point);
        }
      } else {
        segmentation_points[next_state].insert(current_points.begin(),
                                               current_points.end());
      }
    }
  }

  fst.SetFinal(super_final, Weight::One());
  final_state = super_final;
}

}  // namespace fst

#endif  // KALDI_FSTEXT_FSTEXT_UTILS2_MATH_H_
