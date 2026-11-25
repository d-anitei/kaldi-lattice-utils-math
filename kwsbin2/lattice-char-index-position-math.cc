// kwsbin2/lattice-char-index-position-math.cc

// Copyright (c) 2017-2018 Joan Puigcerver <joapuipe@upv.es>
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

#include "base/kaldi-common.h"
#include "base/timer.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils2_math.h"
#include "fstext/fst-info.h"
#include "kwsbin2/utils_math.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

namespace kaldi {

struct AggregateTimingStats {
  std::mutex mutex;
  double overall_elapsed = 0.0;
  double indexing_elapsed = 0.0;
  double segment_wfst_elapsed = 0.0;
  double extract_region_elapsed = 0.0;
  double compute_prix_elapsed = 0.0;
  int64 segment_wfst_calls = 0;
  int64 extract_region_calls = 0;
  int64 compute_prix_calls = 0;
  int64 indexed_expressions = 0;
  int64 lattices = 0;
  int64 discarded_lattices = 0;
};

AggregateTimingStats& GetAggregateTimingStats() {
  static AggregateTimingStats stats;
  return stats;
}

void AddAggregateTimingStats(
    double overall_elapsed,
    double indexing_elapsed,
    double segment_wfst_elapsed,
    double extract_region_elapsed,
    double compute_prix_elapsed,
    int64 segment_wfst_calls,
    int64 extract_region_calls,
    int64 compute_prix_calls,
    int32 indexed_expressions) {
  auto& stats = GetAggregateTimingStats();
  std::lock_guard<std::mutex> lock(stats.mutex);
  stats.overall_elapsed += overall_elapsed;
  stats.indexing_elapsed += indexing_elapsed;
  stats.segment_wfst_elapsed += segment_wfst_elapsed;
  stats.extract_region_elapsed += extract_region_elapsed;
  stats.compute_prix_elapsed += compute_prix_elapsed;
  stats.segment_wfst_calls += segment_wfst_calls;
  stats.extract_region_calls += extract_region_calls;
  stats.compute_prix_calls += compute_prix_calls;
  stats.indexed_expressions += indexed_expressions;
  stats.lattices += 1;
}

void AddDiscardedLattice() {
  auto& stats = GetAggregateTimingStats();
  std::lock_guard<std::mutex> lock(stats.mutex);
  stats.discarded_lattices += 1;
}

void LogAggregateTimingStats() {
  auto& stats = GetAggregateTimingStats();
  std::lock_guard<std::mutex> lock(stats.mutex);

  KALDI_LOG << "Aggregate timing summary: lattices=" << stats.lattices
            << ", discarded_lattices=" << stats.discarded_lattices
            << ", indexed_expressions=" << stats.indexed_expressions
            << ", overall=" << stats.overall_elapsed
            << "s, indexing=" << stats.indexing_elapsed
            << "s, SegmentWFST=" << stats.segment_wfst_elapsed
            << "s (" << stats.segment_wfst_calls << " calls)"
            << ", ExtractRegionFst=" << stats.extract_region_elapsed
            << "s (" << stats.extract_region_calls << " calls)"
            << ", ComputePrIxFST=" << stats.compute_prix_elapsed
            << "s (" << stats.compute_prix_calls << " calls)";

  if (stats.indexed_expressions > 0) {
    KALDI_LOG << "Aggregate per-expression estimates: overall="
              << (stats.overall_elapsed / stats.indexed_expressions)
              << "s, indexing="
              << (stats.indexing_elapsed / stats.indexed_expressions)
              << "s, SegmentWFST="
              << (stats.segment_wfst_elapsed / stats.indexed_expressions)
              << "s, ExtractRegionFst="
              << (stats.extract_region_elapsed / stats.indexed_expressions)
              << "s, ComputePrIxFST="
              << (stats.compute_prix_elapsed / stats.indexed_expressions)
              << "s.";
  }
}

void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam) {
  Timer timer;
  const int64 narcs = NumArcs(*clat);
  const int64 nstates = clat->NumStates();
  // Acoustic scale
  if (acoustic_scale != 1.0 || graph_scale != 1.0) {
    ScaleLattice(fst::LatticeScale(graph_scale, acoustic_scale), clat);
  }
  // Word insertion penalty
  if (insertion_penalty != 0.0) {
    AddWordInsPenToCompactLattice(insertion_penalty, clat);
  }
  // Lattice pruning
  if (beam != std::numeric_limits<BaseFloat>::infinity()) {
    PruneLattice(beam, clat);
  }
  const int64 pruned_narcs = NumArcs(*clat);
  const int64 pruned_nstates = clat->NumStates();
  KALDI_VLOG(1) << "Lattice " << key << ": pruned #states from "
                << nstates << " to " << pruned_nstates << " and #arcs from "
                << narcs << " to " << pruned_narcs;
  if (clat->Start() != fst::kNoStateId) {
    // If needed, sort the compact lattice in topological order
    TopSortCompactLatticeIfNeeded(clat);
  }
  KALDI_VLOG(1) << "Lattice " << key << ": Preprocessing done in "
                << timer.Elapsed() << " seconds.";
}

bool StartsWith(const std::string& str, const std::string& prefix) {
  return str.compare(0, prefix.size(), prefix) == 0;
}

template <typename Label>
Label ResolveSymbolLabel(const fst::SymbolTable& symbols,
                         const std::string& symbol) {
  const auto label = symbols.Find(symbol);
  if (label == fst::kNoLabel) {
    KALDI_ERR << "Symbol table does not contain required symbol: " << symbol;
  }
  return static_cast<Label>(label);
}

template <typename Label>
std::vector<Label> ResolveSymbolVector(
    const fst::SymbolTable& symbols,
    const std::string& symbols_str) {
  std::vector<std::string> symbol_names;
  kaldi::SplitStringToVector(symbols_str, " ", true, &symbol_names);

  std::vector<Label> labels;
  labels.reserve(symbol_names.size());
  for (const auto& symbol : symbol_names) {
    labels.push_back(ResolveSymbolLabel<Label>(symbols, symbol));
  }
  return labels;
}

template <typename Label>
std::unordered_set<Label> ToLabelSet(const std::vector<Label>& label_vector) {
  std::unordered_set<Label> labels;
  labels.insert(label_vector.begin(), label_vector.end());
  return labels;
}

template <typename Label>
std::unordered_set<Label> ResolvePrefixedSymbols(
    const fst::SymbolTable& symbols,
    const std::string& prefixes_str) {
  std::vector<std::string> prefixes;
  kaldi::SplitStringToVector(prefixes_str, " ", true, &prefixes);
  if (prefixes.empty()) {
    KALDI_ERR << "At least one operator segmentation prefix is required.";
  }

  std::unordered_set<Label> labels;
  for (fst::SymbolTableIterator sit(symbols); !sit.Done(); sit.Next()) {
    const auto symbol = sit.Symbol();
    for (const auto& prefix : prefixes) {
      if (StartsWith(symbol, prefix)) {
        labels.insert(static_cast<Label>(sit.Value()));
        break;
      }
    }
  }

  if (labels.empty()) {
    KALDI_ERR << "No symbols matched --operator-segment-prefixes="
              << prefixes_str;
  }
  return labels;
}

template <typename Label>
std::unordered_map<Label, Label> BuildCloseToOpenDelimiterMap(
    const std::vector<Label>& open_labels,
    const std::vector<Label>& close_labels) {
  if (open_labels.size() != close_labels.size()) {
    KALDI_ERR << "--open-delimiters and --close-delimiters must have the same "
              << "number of symbols.";
  }

  std::unordered_map<Label, Label> close_to_open;
  for (size_t i = 0; i < open_labels.size(); ++i) {
    close_to_open[close_labels[i]] = open_labels[i];
  }
  return close_to_open;
}

template <typename Arc>
void ComputePrIxFST(
  std::vector<std::tuple<std::string, int32, int32, int32, double>>& nbest,
  fst::VectorFst<Arc>& fst_subexp,
  std::vector< std::tuple<int32, int32> >& label_to_segm,
  const std::string& key,
  const double total_cost,
  double max_neg_log_posterior,
  int32 pos,
  BaseFloat determinize_delta_,
  int32 num_nbest_
){

  // Remove epsilon symbols O(V^2 + V * E)
  fst::RmEpsilon(&fst_subexp);
  // Remove unnecessary states/arcs
  fst::Connect(&fst_subexp);
  // Push weights to the toward the initial state. This speeds up n-best list
  // retrieval.
  fst::Push<fst::LogArc>(&fst_subexp, fst::REWEIGHT_TO_INITIAL, fst::kPushWeights);    

  fst::SymbolToPathSegmentationFst(&fst_subexp, label_to_segm);

  // Determinize in the log semiring to sum weights of all paths of the
  // same (word, position).
  fst::VectorFst<fst::LogArc> det_log;
  fst::Determinize(
      fst::ProjectFst<fst::LogArc>(fst_subexp, fst::PROJECT_INPUT),
      &det_log,
      fst::DeterminizeOptions<fst::LogArc>(
          determinize_delta_,
          fst::LogArc::Weight::Zero(),
          fst::kNoStateId,
          0,
          fst::DETERMINIZE_FUNCTIONAL));
  // KALDI_VLOG(1) << "Lattice " << key << ": "
  //               << fst::ComputeNumberOfPaths(&det_log) << " pseudo-words.";

  // Determinize in the tropical semiring to obtain the best segmentation
  // for each (word, position).
  using AM1 = fst::WeightConvertMapper<fst::LogArc, fst::StdArc>;
  fst::VectorFst<fst::StdArc> det_std;
  fst::Determinize(
      fst::ArcMapFst<fst::LogArc, fst::StdArc, AM1>(fst_subexp, AM1()),
      &det_std,
      fst::DeterminizeOptions<fst::StdArc>(
          determinize_delta_,
          fst::StdArc::Weight::Zero(),
          fst::kNoStateId,
          0,
          fst::DETERMINIZE_DISAMBIGUATE));

  fst_subexp.DeleteStates();   // Free space.

  // Remove weights from the tropical fst
  fst::ArcMap(&det_std, fst::RmWeightMapper<fst::StdArc, fst::StdArc>());

  // Compose the two deterministic fsts to get a new deterministic fst with
  // the best segmentation and the total cost for each pseudo-word.
  fst::VectorFst<fst::StdArc> det_subpaths;
  fst::Compose(fst::ArcMapFst<fst::LogArc, fst::StdArc, AM1>(det_log, AM1()),
               det_std, &det_subpaths);
  det_log.DeleteStates();   // Free space.
  det_std.DeleteStates();   // Free space.

  // Obtain n-best word segmentations
  fst::VectorFst<fst::StdArc> nbest_fst;
  fst::ShortestPath(det_subpaths, &nbest_fst, num_nbest_);

  // Process n-best list to prepare the output
  std::vector<fst::VectorFst<fst::StdArc>> nbest_fsts;
  ConvertNbestToVector(nbest_fst, &nbest_fsts);
  for (const auto& nbfst : nbest_fsts) {
    std::vector<fst::StdArc::Label> nb_isyms, nb_osyms;
    fst::StdArc::Weight nb_cost;
    GetLinearSymbolSequence(nbfst, &nb_isyms, &nb_osyms, &nb_cost);

    const auto pseudoword =
        LabelSequenceToString(nb_isyms, /* skip_epsilon= */ true);

    if (pseudoword.empty() || nb_osyms.empty()) {
      continue;
    }
    const auto t0 = nb_osyms.front() - 1;
    const auto t1 = nb_osyms.back() - 1;
    const double log_posterior = total_cost - nb_cost.Value();
    if (max_neg_log_posterior != std::numeric_limits<double>::infinity() &&
        !(log_posterior > -max_neg_log_posterior)) {
      continue;
    }
    nbest.emplace_back(pseudoword, pos, t0, t1, log_posterior);
  }
}

using ScopeCloseArcs = std::set<std::tuple<int32, int32, int32>>;

template <typename Arc, typename W>
bool ExtractRegionFst(
    const fst::Fst<Arc>& fst,
    typename Arc::StateId start_state,
    const std::set<typename Arc::StateId>& end_states,
    const std::vector<W>& fw, const std::vector<W>& bw,
    int is_binary_type,
    const ScopeCloseArcs* scope_closes,
    const std::set<typename Arc::StateId>& final_states,
    const std::unordered_set<typename Arc::Label>& operator_segmentation_labels,
    fst::VectorFst<Arc>* region_fst) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  KALDI_ASSERT(region_fst != nullptr);
  region_fst->DeleteStates();
  if (end_states.empty()) return false;

  const auto max_end_state = *prev(end_states.end());
  typedef std::tuple<StateId, int> StateTuple;
  std::map<StateId, StateId> state_map;
  std::set<StateId> final_links_added;
  std::map<StateTuple, bool> visited;

  const StateId output_start = region_fst->AddState();
  const StateId output_final = region_fst->AddState();
  state_map[start_state] = output_start;
  region_fst->SetStart(output_start);
  region_fst->SetFinal(output_final, Weight::One());

  std::queue<StateId> Q;
  Q.push(start_state);

  auto map_state = [&](StateId orig_state) -> StateId {
    auto it = state_map.find(orig_state);
    if (it != state_map.end()) return it->second;
    const StateId mapped_state = region_fst->AddState();
    state_map[orig_state] = mapped_state;
    return mapped_state;
  };

  while (!Q.empty()) {
    const auto n = Q.front();
    Q.pop();

    for (fst::ArcIterator<fst::Fst<Arc>> ait(fst, n); !ait.Done(); ait.Next()) {
      const Arc& arc = ait.Value();
      const StateId v = arc.nextstate;

      if (v > max_end_state) continue;
      const auto t = std::make_tuple(n, arc.ilabel);
      if (!visited.emplace(t, true).second) continue;
      const StateId output_src = map_state(n);
      auto weight = arc.weight;
      if (n == start_state) {
        weight = fst::Times(fw[n], weight);
      }

      if (end_states.count(v) != 0 && scope_closes != nullptr &&
          scope_closes->count(std::make_tuple(n, arc.ilabel, v)) != 0) {
        // Consume the matching delimiter and its weight without emitting it.
        region_fst->AddArc(output_src, Arc(0, arc.olabel,
                                         fst::Times(weight, bw[v]),
                                         output_final));
        continue;
      }

      if (end_states.count(v) != 0 && scope_closes == nullptr) {
        const auto close_weight = fst::Times(weight, bw[v]);
        if (n == start_state) {
          if (operator_segmentation_labels.count(arc.ilabel) == 0) {
            region_fst->AddArc(
                output_src,
                Arc(arc.ilabel, arc.olabel, close_weight, output_final));
          }
          continue;
        } else if (!is_binary_type) {
          region_fst->AddArc(
              output_src,
              Arc(0, arc.olabel, close_weight, output_final));
        }
      }

      const StateId output_dst = map_state(v);
      region_fst->AddArc(
          output_src,
          Arc(arc.ilabel, arc.olabel, weight, output_dst));

      if (final_states.count(v) != 0 && is_binary_type &&
          final_links_added.insert(v).second) {
        region_fst->AddArc(
            output_dst,
            Arc(0, 0, fst.Final(v), output_final));
      }
      Q.push(v);
    }
  }

  return region_fst->NumArcs(output_start) > 0;
}

template <typename StateId>
bool IsSubsetOf(
    const std::set<StateId>& subset,
    const std::vector<StateId>& superset) {
  return std::includes(
      superset.begin(), superset.end(), subset.begin(), subset.end());
}

template <typename StateId>
bool HasVisitedSegment(
    const std::unordered_map<StateId, std::vector<std::vector<StateId>>>&
        visited_segments,
    StateId start_state,
    const std::set<StateId>& end_states) {
  auto it = visited_segments.find(start_state);
  if (it == visited_segments.end()) return false;
  for (const auto& visited_end_states : it->second) {
    if (IsSubsetOf(end_states, visited_end_states)) return true;
  }
  return false;
}

template <typename StateId>
void MarkVisitedSegment(
    std::unordered_map<StateId, std::vector<std::vector<StateId>>>&
        visited_segments,
    StateId start_state,
    const std::set<StateId>& end_states) {
  auto& visited = visited_segments[start_state];
  visited.emplace_back(end_states.begin(), end_states.end());
}

template <typename Weight>
bool StartStatePosteriorMeetsThreshold(
    int32 start_state,
    const std::vector<Weight>& fw,
    const std::vector<Weight>& bw,
    double total_cost,
    double max_neg_log_posterior) {
  if (max_neg_log_posterior == std::numeric_limits<double>::infinity())
    return true;
  const double state_cost =
      fw[start_state].Value() + bw[start_state].Value() - total_cost;
  return state_cost <= max_neg_log_posterior;
}


class LatticeScorerTask {
 public:
  typedef std::unordered_map<fst::StdArc::Label, fst::StdArc::Label> LabelMap;
  typedef std::unordered_set<fst::StdArc::Label> GroupSet;
  typedef std::tuple<std::string, int32, int32, int32, double> output_tuple;
  typedef TableWriter <
  BasicTupleVectorHolder<std::string, int32, int32, int32, double>>
      ScoreWriter;

  LatticeScorerTask(const std::string &key, const CompactLattice &clat,
                    const GroupSet& open_scope_labels,
                    const LabelMap& close_to_open_scope_label,
                    const GroupSet& operator_segmentation_labels,
                    const int32 num_nbest,
                    BaseFloat lattice_timeout,
                    BaseFloat min_posterior,
                    BaseFloat determinize_delta,
                    BaseFloat acoustic_scale, BaseFloat graph_scale,
                    BaseFloat insertion_penalty, BaseFloat beam,
                    ScoreWriter* score_writer) :
      key_(key), clat_(clat), open_scope_labels_(open_scope_labels),
      close_to_open_scope_label_(close_to_open_scope_label),
      operator_segmentation_labels_(operator_segmentation_labels),
      num_nbest_(num_nbest),
      lattice_timeout_(lattice_timeout),
      min_posterior_(min_posterior),
      determinize_delta_(determinize_delta), acoustic_scale_(acoustic_scale),
      graph_scale_(graph_scale), insertion_penalty_(insertion_penalty),
      beam_(beam), score_writer_(score_writer) {}

  ~LatticeScorerTask() {
    if (!discarded_) {
      score_writer_->Write(key_, nbest_segment_);
    } else {
      AddDiscardedLattice();
      KALDI_WARN << "Lattice " << key_ << ": discarded after exceeding "
                 << lattice_timeout_ << " seconds.";
    }
    KALDI_VLOG(1) << "Lattice " << key_ << ": done in " << elapsed_
                  << " seconds.";
    KALDI_VLOG(1) << "Lattice " << key_
                  << ": timing summary: overall=" << overall_elapsed_
                  << "s, indexing=" << elapsed_
                  << "s, SegmentWFST=" << segment_wfst_elapsed_
                  << "s (" << segment_wfst_calls_ << " calls)"
                  << ", ExtractRegionFst=" << extract_region_elapsed_
                  << "s (" << extract_region_calls_ << " calls)"
                  << ", ComputePrIxFST=" << compute_prix_elapsed_
                  << "s (" << compute_prix_calls_ << " calls)"
                  << ", indexed_expressions=" << pos_count_ << ".";
    if (!discarded_) {
      AddAggregateTimingStats(
          overall_elapsed_, elapsed_, segment_wfst_elapsed_,
          extract_region_elapsed_, compute_prix_elapsed_, segment_wfst_calls_,
          extract_region_calls_, compute_prix_calls_, pos_count_);
    }
  }

  void operator()() {
    Timer overall_timer;
    ProcessLattice(key_, &clat_, acoustic_scale_, graph_scale_,
                   insertion_penalty_, beam_);

    if (clat_.Start() == fst::kNoStateId) {
      elapsed_ = 0;
      overall_elapsed_ = overall_timer.Elapsed();
      return;
    }

    // Timer does not include preprocessing time.
    Timer timer;

    // Convert CompactLattice to LogFst with character labels on the input side
    // and frame segments on the output side.
    fst::VectorFst<fst::LogArc> tmp_fst;
    std::vector< std::tuple<int32, int32> > label_to_segm;
    fst::CompactLatticeToWordCountSegmentFst(clat_, &tmp_fst, &label_to_segm);
    clat_.DeleteStates();  // Free memory

    // Compute forward/backward
    fst::vector<fst::LogArc::Weight> fw, bw;
    fst::ShortestDistance(tmp_fst, &fw, /* reverse= */ false);
    fst::ShortestDistance(tmp_fst, &bw, /* reverse= */ true);

    const double total_cost = bw[tmp_fst.Start()].Value();

    std::unordered_map<int32, std::set<int32>> segments;
    std::unordered_map<int32, std::set<int32>> segmentsBin;
    std::unordered_map<int32, ScopeCloseArcs> scope_closes;
    int32 final_state;
    std::set<int32> final_states;

    // Collect delimiter scopes and operator split points from the lattice.
    {
      Timer stage_timer;
      fst::SegmentWFST(segments, segmentsBin, scope_closes, tmp_fst,
                       final_state,
                       final_states, open_scope_labels_,
                       close_to_open_scope_label_,
                       operator_segmentation_labels_);
      segment_wfst_elapsed_ += stage_timer.Elapsed();
      segment_wfst_calls_++;
    }
    if (CheckTimeout(overall_timer, timer)) return;

    int pos = 0;

    std::unordered_map<int32, std::set<int32>> scope_segments = segments;
    scope_segments[tmp_fst.Start()] = final_states;
    std::vector<int32> sorted_scope_keys;
    sorted_scope_keys.reserve(scope_segments.size());
    for (const auto& entry : scope_segments) {
      sorted_scope_keys.push_back(entry.first);
    }
    std::sort(sorted_scope_keys.begin(), sorted_scope_keys.end());
    std::unordered_map<int32, std::vector<std::vector<int32>>> visited_segments;
    std::unordered_map<int32, std::vector<std::vector<int32>>>
        visited_complete_segments;
    const double max_neg_log_posterior =
        min_posterior_ > 0.0
            ? -std::log(min_posterior_)
            : std::numeric_limits<double>::infinity();

    // Single loop over delimiter scopes and operator split points.
    for (const auto beginState : sorted_scope_keys) {
      const auto& endStates = scope_segments.at(beginState);
      if (endStates.empty()) continue;
      if (!StartStatePosteriorMeetsThreshold(beginState, fw, bw, total_cost,
                                             max_neg_log_posterior)) {
        continue;
      }

      const bool is_root_scope = beginState == tmp_fst.Start();
      const auto closing_arcs = scope_closes.find(beginState);
      const auto* scope_boundaries =
          closing_arcs == scope_closes.end() ? nullptr : &closing_arcs->second;
      const bool emit_scope_segment =
          !is_root_scope || segments.count(beginState) != 0;
      auto& scope_visited =
          is_root_scope ? visited_complete_segments : visited_segments;

      if (emit_scope_segment &&
          !HasVisitedSegment(scope_visited, beginState, endStates)) {
        fst::VectorFst<fst::LogArc> fst_subexp;
        bool extracted = false;
        {
          Timer stage_timer;
          extracted = ExtractRegionFst(tmp_fst, beginState, endStates, fw, bw,
                                       0, scope_boundaries, final_states,
                                       operator_segmentation_labels_,
                                       &fst_subexp);
          extract_region_elapsed_ += stage_timer.Elapsed();
          extract_region_calls_++;
        }
        if (extracted) {
          MarkVisitedSegment(scope_visited, beginState, endStates);
          {
            Timer stage_timer;
            ComputePrIxFST(nbest_segment_, fst_subexp, label_to_segm, key_,
                           total_cost, max_neg_log_posterior, pos,
                           determinize_delta_, num_nbest_);
            compute_prix_elapsed_ += stage_timer.Elapsed();
            compute_prix_calls_++;
          }
          if (CheckTimeout(overall_timer, timer)) return;
          pos++;
        }
      }

      auto bin_it = segmentsBin.find(beginState);
      if (bin_it == segmentsBin.end()) continue;

      std::set<int32> bin_splits;
      bin_splits.insert(beginState);
      bin_splits.insert(bin_it->second.begin(), bin_it->second.end());

      if (bin_splits.size() <= 1) continue;

      for (auto it = std::next(bin_splits.begin()); it != bin_splits.end();
           ++it) {
        const auto delim = *it;

        for (auto it_left = bin_splits.begin(); it_left != it; ++it_left) {
          const auto delim_left = *it_left;
          if (delim_left >= delim - 1) continue;
          if (!StartStatePosteriorMeetsThreshold(
                  delim_left, fw, bw, total_cost,
                  max_neg_log_posterior)) {
            continue;
          }
          const std::set<int32> left_end_states{delim};
          if (HasVisitedSegment(visited_segments, delim_left,
                                left_end_states)) {
            continue;
          }

          fst::VectorFst<fst::LogArc> fst_subexp;
          bool extracted = false;
          {
            Timer stage_timer;
            extracted = ExtractRegionFst(tmp_fst, delim_left, left_end_states,
                                         fw, bw, 0, nullptr, final_states,
                                         operator_segmentation_labels_,
                                         &fst_subexp);
            extract_region_elapsed_ += stage_timer.Elapsed();
            extract_region_calls_++;
          }
          if (extracted) {
            MarkVisitedSegment(visited_segments, delim_left, left_end_states);
            {
              Timer stage_timer;
              ComputePrIxFST(nbest_segment_, fst_subexp, label_to_segm, key_,
                             total_cost, max_neg_log_posterior, pos,
                             determinize_delta_, num_nbest_);
              compute_prix_elapsed_ += stage_timer.Elapsed();
              compute_prix_calls_++;
            }
            if (CheckTimeout(overall_timer, timer)) return;
            pos++;
          }
        }

        if (!StartStatePosteriorMeetsThreshold(delim, fw, bw, total_cost,
                                               max_neg_log_posterior)) {
          continue;
        }

        const bool is_complete = is_root_scope;
        auto& right_visited =
            is_complete ? visited_complete_segments : visited_segments;
        std::set<int32> right_end_states;
        for (const auto& endState : endStates) {
          if (delim < endState - 1) right_end_states.insert(endState);
        }
        if (right_end_states.empty() ||
            HasVisitedSegment(right_visited, delim, right_end_states)) {
          continue;
        }

        fst::VectorFst<fst::LogArc> fst_subexp;
        bool extracted = false;
        {
          Timer stage_timer;
          extracted = ExtractRegionFst(tmp_fst, delim, right_end_states, fw,
                                       bw, is_complete, scope_boundaries,
                                       final_states,
                                       operator_segmentation_labels_,
                                       &fst_subexp);
          extract_region_elapsed_ += stage_timer.Elapsed();
          extract_region_calls_++;
        }
        if (extracted) {
          MarkVisitedSegment(right_visited, delim, right_end_states);
          {
            Timer stage_timer;
            ComputePrIxFST(nbest_segment_, fst_subexp, label_to_segm, key_,
                           total_cost, max_neg_log_posterior, pos,
                           determinize_delta_, num_nbest_);
            compute_prix_elapsed_ += stage_timer.Elapsed();
            compute_prix_calls_++;
          }
          if (CheckTimeout(overall_timer, timer)) return;
          pos++;
        }
      }
    }

    // Add transitions from original final states to the synthetic final state.
    for (const auto& fState : final_states) {
      tmp_fst.AddArc(
          fState,
          fst::LogArc(0, 0, tmp_fst.Final(fState), final_state));
      tmp_fst.SetFinal(fState, fst::LogArc::Weight::Zero());
    }

    {
      Timer stage_timer;
      ComputePrIxFST(nbest_segment_, tmp_fst, label_to_segm, key_,
                     total_cost, max_neg_log_posterior, pos,
                     determinize_delta_, num_nbest_);
      compute_prix_elapsed_ += stage_timer.Elapsed();
      compute_prix_calls_++;
    }
    if (CheckTimeout(overall_timer, timer)) return;
    pos++;

    pos_count_ = pos;
    elapsed_ = timer.Elapsed();
    overall_elapsed_ = overall_timer.Elapsed();
  }

 private:
  bool CheckTimeout(const Timer& overall_timer, const Timer& indexing_timer) {
    if (lattice_timeout_ <= 0.0) return false;
    if (overall_timer.Elapsed() <= lattice_timeout_) return false;
    discarded_ = true;
    nbest_segment_.clear();
    pos_count_ = 0;
    elapsed_ = indexing_timer.Elapsed();
    overall_elapsed_ = overall_timer.Elapsed();
    return true;
  }

  const std::string key_;
  CompactLattice clat_;
  GroupSet open_scope_labels_;
  LabelMap close_to_open_scope_label_;
  GroupSet operator_segmentation_labels_;
  const int32 num_nbest_;
  BaseFloat lattice_timeout_;
  BaseFloat min_posterior_;
  BaseFloat determinize_delta_;
  BaseFloat acoustic_scale_, graph_scale_, insertion_penalty_, beam_;
  std::vector<output_tuple> nbest_segment_;
  ScoreWriter* score_writer_;
  double elapsed_;
  double overall_elapsed_ = 0.0;
  double segment_wfst_elapsed_ = 0.0;
  double extract_region_elapsed_ = 0.0;
  double compute_prix_elapsed_ = 0.0;
  int64 segment_wfst_calls_ = 0;
  int64 extract_region_calls_ = 0;
  int64 compute_prix_calls_ = 0;
  int32 pos_count_ = 0;
  bool discarded_ = false;
};

}  // namespace kaldi


int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    typedef fst::StdArc::Label Label;

    const char* usage =
        "Build a position-level word index from character lattices.\n"
        "\n"
        "Math delimiters are resolved from --symbol-table. Opening and "
        "closing delimiter symbols define nested scopes, while operator "
        "prefixes add segmentation points inside those scopes.\n"
        "\n"
        "Usage: lattice-char-index-position-math --symbol-table=words.txt "
        "[options] <lattice-rspecifier> <index-wspecifier>\n"
        " e.g.: lattice-char-index-position-math --symbol-table=words.txt "
        "ark:lats.ark ark,t:-\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat determinize_delta = fst::kDelta / 8.0F;
    BaseFloat lattice_timeout = 300.0;
    BaseFloat min_posterior = 1.0e-15;
    int nbest = 100;
    std::string symbol_table_rxfilename = "";
    std::string open_delimiters_str = "{ {@";
    std::string close_delimiters_str = "} }@";
    std::string operator_segment_prefixes_str = "@bin @rel @punct";
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam,
                "Pruning beam (applied after acoustic scaling and adding the "
                "insertion penalty).");
    po.Register("nbest", &nbest, "Extract this number of n-best hypotheses.");
    po.Register("determinize-delta", &determinize_delta,
                "Delta threshold used for the determinization algorithm.");
    po.Register("lattice-timeout", &lattice_timeout,
                "Discard lattices whose total processing time exceeds this "
                "number of seconds. Use 0 or a negative value to disable.");
    po.Register("min-posterior", &min_posterior,
                "Skip regions whose start-state posterior is below this "
                "threshold and index entries whose posterior is at or below "
                "it. Use 0 or a negative value to disable both filters. "
                "Default is a conservative 1e-15.");
    po.Register("symbol-table", &symbol_table_rxfilename,
                "Symbol table used to resolve math delimiter labels.");
    po.Register("open-delimiters", &open_delimiters_str,
                "Space-separated opening delimiter symbols.");
    po.Register("close-delimiters", &close_delimiters_str,
                "Space-separated closing delimiter symbols. The order must "
                "match --open-delimiters.");
    po.Register("operator-segment-prefixes", &operator_segment_prefixes_str,
                "Space-separated symbol prefixes whose labels add operator "
                "segmentation points inside delimiter scopes.");
    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);

    po.Read(argc, argv);
    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    if (symbol_table_rxfilename.empty()) {
      KALDI_ERR << "--symbol-table is required for resolving math "
                << "delimiter labels.";
    }
    std::unique_ptr<fst::SymbolTable> symbol_table(
        fst::SymbolTable::ReadText(symbol_table_rxfilename));
    if (!symbol_table) {
      KALDI_ERR << "Could not read symbol table: "
                << symbol_table_rxfilename;
    }
    const auto open_scope_label_vector =
        ResolveSymbolVector<Label>(*symbol_table, open_delimiters_str);
    const auto close_scope_label_vector =
        ResolveSymbolVector<Label>(*symbol_table, close_delimiters_str);
    const auto open_scope_labels = ToLabelSet(open_scope_label_vector);
    const auto close_to_open_scope_label =
        BuildCloseToOpenDelimiterMap(open_scope_label_vector,
                                     close_scope_label_vector);
    const auto operator_segmentation_labels =
        ResolvePrefixedSymbols<Label>(*symbol_table,
                                      operator_segment_prefixes_str);

    Timer timer;
    int num_lats = 0;
    LatticeScorerTask::ScoreWriter score_writer(po.GetArg(2));
    TaskSequencer<LatticeScorerTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(po.GetArg(1));
         !lattice_reader.Done(); lattice_reader.Next(), ++num_lats) {
      const std::string lattice_key = lattice_reader.Key();
      // Schedule task.
      task_sequencer.Run(new LatticeScorerTask(lattice_reader.Key(),
                                               lattice_reader.Value(),
                                               open_scope_labels,
                                               close_to_open_scope_label,
                                               operator_segmentation_labels,
                                               nbest,
                                               lattice_timeout,
                                               min_posterior,
                                               determinize_delta,
                                               acoustic_scale,
                                               graph_scale,
                                               insertion_penalty,
                                               beam,
                                               &score_writer));
      lattice_reader.FreeCurrent();
    }
    task_sequencer.Wait();
    score_writer.Close();
    LogAggregateTimingStats();

    const auto effective_num_threads =
        std::min<size_t>(num_lats, task_sequencer_config.num_threads);
    KALDI_LOG << "Done " << num_lats << " lattices in "  << timer.Elapsed()
              << "s using " << effective_num_threads << " threads.";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}
