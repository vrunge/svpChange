#ifndef SVPCHANGE_GAUSSIAN_FOCUS_STATE_H
#define SVPCHANGE_GAUSSIAN_FOCUS_STATE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace changepoint {

struct GaussianFocusCandidate
{
  double cumulative_sum;
  double boundary;
};

// One side of the Gaussian FOCuS functional-pruning recurrence.
class OneSidedGaussianFocusState
{
public:
  OneSidedGaussianFocusState(double cumulative_sum, double observation_count,
                             bool right_side)
    : cumulative_sum_(cumulative_sum),
      observation_count_(observation_count),
      right_side_(right_side),
      active_end_(1)
  {
    candidates_.reserve(30);
    candidates_.push_back({cumulative_sum_, observation_count_});
  }

  void update(double value)
  {
    observation_count_ += 1.0;
    cumulative_sum_ += value;
    prune_by_monotonicity();
    append_candidate();
  }

  double max_cost() const
  {
    double best = 0.0;
    const double total_cost = observation_count_ > 0.0 ?
      cumulative_sum_ * cumulative_sum_ / observation_count_ : 0.0;

    for (std::size_t i = 0; i < active_end_; ++i) {
      const GaussianFocusCandidate& candidate = candidates_[i];
      const double right_length =
        observation_count_ - candidate.boundary;
      if (candidate.boundary <= 0.0 || right_length <= 0.0) continue;

      const double right_sum =
        cumulative_sum_ - candidate.cumulative_sum;
      const double cost =
        candidate.cumulative_sum * candidate.cumulative_sum /
          candidate.boundary +
        right_sum * right_sum / right_length -
        total_cost;
      best = std::max(best, cost);
    }
    return best;
  }

  bool exceeds_cost(double threshold) const
  {
    const double total_cost = observation_count_ > 0.0 ?
      cumulative_sum_ * cumulative_sum_ / observation_count_ : 0.0;
    if (!std::isfinite(total_cost)) return true;

    for (std::size_t i = 0; i < active_end_; ++i) {
      const GaussianFocusCandidate& candidate = candidates_[i];
      const double right_length =
        observation_count_ - candidate.boundary;
      if (candidate.boundary <= 0.0 || right_length <= 0.0) continue;

      const double right_sum =
        cumulative_sum_ - candidate.cumulative_sum;
      const double cost =
        candidate.cumulative_sum * candidate.cumulative_sum /
          candidate.boundary +
        right_sum * right_sum / right_length -
        total_cost;
      if (!std::isfinite(cost) || cost >= threshold) return true;
    }
    return false;
  }

private:
  void prune_by_monotonicity()
  {
    while (active_end_ > 1) {
      const GaussianFocusCandidate& newest = candidates_[active_end_ - 1];
      const GaussianFocusCandidate& previous = candidates_[active_end_ - 2];
      const double newest_mean =
        (cumulative_sum_ - newest.cumulative_sum) /
        (observation_count_ - newest.boundary);
      const double previous_mean =
        (cumulative_sum_ - previous.cumulative_sum) /
        (observation_count_ - previous.boundary);
      const bool prune = right_side_ ?
        newest_mean <= previous_mean :
        newest_mean >= previous_mean;
      if (!prune) break;
      --active_end_;
    }
  }

  void append_candidate()
  {
    const GaussianFocusCandidate candidate{
      cumulative_sum_, observation_count_
    };
    if (active_end_ < candidates_.size()) {
      candidates_[active_end_] = candidate;
    } else {
      candidates_.push_back(candidate);
    }
    ++active_end_;
  }

  double cumulative_sum_;
  double observation_count_;
  bool right_side_;
  std::vector<GaussianFocusCandidate> candidates_;
  std::size_t active_end_;
};

// Two-sided Gaussian FOCuS state used by GaussianMean.
class GaussianFocusState
{
public:
  GaussianFocusState()
    : right_(0.0, 0.0, true),
      left_(0.0, 0.0, false) {}

  void update(double value)
  {
    right_.update(value);
    left_.update(value);
  }

  double max_cost() const
  {
    return std::max(right_.max_cost(), left_.max_cost());
  }

  bool exceeds_cost(double threshold) const
  {
    return right_.exceeds_cost(threshold) ||
      left_.exceeds_cost(threshold);
  }

private:
  OneSidedGaussianFocusState right_;
  OneSidedGaussianFocusState left_;
};

} // namespace changepoint

#endif
