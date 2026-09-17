#ifndef SVPCHANGE_SVP_TESTS_RANK_H
#define SVPCHANGE_SVP_TESTS_RANK_H

#include "svp_test_base.h"
#include "svp_tests_quantile.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

/// GvarCost cost:
/// C_n = 1/n sum (y_t - mean)^2, statistic = C_n / sigma^2 ~ chi^2_{n-1} under H0.
class VarianceCost : public TestBase
{
public:
  VarianceCost()
    : n_(0),
      mean_(0.0),
      M2_(0.0)   // will store sum (y_t - mean)^2
  {}

  void update(double y) override
  {
    ++n_;
    // Welford's online update
    double delta  = y - mean_;
    mean_        += delta / static_cast<double>(n_);
    double delta2 = y - mean_;
    M2_          += delta * delta2; // accumulated SSE around the running mean
  }

  // Returns chi-square statistic: C_n / sigma^2
  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();

    return M2_/n_;
  }


private:
  int    n_;      // number of points
  double mean_;   // running mean
  double M2_;     // running sum of squared deviations from mean (SSE)
};





/// WilcoxonCost:
/// - Stores the data y_1, ..., y_n via update(y).
/// - update() maintains, for each split t = 1..n-1,
///     A = {y_1, ..., y_t},  B = {y_{t+1}, ..., y_n}
///   the Wilcoxon rank-sum statistic W_t (centered),
///   and returns max_t |W_t|.
///
/// New pairwise rank contributions update all existing splits in O(n),
/// including exact half-weight handling of ties. statistic() is O(1).

class WilcoxonCost : public TestBase
{
public:
  WilcoxonCost() : doubled_statistic_(0) {}

  void update(double y) override
  {
    const std::size_t old_n = values_.size();
    std::int64_t contribution = 0;
    for (std::size_t i = 0; i < old_n; ++i) {
      contribution += static_cast<std::int64_t>(values_[i] > y) -
        static_cast<std::int64_t>(values_[i] < y);

      // Existing split t = i + 1 gains the contribution of y to its
      // right-hand group. The final contribution is the new last split.
      if (i < split_statistics_.size()) {
        split_statistics_[i] += contribution;
      }
    }

    values_.push_back(y);
    if (old_n > 0) {
      split_statistics_.push_back(contribution);
    }

    doubled_statistic_ = 0;
    for (std::int64_t value : split_statistics_) {
      const std::int64_t magnitude = value < 0 ? -value : value;
      doubled_statistic_ = std::max(doubled_statistic_, magnitude);
    }
  }

  double statistic() const override
  {
    if (values_.empty()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return 0.5 * static_cast<double>(doubled_statistic_);
  }

private:
  std::vector<double> values_;
  std::vector<std::int64_t> split_statistics_;
  std::int64_t doubled_statistic_;
};


/**
 * MedianMoodCost:
 *
 * Implements a scan version of Mood's median test for a single changepoint.
 *
 * - Data are added via update(y).
 * - statistic() computes:
 *      * global median m of all y's,
 *      * for each split t, counts below/above m in the two groups,
 *      * builds a 2x2 table and computes chi-square,
 *      * returns the maximum chi-square over t.
 *
 * This is robust (median + only above/below information) and
 * detects changes in location (median).
 */
class MedianMoodCost : public TestBase {
public:
  MedianMoodCost() = default;

  void update(double y) override {
    values_.push_back(y);
    tree_.insert(y);
  }

  double statistic() const override {
    const std::size_t n = values_.size();
    if (n < 2) {
      return 0.0;  // not enough data to split
    }

    // The maintained tree returns the same upper median as nth_element(n / 2)
    // without copying and partitioning the complete segment on every call.
    const double med = tree_.kth(static_cast<int>(n / 2));

    const int total_below = tree_.count_less(med);
    const int total_above = static_cast<int>(n) - total_below -
      tree_.count_equal(med);

    const int N_effective   = total_below + total_above;

    // If everything is exactly equal to the median, no information.
    if (N_effective == 0) {
      return 0.0;
    }

    // --- 3) Scan all splits and compute chi-square ---

    double best_chisq = 0.0;
    int prefix_below = 0;
    int prefix_above = 0;

    for (std::size_t t = 1; t < n; ++t) {
      const double value = values_[t - 1];
      if (value < med) ++prefix_below;
      else if (value > med) ++prefix_above;

      // Group A: indices [0, t-1]
      // Group B: indices [t, n-1]
      int a11 = prefix_below;               // A, below
      int a12 = prefix_above;               // A, above
      int a21 = total_below - a11;         // B, below
      int a22 = total_above - a12;         // B, above

      int nA = a11 + a12;
      int nB = a21 + a22;

      // no "signal" in one side -> skip
      if (nA == 0 || nB == 0) {
        continue;
      }

      if (total_below == 0 || total_above == 0) continue;

      // Algebraically identical Pearson chi-square formula for a 2x2 table.
      const double determinant =
        static_cast<double>(a11) * a22 - static_cast<double>(a12) * a21;
      const double N = static_cast<double>(nA + nB);
      const double denominator = static_cast<double>(nA) * nB *
        total_below * total_above;
      const double chisq = N * determinant * determinant / denominator;

      if (chisq > best_chisq) {
        best_chisq = chisq;
      }
    }

    return best_chisq;
  }

private:
  std::vector<double> values_;
  OrderStatisticTree tree_;
};

#endif
