#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

class RunningMedianMood {
public:
  void update(double value) {
    values_.push_back(value);
    if (upper_.empty() || value >= *upper_.begin()) {
      upper_.insert(value);
    } else {
      lower_.insert(value);
    }
    rebalance();
  }

  double statistic() const {
    const std::size_t n = values_.size();
    if (n < 2) return 0.0;

    const double median = *upper_.begin();
    int total_below = 0;
    int total_above = 0;
    for (double value : values_) {
      total_below += value < median;
      total_above += value > median;
    }
    const int effective_n = total_below + total_above;
    if (effective_n == 0 || total_below == 0 || total_above == 0) {
      return 0.0;
    }

    int prefix_below = 0;
    int prefix_above = 0;
    double best = 0.0;
    for (std::size_t split = 1; split < n; ++split) {
      const double value = values_[split - 1];
      prefix_below += value < median;
      prefix_above += value > median;
      const int n_a = prefix_below + prefix_above;
      const int n_b = effective_n - n_a;
      if (n_a == 0 || n_b == 0) continue;

      const int determinant =
        prefix_below * (total_above - prefix_above) -
        prefix_above * (total_below - prefix_below);
      const double denominator = static_cast<double>(n_a) * n_b *
        total_below * total_above;
      const double value_chisq = static_cast<double>(effective_n) *
        determinant * determinant / denominator;
      best = std::max(best, value_chisq);
    }
    return best;
  }

private:
  void rebalance() {
    const std::size_t target_lower = values_.size() / 2;
    while (lower_.size() > target_lower) {
      auto item = std::prev(lower_.end());
      upper_.insert(*item);
      lower_.erase(item);
    }
    while (lower_.size() < target_lower) {
      auto item = upper_.begin();
      lower_.insert(*item);
      upper_.erase(item);
    }
  }

  std::vector<double> values_;
  std::multiset<double> lower_;
  std::multiset<double> upper_;
};

class RunningWilcoxon {
public:
  void update(double value) {
    const std::size_t old_n = values_.size();
    long long contribution = 0;
    for (std::size_t i = 0; i < old_n; ++i) {
      contribution += static_cast<long long>(values_[i] > value) -
        static_cast<long long>(values_[i] < value);
      if (i < split_statistics_.size()) {
        split_statistics_[i] += contribution;
      }
    }
    values_.push_back(value);
    if (old_n > 0) split_statistics_.push_back(contribution);

    doubled_statistic_ = 0;
    for (long long candidate : split_statistics_) {
      doubled_statistic_ = std::max(
        doubled_statistic_, candidate < 0 ? -candidate : candidate
      );
    }
  }

  double statistic() const {
    return 0.5 * static_cast<double>(doubled_statistic_);
  }

private:
  std::vector<double> values_;
  std::vector<long long> split_statistics_;
  long long doubled_statistic_ = 0;
};

double prefix_max(const Rcpp::NumericVector& data, bool reverse) {
  RunningMedianMood test;
  double best = 0.0;
  const R_xlen_t n = data.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const R_xlen_t index = reverse ? n - i - 1 : i;
    test.update(data[index]);
    if (i > 0) best = std::max(best, test.statistic());
  }
  return best;
}

double wilcoxon_prefix_max(const Rcpp::NumericVector& data, bool reverse) {
  RunningWilcoxon test;
  double best = 0.0;
  const R_xlen_t n = data.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const R_xlen_t index = reverse ? n - i - 1 : i;
    test.update(data[index]);
    if (i > 0) best = std::max(best, test.statistic());
  }
  return best;
}

// [[Rcpp::export]]
double well_log_mood_right_critical_cpp(Rcpp::NumericVector data) {
  return prefix_max(data, false);
}

// [[Rcpp::export]]
double well_log_mood_bidirectional_critical_cpp(Rcpp::NumericVector data) {
  return std::max(prefix_max(data, false), prefix_max(data, true));
}

// [[Rcpp::export]]
double well_log_wilcoxon_right_critical_cpp(Rcpp::NumericVector data) {
  return wilcoxon_prefix_max(data, false);
}
