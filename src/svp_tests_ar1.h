#ifndef SVPCHANGE_SVP_TESTS_AR1_H
#define SVPCHANGE_SVP_TESTS_AR1_H

#include "ar1_focus_state.h"
#include "svp_test_base.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Approximate fixed-rho AR(1) FOCuS likelihood recurrence.
////////////////////////////////////////////////////////////////////////////////

class AR1FocusMeanChange : public TestBase
{
public:
  explicit AR1FocusMeanChange(double rho, double sigma2 = 1.0)
    : sigma2_(sigma2), observation_count_(0),
      info_(rho) {}

  void update(double y) override
  {
    ++observation_count_;
    info_.update(y);
  }

  double statistic() const override
  {
    if (observation_count_ < 4) return 0.0;
    // Unified AR FOCuS returns the conventional 2 log-likelihood ratio;
    // SVP's gamma and exact AR1 test use the half-log-likelihood scale.
    return std::max(0.0, 0.5 * info_.max_statistic() / sigma2_);
  }


private:
  double sigma2_;
  std::size_t observation_count_;
  changepoint::AR1FocusState info_;
};

////////////////////////////////////////////////////////////////////////////////
// Conditional fixed-rho likelihood scan for a change in the marginal AR(1)
// mean, conditional on the first observation.
//
// For a change after observation tau, the transition tau -> tau + 1 has residual
//   x[tau + 1] - rho*x[tau] - mu2 + rho*mu1,
// so it cannot be represented by an ordinary two-mean split of the innovations.
// Prefix sums make each candidate tau O(1); scanning all tau is O(n). The
// initial stationary residual likelihood is omitted, consistently with the
// conditional innovation likelihood used by the R validity reference.
////////////////////////////////////////////////////////////////////////////////

class AR1ExactMeanChange : public TestBase
{
public:
  AR1ExactMeanChange(double rho, double sigma2 = 1.0, bool profile_sigma = false)
    : rho_(rho),
      sigma2_(sigma2),
      profile_sigma_(profile_sigma),
      previous_(0.0),
      has_previous_(false),
      observation_count_(0),
      cached_statistic_(0.0),
      cached_rss0_(0.0),
      cached_rss1_(0.0),
      cached_changepoint_(-1),
      dirty_(false)
  {
    prefix_sum_.push_back(0.0);
    prefix_sq_.push_back(0.0);
  }

  void update(double y) override
  {
    ++observation_count_;
    if (has_previous_) {
      const double innovation = y - rho_ * previous_;
      innovations_.push_back(innovation);
      prefix_sum_.push_back(prefix_sum_.back() + innovation);
      prefix_sq_.push_back(prefix_sq_.back() + innovation * innovation);
    }
    previous_ = y;
    has_previous_ = true;
    dirty_ = true;
  }

  double statistic() const override
  {
    compute();
    return cached_statistic_;
  }

  int changepoint() const
  {
    compute();
    return cached_changepoint_;
  }

  double rss0() const
  {
    compute();
    return cached_rss0_;
  }

  double rss1() const
  {
    compute();
    return cached_rss1_;
  }

private:
  void compute() const
  {
    if (!dirty_) return;

    cached_statistic_ = 0.0;
    cached_rss0_ = 0.0;
    cached_rss1_ = 0.0;
    cached_changepoint_ = -1;

    // Candidate changes leave at least two observations on each side.
    if (observation_count_ < 4) {
      dirty_ = false;
      return;
    }

    const std::size_t n = observation_count_;
    const std::size_t transition_count = n - 1;
    const double total_sum = prefix_sum_[transition_count];
    const double total_sq = prefix_sq_[transition_count];
    cached_rss0_ = std::max(
      0.0,
      total_sq - total_sum * total_sum /
        static_cast<double>(transition_count)
    );

    const double one_minus_rho = 1.0 - rho_;
    double best_rss = std::numeric_limits<double>::infinity();
    int best_tau = -1;

    for (std::size_t tau = 2; tau <= n - 2; ++tau) {
      const std::size_t pre_count = tau - 1;
      const std::size_t post_count = n - tau - 1;
      const double pre_sum = prefix_sum_[pre_count];
      const double transition = innovations_[tau - 1];
      const double post_sum = total_sum - prefix_sum_[tau];

      const double a =
        static_cast<double>(pre_count) * one_minus_rho * one_minus_rho +
        rho_ * rho_;
      const double b = -rho_;
      const double d =
        1.0 +
        static_cast<double>(post_count) * one_minus_rho * one_minus_rho;
      const double v1 = one_minus_rho * pre_sum - rho_ * transition;
      const double v2 = transition + one_minus_rho * post_sum;
      const double determinant = a * d - b * b;

      if (determinant <= std::numeric_limits<double>::epsilon()) continue;

      const double fitted =
        (d * v1 * v1 - 2.0 * b * v1 * v2 + a * v2 * v2) /
        determinant;
      const double rss = std::max(0.0, total_sq - fitted);
      if (rss < best_rss) {
        best_rss = rss;
        best_tau = static_cast<int>(tau);
      }
    }

    if (best_tau < 0) {
      dirty_ = false;
      return;
    }

    cached_rss1_ = best_rss;
    cached_changepoint_ = best_tau;
    if (profile_sigma_) {
      if (cached_rss0_ > 0.0 && best_rss <=
          std::numeric_limits<double>::epsilon() * cached_rss0_) {
        cached_statistic_ = std::numeric_limits<double>::infinity();
      } else if (best_rss > 0.0 && cached_rss0_ > best_rss) {
        // Use the log-likelihood-ratio scale (FOCuS uses this convention).
        cached_statistic_ =
          0.5 * static_cast<double>(transition_count) *
          std::log(cached_rss0_ / best_rss);
      }
    } else {
      // Half the usual 2 log-likelihood ratio, matching Gaussian FOCuS.
      cached_statistic_ =
        std::max(0.0, (cached_rss0_ - best_rss) / (2.0 * sigma2_));
    }
    dirty_ = false;
  }

  double rho_;
  double sigma2_;
  bool profile_sigma_;
  double previous_;
  bool has_previous_;
  std::size_t observation_count_;
  std::vector<double> innovations_;
  std::vector<double> prefix_sum_;
  std::vector<double> prefix_sq_;
  mutable double cached_statistic_;
  mutable double cached_rss0_;
  mutable double cached_rss1_;
  mutable int cached_changepoint_;
  mutable bool dirty_;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
