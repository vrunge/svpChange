#ifndef SVPCHANGE_SVP_TESTS_GAUSSIAN_H
#define SVPCHANGE_SVP_TESTS_GAUSSIAN_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"
#include "svp_tests_gamma.h"
class GaussianMean : public TestBase
{
  public:
    changepoint::GaussianFocusState info;

    GaussianMean() = default;

    void update(double y) override
    {
      info.update(y);
    }

    double statistic() const override
    {
      return info.max_cost();
    }

    bool passes(double gamma) const override
    {
      return !info.exceeds_cost(gamma);
    }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class GaussianVariance : public TestBase
{
  public:
    GammaRate gamma_rate;

    explicit GaussianVariance(double shape = 1.0)
      : gamma_rate(shape) {}

    void update(double y) override
    {
      gamma_rate.update(y * y); // Use squared data
    }

    double statistic() const override
    {
      return gamma_rate.statistic();
    }
};

////////////////////////////////////////////////////////////////////////////////

#endif
