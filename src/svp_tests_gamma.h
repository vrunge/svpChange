#ifndef SVPCHANGE_SVP_TESTS_GAMMA_H
#define SVPCHANGE_SVP_TESTS_GAMMA_H

#include "gamma_focus.h"
#include "svp_test_base.h"
#include <cmath>
// Gamma change-in-rate test
class GammaRate : public TestBase
{
  public:
    GammaFocusInfo info;

    explicit GammaRate(double shape = 1.0)
      : info([shape](double St, int tau, double m0) {
        auto p = std::make_unique<PieceGam>();
        p->St = St;
        p->tau = tau;
        p->m0 = m0;
        p->set_shape(shape);
        return p;
      }, NAN) {}

    void update(double y) override
    {
      info.update(y);
    }

    double statistic() const override
    {
      return info.statistic();
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
