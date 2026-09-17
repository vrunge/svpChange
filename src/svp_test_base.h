#ifndef SVPCHANGE_SVP_TEST_BASE_H
#define SVPCHANGE_SVP_TEST_BASE_H

#include <cmath>

class TestBase
{
  public:
    virtual ~TestBase() = default;
    virtual void update(double y) = 0;
    virtual double statistic() const = 0;

    // Validity-only interface. Tests with an incremental threshold check can
    // override this to stop before computing a full statistic.
    virtual bool passes(double gamma) const
    {
      const double value = statistic();
      return std::isfinite(value) && value < gamma;
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
