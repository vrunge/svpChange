#pragma once
#include "focus.h"
#include <memory>

class TestBase
{
  public:
    virtual ~TestBase() = default;
    virtual void update(double y) = 0;
    virtual double statistic() const = 0;
};

class GaussianMean : public TestBase
{
  public:
    std::unique_ptr<Info> info;

    GaussianMean()
    {
      auto newP = [](double St, int tau, double m0){
        std::unique_ptr<Piece> p = std::make_unique<PieceGau>();
        p->St = St;
        p->tau = tau;
        p->m0 = m0;
        return p;
      };
      info = std::make_unique<Info>(newP, NAN);
    }

    void update(double y) override
    {
      info->update(y);
    }

    double statistic() const override
    {
      return info->statistic();
    }
};
