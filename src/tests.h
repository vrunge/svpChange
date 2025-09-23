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

// Gamma change-in-rate test
class GammaRate : public TestBase
{
  public:
    std::unique_ptr<Info> info;

    GammaRate(double shape = 1.0)
    {
      auto newP = [shape](double St, int tau, double m0){
        auto p = std::make_unique<PieceGam>();
        p->St = St;
        p->tau = tau;
        p->m0 = m0;
        p->set_shape(shape);
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

class GaussianVariance : public TestBase
{
  public:
    std::unique_ptr<GammaRate> gamma_rate;

    GaussianVariance(double shape = 1.0)
    {
      gamma_rate = std::make_unique<GammaRate>(shape);
    }

    void update(double y) override
    {
      gamma_rate->update(y * y); // Use squared data
    }

    double statistic() const override
    {
      return gamma_rate->statistic();
    }
};