/*
 This code is sourced from the R package `gtromano/exfocus.rcpp` for change point detection.
 See https://github.com/gtromano/exfocus.rcpp for more information.

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Copyright (C) 2024 Gaetano Romano

 This and all the code in this package is publicly released
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3
 See DESCRIPTION for further information about the authors.
 */

#ifndef SVPCHANGE_GAMMA_FOCUS_H
#define SVPCHANGE_GAMMA_FOCUS_H

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

struct CUSUM {
  double Sn = 0.0;
  int n = 0.0;
};


struct Piece
{
  double St = 0.0;
  int tau = 0;
  double m0 = 0.0;
  double Mdiff = 0.0;

  virtual ~Piece() = default;

  // eval method has no generic, and it's distribution specific
  virtual double eval (const CUSUM& cs, double x, const double& theta0) const = 0;

  // this is the generic argmax method, that should work for all but gamma
  virtual double argmax (const CUSUM &cs ) const {
    return (cs.Sn - St) / (double)(cs.n - tau);
  }

  virtual void set_shape (double s) {
    throw("This is only to set a gamma shape");
  };

};

struct PieceGam:Piece {
  double shape = 1.0;
  double eval (const CUSUM& cs, double x, const double& theta0) const {
    auto c = (double)(cs.n - tau);
    auto S = (cs.Sn - St);

    if (std::isnan(theta0))
      return -c * shape * log(x) - S * (1/x) + m0;
    else
      return c * shape * log(theta0/x) - S * (1/x - 1/theta0);

  }

  void set_shape(double s) {
    shape = s;
  }

  double argmax (const CUSUM &cs ) const {
    return (cs.Sn - St) / ( shape * (double)(cs.n - tau));
  }

};

// the cost is a list of shared pointers to pieces of type Piece
struct Cost {
  std::vector<std::unique_ptr<Piece>> ps;
  double opt = 0;
  long unsigned int k = 0;

  // Add this constructor:
  Cost(std::vector<std::unique_ptr<Piece>>&& ps_, double opt_, long unsigned int k_)
    : ps(std::move(ps_)), opt(opt_), k(k_) {}
  Cost() = default;
};

// Incremental two-sided FOCuS state used only by the Gamma-rate test.
// The specific name avoids confusion with changepoint::Info.
struct GammaFocusInfo {
  CUSUM cs;
  Cost Ql;
  Cost Qr;
  std::function<std::unique_ptr<Piece>(double, int, double)> newP;
  double theta0;

  void update(const double& y); // Remove adp_max_check
  double statistic() const { return std::max(Ql.opt, Qr.opt); }

  GammaFocusInfo(
      std::function<std::unique_ptr<Piece>(double, int, double)> newP_,
      double theta0_)
    : newP(newP_), theta0(theta0_) {
    std::vector<std::unique_ptr<Piece>> initpsl;
    initpsl.reserve(16);
    initpsl.push_back(newP(0.0, 0, 0.0));
    std::vector<std::unique_ptr<Piece>> initpsr;
    initpsr.reserve(16);
    initpsr.push_back(newP(0.0, 0, 0.0));
    Cost initQl(std::move(initpsl), 0.0, 0);
    Cost initQr(std::move(initpsr), 0.0, 0);
    Ql = std::move(initQl);
    Qr = std::move(initQr);

    CUSUM initcusum;
    cs = initcusum;
  }
};

#endif
