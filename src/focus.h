/*
 This code is sourced from the R package `gtromano/exfocus.rcpp` for change point detection.
 See https://github.com/gtromano/exfocus.rcpp for more information.

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Copyright (C) 2024 Gaetano Romano

 This and all the code in this package is publicly released
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3
 See DESCRIPTION for further information about the authors.
 */

#ifndef ___FCOST_H___
#define ___FCOST_H___

#include <iostream>
//#include <list>
#include <vector>
#include <numeric>      // std::iota
#include <cmath>
//#include <algorithm>
#include <functional>
#include <memory>

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

struct PieceGau:Piece
{
  double eval (const CUSUM& cs, double x, const double& theta0) const
  {
    auto c = (double)(cs.n - tau);
    auto S = (cs.Sn - St);

    // std::cout << "running gaussian" << std::endl;

    if (std::isnan(theta0))
      return - 0.5 * c * x * x + S * x + m0;
    else
      return - 0.5 * c * x * x + S * x;

  }
};

struct PieceBer:Piece
{
  double eval (const CUSUM& cs, double x, const double& theta0) const
  {
    auto c = (double)(cs.n - tau);
    auto S = (cs.Sn - St);

    // std::cout << "running bernoulli" << std::endl;

    if (std::isnan(theta0))
      return S * log(x) + (c - S) * log((1-x)) + m0;
    else
      return S * log(x/theta0) + (c - S) * log((1-x) / (1 - theta0));
  }

  // this is to avoid nans that might be quite annoying in comparisons
  double argmax (const CUSUM &cs ) const
  {
    auto agm = (cs.Sn - St) / (double)(cs.n - tau);
    if (agm == 0) {
      return 0.000000001;
    } else if (agm == 1) {
      return 0.99999999;
    } else {
      return agm;
    }
  }

};

struct PiecePoi:Piece
{
  double eval (const CUSUM& cs, double x, const double& theta0) const {
    auto c = (double)(cs.n - tau);
    auto S = (cs.Sn - St);

    if (std::isnan(theta0))
      return - c * (x) + S * log(x) + m0;
    else
      return -c * (x - theta0) + S * log(x/theta0);

  }

  // this is to avoid nans that might be quite annoying in comparisons
  double argmax (const CUSUM &cs ) const
  {
    auto agm = (cs.Sn - St) / (double)(cs.n - tau);
    if (agm == 0) {
      return 0.000000001;
    } else {
      return agm;
    }
  }
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

struct Info {
  CUSUM cs;
  Cost Ql;
  Cost Qr;
  std::function<std::unique_ptr<Piece>(double, int, double)> newP;
  double theta0;

  void update(const double& y); // Remove adp_max_check
  double statistic() const { return std::max(Ql.opt, Qr.opt); }

  Info(std::function<std::unique_ptr<Piece>(double, int, double)> newP_, double theta0_)
    : newP(newP_), theta0(theta0_) {
    std::vector<std::unique_ptr<Piece>> initpsl;
    for (auto i = 0; i < 50; i++) {
      initpsl.push_back(newP(0.0, 0, 0.0));
    }
    std::vector<std::unique_ptr<Piece>> initpsr;
    for (auto i = 0; i < 50; i++) {
      initpsr.push_back(newP(0.0, 0, 0.0));
    }
    Cost initQl(std::move(initpsl), 0.0, 0);
    Cost initQr(std::move(initpsr), 0.0, 0);
    Ql = std::move(initQl);
    Qr = std::move(initQr);

    CUSUM initcusum;
    cs = initcusum;
  }
};

int get_tau_max (const Cost&, const CUSUM&, const double& , const double&);
// double argmax (const Piece &q, const CUSUM &cs );
// double argmax (const PieceGam &q, const CUSUM &cs );
//
// double eval (const PieceGau &q, const CUSUM &cs, double x, const double &theta0);
// double eval (const PieceBer &q, const CUSUM &cs, double x, const double &theta0);
// void prune (Cost &Q, const CUSUM &cs, const double &theta0, std::function<bool(Piece, Piece)> cond);
//Info focus_step (Info I, const double& y, std::function<std::unique_ptr<Piece>(double, int, double)> newP, const double& thres, const double& theta0, const bool& adp_max_check);

#endif
