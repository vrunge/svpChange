/*
 This code is sourced from the R package `gtromano/exfocus.rcpp` for change point detection.
 See https://github.com/gtromano/exfocus.rcpp for more information.

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Copyright (C) 2024 Gaetano Romano

 This and all the code in this package is publicly released
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3
 See DESCRIPTION for further information about the authors.
 */

# include "focus.h"

/*
 * PRUNING FUNCTION checks and removes quadratics that are no longer optimal
 */

void prune (Cost& Q, const CUSUM& cs, const double& theta0, const bool isRight)
{
  // if we only have one piece, exit
  if (Q.k == 0)
    return ;


  // this sets the condition
  std::function<bool(const std::unique_ptr<Piece>&, const std::unique_ptr<Piece>&)> cond;

  if (std::isnan(theta0))
  {
    if (isRight)
    {
      cond = [cs, theta0](const std::unique_ptr<Piece>& q1, const std::unique_ptr<Piece>& q2)
      {
        return q1->argmax(cs) <= q2->argmax(cs);
      };
    } else
    {
      cond = [cs, theta0](const std::unique_ptr<Piece>& q1, const std::unique_ptr<Piece>& q2)
      {
        return q1->argmax(cs) >= q2->argmax(cs);
      };
    }
  } else {
    if (isRight)
    {
      cond = [cs, theta0](const std::unique_ptr<Piece>& q1, const std::unique_ptr<Piece>& q2)
      {
        return q1->argmax(cs) <= std::max(theta0, q2->argmax(cs));
      };
    } else
    {
      cond = [cs, theta0](const std::unique_ptr<Piece>& q1, const std::unique_ptr<Piece>& q2)
      {
        return q1->argmax(cs) >= std::min(theta0, q2->argmax(cs));
      };
    }

  }


  // reverse iterator at the end of the list
  // auto qi = Q.ps.rbegin();


  // compare the last two quadratics
  while (cond(Q.ps[Q.k], Q.ps[Q.k - 1])) {

    // std::cout << "RUNNING!" << std::endl;

    // prune the last quadratic (we just drop the k value of one and we overwrite the other quadratics)
    Q.k--;

    // std::cout << "PRUNED!" << std::endl;

    if (Q.k == 0){
      break;
    }

    // std::cout << "This should be fine..." << (*qi)->argmax(cs) << std::endl;
    // std::cout << "I think that here's the bug..." << (*std::next(qi, 1))->argmax(cs) << std::endl;
  }
}


/*
 * CHECK MAXIMA FUNCTIONS
 */

// get max of a single piece
//template <typename T>
double get_max (const std::unique_ptr<Piece>& q, const CUSUM& cs, const double& theta0) {
  //std::cout << q.eval(cs, argmax(q, cs), theta0) << std::endl;
  return q->eval(cs, q->argmax(cs), theta0);
}

double get_max_all (const Cost& Q, const CUSUM& cs, const double& theta0, const double& m0val) {

  double max = -INFINITY; //std::numeric_limits<double>::lowest();

  for (long unsigned int i = 0; i <= Q.k; i++)
  {
    max = std::max( max, get_max(Q.ps[i], cs, theta0) - m0val );
    // std::cout << "tau: " << Q.ps[i]->tau << " st: " << Q.ps[i]->St << " m0: " << Q.ps[i]->m0 << " max-m0val: "<< max<< " | \n";
  }
  // std::cout << std::endl;

  return max;
}

int get_tau_max (const Cost& Q, const CUSUM& cs, const double& theta0, const double& m0val)
{
  auto tau = 0;
  double max = -INFINITY;

  for (long unsigned int i = 0; i <= Q.k; i++)
  {
    if( max< get_max(Q.ps[i], cs, theta0) - m0val )
    {
      tau = Q.ps[i]->tau;
    }
  }

  return tau;
}


// to write the adaptive maxima checking


/*
 * focus recursion, one iteration
 */

void Info::update(const double& y)
{
  cs.n++;
  cs.Sn += y;

  auto m0val = 0.0;
  if (std::isnan(theta0))
    m0val = get_max(Qr.ps.front(), cs, theta0);

  prune(Qr, cs, theta0, true);
  prune(Ql, cs, theta0, false);

  Qr.opt = get_max_all(Qr, cs, theta0, m0val);
  Ql.opt = get_max_all(Ql, cs, theta0, m0val);

  if (Qr.k < (Qr.ps.size() - 1))
  {
    Qr.k++;
    Qr.ps[Qr.k]->St = cs.Sn;
    Qr.ps[Qr.k]->tau = cs.n;
    Qr.ps[Qr.k]->m0 = m0val;
  } else
  {
    Qr.ps.push_back(std::move(newP(cs.Sn, cs.n, m0val)));
    Qr.k = Qr.ps.size() - 1;
  }

  if (Ql.k < (Ql.ps.size() - 1))
  {
    Ql.k++;
    Ql.ps[Ql.k]->St = cs.Sn;
    Ql.ps[Ql.k]->tau = cs.n;
    Ql.ps[Ql.k]->m0 = m0val;
  } else
  {
    Ql.ps.push_back(std::move(newP(cs.Sn, cs.n, m0val)));
    Ql.k = Ql.ps.size() - 1;
  }
}
