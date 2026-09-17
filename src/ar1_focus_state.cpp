#include "ar1_focus_state.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace changepoint {
namespace {

struct Triple {
  int tau = 0;
  double m0 = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double l = 0.0;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double E = 0.0;
  double f = 0.0;
};

struct State {
  State() { triples.reserve(32); }

  std::vector<Triple> triples;
  double max_val = -1.0;
};

struct Ar1State {
  explicit Ar1State(double rho_) : rho(rho_), one_minus_rho(1.0 - rho_) {}

  double rho;
  double one_minus_rho;
  int observation_count = 0;
  int innovation_count = 0;
  double previous_observation = 0.0;
  double innovation_sum = 0.0;
  double innovation_sum_square = 0.0;
  std::array<double, 3> innovations{{0.0, 0.0, 0.0}};
  int buffered = 0;
  State right;
  State left;
  State right_negative;
  State left_negative;
};

inline void push_innovation(Ar1State& state, double value) {
  if (state.buffered < 3) {
    state.innovations[static_cast<std::size_t>(state.buffered++)] = value;
  } else {
    state.innovations[0] = state.innovations[1];
    state.innovations[1] = state.innovations[2];
    state.innovations[2] = value;
  }
  state.innovation_sum += value;
  state.innovation_sum_square += value * value;
  ++state.innovation_count;
}

inline Triple introduce_triple(int tau, double m0, double m1, double m2,
                               double intersection, double rho,
                               double one_minus_rho, double last_innovation,
                               double sum_square) {
  Triple out;
  out.tau = tau;
  out.m0 = m0;
  out.m1 = m1;
  out.m2 = m2;
  out.l = intersection;
  out.A = -0.5 * (static_cast<double>(tau) * one_minus_rho * one_minus_rho +
                  rho * rho);
  out.B = -0.5;
  out.C = rho;
  out.D = one_minus_rho * m1 - rho * last_innovation;
  out.E = last_innovation;
  out.f = -0.5 * sum_square;
  return out;
}

inline double intersection(const Triple& triple, int current_n,
                           const Ar1State& detector, double sign) {
  const double last = sign * detector.innovations[2];
  const double boundary_delta = (triple.m2 - triple.m1) - last;
  const double tail_sum = sign * detector.innovation_sum - triple.m2;
  const double numerator = 2.0 *
    (boundary_delta + detector.one_minus_rho * tail_sum);
  const double denominator = detector.one_minus_rho * detector.one_minus_rho *
    static_cast<double>((current_n - 1) - triple.tau);
  if (denominator == 0.0) return 0.0;
  return std::max(0.0, numerator / denominator);
}

inline void prune(State& state, int current_n, const Ar1State& detector,
                  bool right_side, double sign) {
  if (state.triples.empty()) return;

  std::size_t active = state.triples.size();
  state.triples[active - 1].l =
    intersection(state.triples[active - 1], current_n, detector, sign);
  const double boundary = right_side
    ? -std::numeric_limits<double>::infinity()
    : std::numeric_limits<double>::infinity();

  while (active > 0) {
    const double current = state.triples[active - 1].l;
    const double previous = active == 1 ? boundary : state.triples[active - 2].l;
    const bool remove = right_side ? current <= previous : current >= previous;
    if (!remove) break;
    --active;
    if (active > 0) {
      state.triples[active - 1].l =
        intersection(state.triples[active - 1], current_n, detector, sign);
    }
  }
  state.triples.resize(active);
}

inline double optimized_loglikelihood(const Triple& triple) {
  const double denominator = 4.0 * triple.A * triple.B - triple.C * triple.C;
  if (std::abs(denominator) < 1e-12) {
    return std::numeric_limits<double>::lowest();
  }
  const double inverse = 1.0 / denominator;
  const double mu0 = -inverse *
    (2.0 * triple.B * triple.D - triple.C * triple.E);
  const double mu1 = -inverse *
    (-triple.C * triple.D + 2.0 * triple.A * triple.E);
  return triple.A * mu0 * mu0 + triple.B * mu1 * mu1 +
    triple.C * mu0 * mu1 + triple.D * mu0 + triple.E * mu1 + triple.f;
}

inline void compute_max(State& state, int current_n,
                        const Ar1State& detector, double sign) {
  state.max_val = -1.0;

  const double a0 = -0.5 * static_cast<double>(current_n) *
    detector.one_minus_rho * detector.one_minus_rho;
  const double b0 = detector.one_minus_rho * sign * detector.innovation_sum;
  const double c0 = -0.5 * detector.innovation_sum_square;
  const double null_mean = -(b0 / (2.0 * a0));
  const double null_loglikelihood =
    a0 * null_mean * null_mean + b0 * null_mean + c0;

  for (const Triple& triple : state.triples) {
    // The recurrence tau is one less than the observation split. The exact
    // AR1 validity test leaves at least two observations on each side.
    if (triple.tau < 1 || triple.tau > current_n - 2) continue;
    const double value = 2.0 *
      (optimized_loglikelihood(triple) - null_loglikelihood);
    if (value > state.max_val) {
      state.max_val = value;
    }
  }
}

inline void update_side(State& state, Ar1State& detector, bool right_side,
                        double sign) {
  const int current_n = detector.innovation_count;
  const double last = sign * detector.innovations[2];

  if (current_n == 3) {
    const double first = sign * detector.innovations[0];
    const double second = sign * detector.innovations[1];
    state.triples.clear();
    state.triples.push_back(introduce_triple(
      1, 0.0, first, first + second, 0.0, detector.rho,
      detector.one_minus_rho, second,
      detector.innovation_sum_square - last * last
    ));
  }

  prune(state, current_n, detector, right_side, sign);

  const double signed_sum = sign * detector.innovation_sum;
  const double m0 = signed_sum -
    sign * detector.innovations[1] - sign * detector.innovations[2];
  const double m1 = signed_sum - sign * detector.innovations[2];
  const double m2 = signed_sum;
  Triple newest = introduce_triple(
    current_n - 1, m0, m1, m2, 2.0 * (m2 - m1), detector.rho,
    detector.one_minus_rho, last, detector.innovation_sum_square
  );
  state.triples.push_back(newest);

  const double half_weight =
    0.5 * detector.one_minus_rho * detector.one_minus_rho;
  for (std::size_t i = 0; i + 1 < state.triples.size(); ++i) {
    state.triples[i].B -= half_weight;
    state.triples[i].E += detector.one_minus_rho * last;
    state.triples[i].f -= 0.5 * last * last;
  }

  compute_max(state, current_n, detector, sign);
}

}  // namespace

void update_ar1_focus_state(double observation,
                            double rho,
                            void*& opaque_state,
                            double& out_max_stat) {
  Ar1State* state = reinterpret_cast<Ar1State*>(opaque_state);
  if (state == nullptr) {
    state = new Ar1State(rho);
    opaque_state = state;
  }

  ++state->observation_count;
  if (state->observation_count == 1) {
    state->previous_observation = observation;
    out_max_stat = -1.0;
    return;
  }

  const double innovation = observation - state->rho * state->previous_observation;
  state->previous_observation = observation;
  push_innovation(*state, innovation);

  if (state->innovation_count < 3) {
    out_max_stat = -1.0;
    return;
  }

  update_side(state->right, *state, true, 1.0);
  update_side(state->left, *state, false, 1.0);
  update_side(state->right_negative, *state, true, -1.0);
  update_side(state->left_negative, *state, false, -1.0);

  const State* best = &state->right;
  const State* candidates[] = {
    &state->left, &state->right_negative, &state->left_negative
  };
  for (const State* candidate : candidates) {
    if (candidate->max_val > best->max_val) best = candidate;
  }
  out_max_stat = best->max_val;
}

void destroy_ar1_focus_state(void* opaque_state) {
  delete reinterpret_cast<Ar1State*>(opaque_state);
}

}  // namespace changepoint
