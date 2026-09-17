#ifndef SVPCHANGE_AR1_FOCUS_STATE_H
#define SVPCHANGE_AR1_FOCUS_STATE_H

namespace changepoint {

void update_ar1_focus_state(double observation,
                            double rho,
                            void*& opaque_state,
                            double& max_statistic);

void destroy_ar1_focus_state(void* opaque_state);

// Owns the opaque state of the fixed-rho AR(1) FOCuS recurrence.
class AR1FocusState
{
public:
  explicit AR1FocusState(double rho)
    : rho_(rho),
      max_statistic_(-1.0),
      opaque_state_(nullptr) {}

  ~AR1FocusState()
  {
    destroy_ar1_focus_state(opaque_state_);
  }

  AR1FocusState(const AR1FocusState&) = delete;
  AR1FocusState& operator=(const AR1FocusState&) = delete;

  AR1FocusState(AR1FocusState&& other) noexcept
    : rho_(other.rho_),
      max_statistic_(other.max_statistic_),
      opaque_state_(other.opaque_state_)
  {
    other.opaque_state_ = nullptr;
  }

  AR1FocusState& operator=(AR1FocusState&& other) noexcept
  {
    if (this == &other) return *this;
    destroy_ar1_focus_state(opaque_state_);
    rho_ = other.rho_;
    max_statistic_ = other.max_statistic_;
    opaque_state_ = other.opaque_state_;
    other.opaque_state_ = nullptr;
    return *this;
  }

  void update(double observation)
  {
    update_ar1_focus_state(
      observation, rho_, opaque_state_, max_statistic_
    );
  }

  double max_statistic() const
  {
    return max_statistic_;
  }

private:
  double rho_;
  double max_statistic_;
  void* opaque_state_;
};

} // namespace changepoint

#endif
