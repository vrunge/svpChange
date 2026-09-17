#ifndef SVPCHANGE_SVP_TESTS_QUANTILE_H
#define SVPCHANGE_SVP_TESTS_QUANTILE_H

#include "svp_test_base.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

// ----------------------------------------------------------
// Order-statistics tree (treap) for doubles, with multiplicities
// ----------------------------------------------------------

class OrderStatisticTree
{
private:
  struct Node {
    double key;
    int priority;   // heap key
    int count;      // multiplicity of this key
    int size;       // total size of subtree (with multiplicities)
    Node* left;
    Node* right;

    Node(double k, int pr)
      : key(k), priority(pr), count(1), size(1),
        left(nullptr), right(nullptr) {}
  };

  Node* root_;
  std::uint64_t rng_state_;

  static int get_size(Node* n) {
    return n ? n->size : 0;
  }

  static void recalc(Node* n) {
    if (!n) return;
    n->size = n->count + get_size(n->left) + get_size(n->right);
  }

  Node* rotate_right(Node* y) {
    Node* x = y->left;
    Node* T2 = x->right;
    x->right = y;
    y->left = T2;
    recalc(y);
    recalc(x);
    return x;
  }

  Node* rotate_left(Node* x) {
    Node* y = x->right;
    Node* T2 = y->left;
    y->left = x;
    x->right = T2;
    recalc(x);
    recalc(y);
    return y;
  }

  int random_priority() {
    // xorshift64*: fast deterministic priorities without constructing a
    // hardware-seeded random engine for every SVP candidate.
    rng_state_ ^= rng_state_ >> 12;
    rng_state_ ^= rng_state_ << 25;
    rng_state_ ^= rng_state_ >> 27;
    return static_cast<int>((rng_state_ * 2685821657736338717ULL) >> 33);
  }

  Node* insert(Node* node, double key) {
    if (!node)
      return new Node(key, random_priority());

    if (key == node->key) {
      node->count += 1;
    } else if (key < node->key) {
      node->left = insert(node->left, key);
      if (node->left->priority > node->priority)
        node = rotate_right(node);
    } else {
      node->right = insert(node->right, key);
      if (node->right->priority > node->priority)
        node = rotate_left(node);
    }
    recalc(node);
    return node;
  }

  double kth(Node* node, int k) const {
    // k is 0-based, assume 0 <= k < size(root)
    int left_size = get_size(node->left);
    if (k < left_size) {
      return kth(node->left, k);
    }
    if (k < left_size + node->count) {
      return node->key;
    }
    return kth(node->right, k - left_size - node->count);
  }

  void clear(Node* node) {
    if (!node) return;
    clear(node->left);
    clear(node->right);
    delete node;
  }


public:
  OrderStatisticTree()
    : root_(nullptr),
      rng_state_(0x9e3779b97f4a7c15ULL) {}

  OrderStatisticTree(const OrderStatisticTree&) = delete;
  OrderStatisticTree& operator=(const OrderStatisticTree&) = delete;

  OrderStatisticTree(OrderStatisticTree&& other) noexcept
    : root_(other.root_), rng_state_(other.rng_state_) {
    other.root_ = nullptr;
  }

  OrderStatisticTree& operator=(OrderStatisticTree&& other) noexcept {
    if (this != &other) {
      clear(root_);
      root_ = other.root_;
      rng_state_ = other.rng_state_;
      other.root_ = nullptr;
    }
    return *this;
  }

  ~OrderStatisticTree() {
    clear(root_);
  }

  void insert(double key) {
    root_ = insert(root_, key);
  }


  int count_less(double key) const {
    int result = 0;
    Node* node = root_;
    while (node) {
      if (key <= node->key) {
        node = node->left;
      } else {
        result += get_size(node->left) + node->count;
        node = node->right;
      }
    }
    return result;
  }

  int count_equal(double key) const {
    Node* node = root_;
    while (node) {
      if (key < node->key) node = node->left;
      else if (key > node->key) node = node->right;
      else return node->count;
    }
    return 0;
  }

  double kth(int k) const {
    // No range check here; caller must ensure 0 <= k < size()
    return kth(root_, k);
  }
};

////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------
// Exact quantile-based cost: C^{quant}(y_{a..b}; x) = q_{1-x} - q_x
// with O(log n) update via order-statistics tree
// ----------------------------------------------------------

class QuantileCostExact : public TestBase
{
public:
  // x in (0, 0.5] typically
  QuantileCostExact(double x)
    : x_(x),
      tree_(),
      n_(0),
      q_low_(std::numeric_limits<double>::quiet_NaN()),
      q_high_(std::numeric_limits<double>::quiet_NaN())
  {}

  void update(double y) override
  {
    tree_.insert(y);
    ++n_;

    if (n_ == 1) {
      // Single point: both quantiles are that point
      q_low_ = q_high_ = y;
      return;
    }

    // 0-based indices in [0, n_-1]
    int idx_low  = static_cast<int>(std::floor(x_ * (n_ - 1)));
    int idx_high = static_cast<int>(std::floor((1.0 - x_) * (n_ - 1)));

    if (idx_low < 0) idx_low = 0;
    if (idx_high < 0) idx_high = 0;
    if (idx_low >= n_) idx_low = n_ - 1;
    if (idx_high >= n_) idx_high = n_ - 1;

    q_low_  = tree_.kth(idx_low);
    q_high_ = tree_.kth(idx_high);
  }

  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();

    return q_high_ - q_low_;
  }

private:
  double x_;
  OrderStatisticTree tree_;
  int n_;
  double q_low_;
  double q_high_;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////





// ----------------------------------------------------------
// Streaming quantile estimator (P² algorithm) for one p
// ----------------------------------------------------------

class P2Quantile
{
public:
  explicit P2Quantile(double p) : p_(p)
  {
    reset();
  }

  void reset()
  {
    count_ = 0;
  }

  void add(double x)
  {
    // Collect first 5 points
    if (count_ < 5) {
      q_[count_] = x;
      ++count_;
      if (count_ == 5) {
        std::sort(q_, q_ + 5);

        for (int i = 0; i < 5; ++i)
          n_[i] = i;

        ns_[0]  = 0.0;
        ns_[1]  = 2.0 * p_;
        ns_[2]  = 4.0 * p_;
        ns_[3]  = 2.0 + 2.0 * p_;
        ns_[4]  = 4.0;
        dns_[0] = 0.0;
        dns_[1] = p_ / 2.0;
        dns_[2] = p_;
        dns_[3] = (1.0 + p_) / 2.0;
        dns_[4] = 1.0;
      }
      return;
    }

    // General case: update markers
    int k;
    if (x < q_[0]) {
      q_[0] = x;
      k = 0;
    } else if (x < q_[1]) {
      k = 0;
    } else if (x < q_[2]) {
      k = 1;
    } else if (x < q_[3]) {
      k = 2;
    } else if (x < q_[4]) {
      k = 3;
    } else {
      q_[4] = x;
      k = 3;
    }

    // Update positions
    for (int i = k + 1; i < 5; ++i)
      ++n_[i];
    for (int i = 0; i < 5; ++i)
      ns_[i] += dns_[i];

    // Adjust interior markers
    for (int i = 1; i <= 3; ++i) {
      double d = ns_[i] - n_[i];
      if ((d >= 1.0 && n_[i + 1] - n_[i] > 1) ||
          (d <= -1.0 && n_[i - 1] - n_[i] < -1)) {
        int dInt = (d > 0.0) ? 1 : -1;
        double qs = parabolic(i, dInt);
        if (q_[i - 1] < qs && qs < q_[i + 1])
          q_[i] = qs;
        else
          q_[i] = linear(i, dInt);
        n_[i] += dInt;
      }
    }

    ++count_;
  }

  double quantile() const
  {
    if (count_ == 0)
      return std::numeric_limits<double>::quiet_NaN();

    // For the first few points, just use the empirical quantile
    if (count_ <= 5) {
      double tmp[5];
      for (int i = 0; i < count_; ++i)
        tmp[i] = q_[i];
      std::sort(tmp, tmp + count_);

      double pos  = (count_ - 1) * p_;
      int index   = static_cast<int>(std::round(pos));
      if (index < 0) index = 0;
      if (index >= count_) index = count_ - 1;
      return tmp[index];
    }

    // Steady state: p-quantile is q_[2]
    return q_[2];
  }

private:
  double parabolic(int i, int d) const
  {
    return q_[i] + d / static_cast<double>(n_[i + 1] - n_[i - 1]) *
      ( (n_[i] - n_[i - 1] + d) * (q_[i + 1] - q_[i]) / static_cast<double>(n_[i + 1] - n_[i]) +
      (n_[i + 1] - n_[i] - d) * (q_[i] - q_[i - 1]) / static_cast<double>(n_[i] - n_[i - 1]) );
  }

  double linear(int i, int d) const
  {
    return q_[i] + d * (q_[i + d] - q_[i]) / static_cast<double>(n_[i + d] - n_[i]);
  }

  double p_;       // target quantile in (0,1)
  int    n_[5];    // marker positions
  double ns_[5];   // desired positions
  double dns_[5];  // position increments
  double q_[5];    // marker heights
  int    count_{}; // number of observations seen
};

// ----------------------------------------------------------
// Quantile-based cost: C^{quant}(y_{a..b}; x) = q_{1-x} - q_x
// ----------------------------------------------------------

class QuantileCost : public TestBase
{
public:
  // x is the lower quantile level in (0, 0.5], typically
  explicit QuantileCost(double x)
    : lower_(x),
      upper_(1.0 - x),
      q_low_(std::numeric_limits<double>::quiet_NaN()),
      q_high_(std::numeric_limits<double>::quiet_NaN()),
      n_(0)
  {}

  void update(double y) override
  {
    lower_.add(y);
    upper_.add(y);
    ++n_;

    // Cache current quantiles so statistic() is trivial
    q_low_  = lower_.quantile();
    q_high_ = upper_.quantile();
  }

  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();

    // O(1): just use the cached values
    return q_high_ - q_low_;
  }

private:
  P2Quantile lower_; // estimator for q_x
  P2Quantile upper_; // estimator for q_{1-x}
  double q_low_;
  double q_high_;
  int    n_;        // number of points seen
};




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
