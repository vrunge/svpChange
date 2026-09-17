#ifndef SVPCHANGE_VALIDITY_TESTS_H
#define SVPCHANGE_VALIDITY_TESTS_H

// Compatibility umbrella for the built-in incremental validity tests.
// Keeping this header lets external simulation code continue to include
// validity_tests.h while the implementations remain organized by test family.
#include "svp_test_base.h"
#include "svp_tests_gamma.h"
#include "svp_tests_gaussian.h"
#include "svp_tests_ar1.h"
#include "svp_tests_quantile.h"
#include "svp_tests_rank.h"

// Compatibility aliases for code that included the former monolithic header.
using AR1MeanChange = AR1FocusMeanChange;
using varCost = VarianceCost;

#endif
