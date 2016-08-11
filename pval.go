package lsq

import "math"

// PvalueNormalDist(t):
//
// To go from t-value to 2-sided p-value,
// we'll use the normal approximation to Student's-T distribution, meaning
// that it will be slightly different (in the typically insignificant 3rd decimal place) from
// what R computes. R uses the actual T-distribution, and goes to great, great lengths
// to be exact -- See for instance R-3.1.1/src/nmath/pt.c for the gorey details.
// The normal approximation to the T-distribution is less exact for N < 30,
// but still the differences are unlikely to be significant.
//
// We are assuming a 2-sided p-value to be conservative. 1-sided p-values
// would be half of the returned p-value.
//
// Reference: https://en.wikipedia.org/wiki/Error_function
//
func PvalueAssumingNormalDistribution(t float64) float64 {
	return math.Erfc(math.Abs(t) / math.Sqrt2)
}
