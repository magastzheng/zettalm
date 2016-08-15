package lsq

import (
	"github.com/glycerine/gostat"
	"github.com/glycerine/gostat/fn"
	"math"
)

// return  P[ abs(x) > T], the right-hand tail probability, where
// T ~ t_{n}  (Student's t-distribution with n degrees of freedom).
//
// Note: this is the CENTRAL t distribution. Contast it with the NON-central t distribution.
//
// If x is negative, we'll convert it to abs(x) (make it positive) before providing the
// upper (right-hand) tail probability corresponding to the (now) positive point.
//
// Double the returned value to get a 2-sided test p-value.
//
func Pt(x float64, n float64) float64 {
	if math.IsNaN(x) || math.IsNaN(n) {
		return math.NaN()
	}

	if x == 0 {
		return 0.5 // distribution is symmetric around 0
	}

	if math.IsInf(x, 0) {
		return 0 // either +Inf or -Inf
	}

	// make x always positive
	if x < 0 {
		x = -x
	}

	var val, nx float64
	if n <= 0.0 {
		panic("degrees of freedom for the student's t-distribution cannot be <= 0")
	}

	if math.IsInf(n, 0) {
		// normal is now exact...
		return PvalueAssumingNormalDistribution(x)
	}

	var lval float64
	nx = 1 + (x/n)*x
	if nx > 1e100 { // <=>  x*x > 1e100 * n
		// really only matters in log-space. Keep around in case we
		// want to compute log of tail probability later.
		//
		// To address the danger of underflow, use Abramowitz & Stegun 26.5.4
		//  pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
		//   with z = 1/nx,  a = n/2,  b= 1/2 :
		//
		lval = -0.5*n*(2*math.Log(math.Abs(x))-math.Log(n)) - fn.LnBeta(0.5*n, 0.5) - math.Log(0.5*n)
		val = math.Exp(lval)
	} else {
		if n > x*x {
			betaCDF := gostat.Beta_CDF(0.5, n/2)
			val = 0.5 - betaCDF(x*x/(n+x*x)) + 0.5
		} else {
			betaCDF := gostat.Beta_CDF(n/2, 0.5)
			val = betaCDF(1 / nx)
		}
	}
	val /= 2 // default to 1-sided p-value. Double it for 2-sided.
	return val
}
