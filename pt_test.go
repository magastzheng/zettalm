package lsq

import (
	"fmt"
	"math"
	"testing"

	cv "github.com/smartystreets/goconvey/convey"
)

func TestNormalPvalue(t *testing.T) {
	const eps = 1e-5
	cv.Convey("Given the PvalueAssumingNormalDistribution() function, it should return Normal distribution tail probabilities (for the upper-tail)", t, func() {
		cv.So(WithinEpsilon(PvalueAssumingNormalDistribution(1.96), .05, eps), cv.ShouldBeTrue)
		cv.So(WithinEpsilon(PvalueAssumingNormalDistribution(2.0), 0.0455, eps), cv.ShouldBeTrue)
	})
}

func TestPt(t *testing.T) {
	const eps = 1e-5
	cv.Convey("Given the Pt function, it should return reasonable Student's T distrubution tail probabilities", t, func() {
		//
		// from R:
		//
		//	> tval=c(1.96, 7.49)
		//	> rdf= 45
		//	> pt(abs(tval), rdf, lower.tail = FALSE)
		// [1] 2.810293e-02 9.616653e-10
		//	> rdf=100
		//	> pt(abs(tval), rdf, lower.tail = FALSE)
		// [1] 2.638945e-02 1.395068e-11

		// n > x*x path
		cv.So(WithinEpsilon(Pt(1.96, 45), 2.810293e-02, eps), cv.ShouldBeTrue)
		cv.So(WithinEpsilon(Pt(7.49, 45), 9.616653e-10, eps), cv.ShouldBeTrue)
		cv.So(WithinEpsilon(Pt(1.96, 100), 2.638945e-02, eps), cv.ShouldBeTrue)
		cv.So(WithinEpsilon(Pt(7.49, 100), 1.395068e-11, eps), cv.ShouldBeTrue)

		// n < x*x path
		// n = 5; x=3 => x*x = 9
		//
		// from R:
		//> pt(3, 5, ncp=F,lower.tail=F)
		//[1] 0.01504962
		//> pt(8, 5, ncp=F,lower.tail=F)
		//[1] 0.0002464533
		//>
		cv.So(WithinEpsilon(Pt(3, 5), 0.01504962, eps), cv.ShouldBeTrue)
		cv.So(WithinEpsilon(Pt(8, 5), 0.0002464533, eps), cv.ShouldBeTrue)

		// nx > 1e100
		cv.So(WithinEpsilon(Pt(3e55, 6), 0, eps), cv.ShouldBeTrue)
	})
}

func WithinEpsilon(x float64, y float64, epsilon float64) bool {
	if math.Abs(x-y) <= epsilon {
		return true
	}
	fmt.Printf("\nWithinEpsilon failing! difference: abs(%v - %v) = %v\n", x, y, math.Abs(x-y))
	return false
}
