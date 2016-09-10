package lsq

import (
	"fmt"
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

func TestDotProductOnRD(t *testing.T) {

	// 50 rows
	// 9 columns, including row count in first column
	big, err := readData("bigger.dat")
	if err != nil {
		panic(err)
	}

	nxvar := 9
	nyvar := 0
	m := NewMillerLSQ(nxvar, nyvar)
	for i := range big.Rows {
		m.Includ(1.0, big.Rows[i], []float64{}, NAN_TO_ZERO)
	}
	square := m.ExtractSquareR(true)
	strSquare := square.String()
	fmt.Printf("square = %v\n", strSquare)

	var knowngood_ij, myRdot_ij float64

	cv.Convey("Given a QR-decomposition in LSQ from bigger.dat", t, func() {
		cv.Convey("The choleskry factorization (R&D) should match our expected (positive control) QR from R", func() {
			cv.So(strSquare, cv.ShouldEqual, `matrix(ncol=10, nrow=10, byrow=TRUE, data=c(
 7.071068, 180.312229, 34.793767, 46.643465, 76.305057, 34.961242, 93.183337, 758.260075, 2468.285563, -1707.858316, 
 0.000000, 102.041658, 10.027417, -2.764113, 13.237656, 84.967077, 74.897479, 2.682443, -123.836867, 128.270014, 
 0.000000, 0.000000, 26.198489, -1.030243, 77.458759, 9.170345, -7.728494, 7.792060, -2042.066148, 2049.417947, 
 0.000000, 0.000000, 0.000000, 20.316808, -21.153666, -32.256000, 32.703317, 19.556507, 4025.741664, -4006.441254, 
 0.000000, 0.000000, 0.000000, 0.000000, 188.494989, 86.021279, -12.213250, 18.224347, -2325.028732, 2345.378941, 
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 578.528125, 16.078763, 3.430751, 274.795346, -271.424121, 
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 248.025968, 4.222613, -162.430327, 168.737549, 
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 71.686978, 1719.991349, -1650.110708, 
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 13218.359528, -13219.227587, 
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 13.847606))`)
		})
		cv.Convey("and we should compute the correct dot product of the i,j columns for the implied matrix represented by R and D", func() {
			for i := 0; i < nxvar; i++ {
				for j := 0; j < nxvar; j++ {
					knowngood_ij = square.ColumnDotProduct(i, j)
					myRdot_ij = m.DotColumns(i, j)
					cv.So(myRdot_ij, cv.ShouldEqual, knowngood_ij)
				}
			}
		})
	})
}
