package lsq

import (
	"fmt"
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

/*
func TestUpperTriMatrix2x2(t *testing.T) {

	m := NewUpperTriMatrix(2)
	m.Set(0, 0, 1)
	m.Set(0, 1, 3)
	m.Set(1, 0, 32) // should override the 3
	m.Set(1, 1, 2)
	cv.Convey("Given a 2x2 UpperTriMatrix, Set() to it should work and act symmetric, and D normalization for storage should be in effect", t, func() {
		cv.So(len(m.D), cv.ShouldResemble, 2)
		cv.So(len(m.R), cv.ShouldResemble, 1)
		cv.So(m.D, cv.ShouldResemble, []float64{1, 4})
		cv.So(m.R[0], cv.ShouldEqual, 8) // 32/4 == 8
	})
}

func TestUpperTriMatrix3x3(t *testing.T) {

	m := NewUpperTriMatrix(3)
	m.Set(0, 1, 64)
	m.Set(0, 2, 32)
	m.Set(1, 2, 16)
	m.Set(1, 1, 2)
	m.Set(2, 2, 4)

	cv.Convey("Given a 3x3 UpperTriMatrix, Set()/Get() should work, including D normalization", t, func() {
		cv.So(len(m.D), cv.ShouldResemble, 3)
		cv.So(len(m.R), cv.ShouldResemble, 3)
		cv.So(m.D, cv.ShouldResemble, []float64{1, 4, 16})
		cv.So(m.R[0], cv.ShouldEqual, 16) // 64/(2*2)
		cv.So(m.R[1], cv.ShouldEqual, 2)  // 32/(4*4)
		cv.So(m.R[2], cv.ShouldEqual, 1)  // 16/(4*4)
	})
}

func TestUpperTriMatrix4x4(t *testing.T) {

	m := NewUpperTriMatrix(4)
	m.Set(0, 1, 64)
	m.Set(0, 2, 32)
	m.Set(0, 3, 128)

	m.Set(1, 2, 16)
	m.Set(1, 3, 256)

	m.Set(2, 3, 512)

	m.Set(0, 0, 64)
	m.Set(1, 1, 2)
	m.Set(2, 2, 4)
	m.Set(3, 3, 8)

	cv.Convey("Given a 4x4 UpperTriMatrix, Set()/Get() with the D(diagonal) normalization should work", t, func() {
		cv.So(len(m.D), cv.ShouldResemble, 4)
		cv.So(len(m.R), cv.ShouldResemble, 6)
		cv.So(m.D, cv.ShouldResemble, []float64{4096, 4, 16, 64})
		cv.So(m.R, cv.ShouldResemble, []float64{16, 2, 2, 1, 4, 8})
	})
}

func TestUpperTriMatrix10x10(t *testing.T) {

	cv.Convey("Given a 10x10 UpperTriMatrix, Set() to it should work and act symmetric", t, func() {
		R := NewUpperTriMatrix(10)
		for i := 0; i < 10; i++ {
			for j := 0; j < 10; j++ {
				if j < i {
					continue
				}
				v := float64(10*(i) + j + 1)
				R.R[i*R.Ncol+j] = v
			}
		}
		cv.So(R.R, cv.ShouldResemble, []float64{2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 24, 25, 26, 27, 28, 29, 30, 35, 36, 37, 38, 39, 40, 46, 47, 48, 49, 50, 57, 58, 59, 60, 68, 69, 70, 79, 80, 90})
		cv.So(R.D, cv.ShouldResemble, []float64{1, 12, 23, 34, 45, 56, 67, 78, 89, 100})
		cv.So(R.Row_ptr, cv.ShouldResemble, []int{0, 9, 17, 24, 30, 35, 39, 42, 44, 46})
		//fmt.Printf("R.Rrow_ptr = %#v\n", R.Row_ptr)
		//fmt.Printf("R.R = %#v\n", R.R)
		//fmt.Printf("R.D = %#v\n", R.D)

		for i := 0; i < 10; i++ {
			for j := 0; j < 10; j++ {
				if j > i {
					continue
				}
				v := float64(10*(i) + j + 1)
				R.R[i, j] =v
			}
		}
		//		cv.So(R.R, cv.ShouldResemble, []float64{11, 21, 31, 41, 51, 61, 71, 81, 91, 22, 32, 42, 52, 62, 72, 82, 92, 33, 43, 53, 63, 73, 83, 93, 44, 54, 64, 74, 84, 94, 55, 65, 75, 85, 95, 66, 76, 86, 96, 77, 87, 97, 88, 98, 99})
		cv.So(R.R, cv.ShouldResemble, []float64{0.9166666666666666, 0.9130434782608695, 0.9117647058823529, 0.9111111111111111, 0.9107142857142857, 0.9104477611940298, 0.9102564102564102, 0.9101123595505618, 0.91, 0.9565217391304348, 0.9411764705882353, 0.9333333333333333, 0.9285714285714286, 0.9253731343283582, 0.9230769230769231, 0.9213483146067416, 0.92, 0.9705882352941176, 0.9555555555555556, 0.9464285714285714, 0.9402985074626866, 0.9358974358974359, 0.9325842696629213, 0.93, 0.9777777777777777, 0.9642857142857143, 0.9552238805970149, 0.9487179487179487, 0.9438202247191011, 0.94, 0.9821428571428571, 0.9701492537313433, 0.9615384615384616, 0.9550561797752809, 0.95, 0.9850746268656716, 0.9743589743589743, 0.9662921348314607, 0.96, 0.9871794871794872, 0.9775280898876404, 0.97, 0.9887640449438202, 0.98, 0.99})
		//fmt.Printf("R.R = %#v\n", R.R)
		//fmt.Printf("R.D = %#v\n", R.D)
	})

}
*/

func TestUpperCholesky2x2(t *testing.T) {

	// R for 2x2
	m := NewUpperTriMatrix(2)

	src := []float64{21, 8, 8, 38}

	m.CholeskyFactor(src)
	s := m.String()

	cv.Convey("Given a 2x2 matrix, the Cholesky decomposition should be what we expect", t, func() {
		cv.So(s, cv.ShouldEqual,
			`matrix(ncol=2, byrow=TRUE, data=c(
 4.582576e+00, 1.745743e+00, 
 0.000000e+00, 5.912054e+00))`)
	})
	fmt.Printf("\n 2x2 Cholesky factor is \n'%s'\n\n", s)
}

func TestUpperCholesky3x3(t *testing.T) {

	m := NewUpperTriMatrix(3)

	src := []float64{
		25, 15, -5,
		15, 18, 0,
		-5, 0, 11,
	}

	m.CholeskyFactor(src)
	s := m.String()

	cv.Convey("Given a 3x3 matrix, the Cholesky decomposition should be what we expect", t, func() {
		cv.So(s, cv.ShouldEqual,
			`matrix(ncol=3, byrow=TRUE, data=c(
 5.000000e+00, 3.000000e+00, -1.000000e+00, 
 0.000000e+00, 3.000000e+00, 1.000000e+00, 
 0.000000e+00, 0.000000e+00, 3.000000e+00))`)
	})

	fmt.Printf("\n 3x3 Cholesky factor is \n'%s'\n\n", s)
}

func TestUpperCholesky4x4(t *testing.T) {
	m := NewUpperTriMatrix(4)

	src := []float64{
		18, 22, 54, 42,
		22, 70, 86, 62,
		54, 86, 174, 134,
		42, 62, 134, 106,
	}

	m.CholeskyFactor(src)
	s := m.String()

	cv.Convey("Given a 4x4 matrix, the Cholesky decomposition should be what we expect", t, func() {
		cv.So(s, cv.ShouldEqual,
			`matrix(ncol=4, byrow=TRUE, data=c(
 4.242641e+00, 5.185450e+00, 1.272792e+01, 9.899495e+00, 
 0.000000e+00, 6.565905e+00, 3.046038e+00, 1.624554e+00, 
 0.000000e+00, 0.000000e+00, 1.649742e+00, 1.849711e+00, 
 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.392621e+00))`)
	})

	fmt.Printf("\n 4x4 Cholesky factor is \n'%s'\n\n", s)

}

func TestAt(t *testing.T) {

	nc := 3
	m := NewUpperTriMatrix(3)

	src := []float64{
		25, 15, -5,
		15, 18, 0,
		-5, 0, 11,
	}

	m.CholeskyFactor(src)
	s := m.String()

	cv.Convey("Given a 3x3 matrix, the Cholesky decomposition should be what we expect", t, func() {
		cv.So(s, cv.ShouldEqual,
			`matrix(ncol=3, byrow=TRUE, data=c(
 5.000000e+00, 3.000000e+00, -1.000000e+00, 
 0.000000e+00, 3.000000e+00, 1.000000e+00, 
 0.000000e+00, 0.000000e+00, 3.000000e+00))`)
	})

	expected := []float64{5, 3, -1, 0, 3, 1, 0, 0, 3}

	pos := 0
	cv.Convey("Given a 3x3 matrix in UpperTriMatrix, At(i,j) should match what we expect", t, func() {
		for i := 0; i < nc; i++ {
			for j := 0; j < nc; j++ {
				//fmt.Printf("\nm.At(%d, %d) = %v\n", i, j, m.At(i, j))
				cv.So(m.At(i, j), cv.ShouldEqual, expected[pos])
				pos++
			}
		}
	})
}
