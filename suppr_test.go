package lsq

import (
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

func TestSupplementDeSupplement(t *testing.T) {

	ncDesup := 3
	ncSup := ncDesup + 1
	smaller := NewUpperTriMatrix(ncDesup)
	smaller2 := NewUpperTriMatrix(ncDesup)
	larger := NewUpperTriMatrix(ncSup)

	for i := 0; i < ncDesup; i++ {
		for j := i; j < ncDesup; j++ {
			smaller.Set(i, j, float64(i*100+j))
		}
	}

	rhs := make([]float64, ncDesup)
	for i := 0; i < ncDesup; i++ {
		rhs[i] = float64(1000 + i)
	}
	sserr := float64(99)

	sup := SupplementR(smaller, rhs, sserr, NoRecycle)
	desup2, rhs2, sserr2 := sup.DeSupplementR(NoRecycle)

	cv.Convey("Given an R, Rhs, Sserr, SupplementR() should combine them into a R.Ncol+1 x R.Ncol+1 UpperTriMatrix, and DeSupplementR() should be the inverse.", t, func() {
		cv.So(len(desup2.R), cv.ShouldEqual, len(smaller.R))
		cv.So(len(desup2.D), cv.ShouldEqual, len(smaller.D))
		cv.So(desup2.Ncol, cv.ShouldEqual, smaller.Ncol)
		cv.So(desup2.R_dim, cv.ShouldEqual, smaller.R_dim)
		cv.So(desup2.Row_ptr, cv.ShouldResemble, smaller.Row_ptr)

		cv.So(sserr2, cv.ShouldEqual, sserr)

		for i := 0; i < ncDesup; i++ {
			for j := i; j < ncDesup; j++ {
				cv.So(desup2.At(i, j), cv.ShouldEqual, smaller.At(i, j))
			}
		}

		cv.So(rhs2, cv.ShouldResemble, rhs)
	})

	// recylce now
	sup = SupplementR(smaller, rhs, sserr, larger)
	desup2, rhs2, sserr2 = sup.DeSupplementR(smaller2)

	cv.Convey("Given an R, Rhs, Sserr, SupplementR() should combine them into a R.Ncol+1 x R.Ncol+1 UpperTriMatrix, and DeSupplementR() should be the inverse--with recycling this time", t, func() {

		cv.Convey("We should have recycled the larger matrix for sup", func() {
			cv.So(sup, cv.ShouldEqual, larger)
		})

		cv.Convey("We should have recycled the smaller2 matrix for desup2", func() {
			cv.So(desup2, cv.ShouldEqual, smaller2)
		})

		cv.So(len(desup2.R), cv.ShouldEqual, len(smaller.R))
		cv.So(len(desup2.D), cv.ShouldEqual, len(smaller.D))
		cv.So(desup2.Ncol, cv.ShouldEqual, smaller.Ncol)
		cv.So(desup2.R_dim, cv.ShouldEqual, smaller.R_dim)
		cv.So(desup2.Row_ptr, cv.ShouldResemble, smaller.Row_ptr)

		cv.So(sserr2, cv.ShouldEqual, sserr)

		for i := 0; i < ncDesup; i++ {
			for j := i; j < ncDesup; j++ {
				cv.So(desup2.At(i, j), cv.ShouldEqual, smaller.At(i, j))
			}
		}

		cv.So(rhs2, cv.ShouldResemble, rhs)
	})

}
