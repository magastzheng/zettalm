package lsq

import (
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

func TestDframeAddFirstOnesColumn(t *testing.T) {
	x, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	cv.Convey("When smallfuel.dat is loaded into X, it should read what we expect", t, func() {
		cv.So(x.String(), cv.ShouldEqual, `DataFrame [10 x 4] = 
 RoadMls   FuelCon   DLic      Fuel_Pop 
 1.976  557.000  52.500  541.000 
 1.250  404.000  57.200  524.000 
 1.586  259.000  58.000  561.000 
 2.351  2396.000  52.900  414.000 
 0.431  397.000  54.400  410.000 
 1.333  1408.000  57.100  457.000 
 11.868  6312.000  45.100  344.000 
 2.138  3439.000  55.300  467.000 
 8.577  5528.000  52.900  464.000 
 8.507  5375.000  55.200  498.000 
`)
	})

	x.AddFirstOnesColumn()
	cv.Convey("When smallfuel.dat is loaded into X, *and* we call AddFirstOnesColumnit should read what we expect", t, func() {
		cv.So(x.String(), cv.ShouldEqual, `DataFrame [10 x 5] = 
 const  RoadMls   FuelCon   DLic      Fuel_Pop 
 1.000  1.976  557.000  52.500  541.000 
 1.000  1.250  404.000  57.200  524.000 
 1.000  1.586  259.000  58.000  561.000 
 1.000  2.351  2396.000  52.900  414.000 
 1.000  0.431  397.000  54.400  410.000 
 1.000  1.333  1408.000  57.100  457.000 
 1.000  11.868  6312.000  45.100  344.000 
 1.000  2.138  3439.000  55.300  467.000 
 1.000  8.577  5528.000  52.900  464.000 
 1.000  8.507  5375.000  55.200  498.000 
`)
	})
}
