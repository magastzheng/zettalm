package lsq

import (
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

// on actual data
func TestSerializationOnFuelData(t *testing.T) {

	nvar := 7
	preDisk := NewMillerLSQ(nvar, 1)

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("in fuelconst.dat, df.Rows[0] = %v\n", df.Rows[0])
	last := df.Ncol - 1
	for i := range df.Rows {
		preDisk.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:], NAN_TO_ZERO)
	}

	// now save and restore
	fn := "fuel_fit.gob"
	preDisk.GobEncode(fn)
	//preDisk.LsqToCapnpFile(fn, nil, nil, nil, nil)

	m := ReadLSQGobFile(fn)
	//m, _, _, err := LsqFromCapnpFile(fn)
	if err != nil {
		panic(err)
	}

	mat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)

	emat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormatE)
	//fmt.Printf("R is \n%v\n", mat)

	cv.Convey("Given an LSQ restored from disk (with the fuelcons.dat data loaded, after 48 Includ() calls but before any other calls)", t, func() {

		cv.Convey("Then the previously seen lsq.Rows should be 48", func() {
			cv.So(m.Nobs, cv.ShouldEqual, 48)
		})

		cv.Convey("Then the R matrix should match the lsq.f90 version output", func() {
			cv.So(mat, cv.ShouldEqual, `
[ 1.000000 4296.916667 7.668333 2362.083333 4.241833 5.565417 2253.041667 57.033333 
 0.000000 1.000000 -0.000031 0.532401 0.000053 0.000508 0.463974 -0.000456 
 0.000000 0.000000 1.000000 -87.904537 0.044881 -1.603833 -215.545385 -2.036955 
 0.000000 0.000000 0.000000 1.000000 -0.000000 -0.000707 1.270032 0.007452 
 0.000000 0.000000 0.000000 0.000000 1.000000 -1.347978 -260.251651 3.888110 
 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 54.974483 0.547643 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 -0.006200 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000]`)

			// check more precisely for those small entries like row 4, column 5
			cv.So(emat, cv.ShouldEqual, `
[ 1.000000e+00 4.296917e+03 7.668333e+00 2.362083e+03 4.241833e+00 5.565417e+00 2.253042e+03 5.703333e+01 
 0.000000e+00 1.000000e+00 -3.139515e-05 5.324008e-01 5.295497e-05 5.081340e-04 4.639739e-01 -4.564454e-04 
 0.000000e+00 0.000000e+00 1.000000e+00 -8.790454e+01 4.488092e-02 -1.603833e+00 -2.155454e+02 -2.036955e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -4.463739e-09 -7.071586e-04 1.270032e+00 7.451966e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -1.347978e+00 -2.602517e+02 3.888110e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 5.497448e+01 5.476433e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -6.199548e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00]`)
		})

		u := m.XStats.Mean()
		sd := m.XStats.Sd()

		//fmt.Printf("u = %v\n", u)
		//fmt.Printf("sd = %v\n", sd)

		kgMean := []float64{4296.916666666666, 7.668333333333333, 2362.0833333333335, 4.241833333333332, 5.565416666666667, 2253.0416666666665, 57.033333333333324}
		kgSd := []float64{4441.10870206372, 0.9507697516051801, 2384.9691182561382, 0.5736237677697418, 3.491507166078876, 2124.270577772204, 5.547026549972453}
		cv.Convey("The mean and stddev should be what we expect", func() {
			//fmt.Printf("u = %v\n", u)
			//fmt.Printf("kgMean = %v\n", kgMean)
			cv.So(EpsSliceEqual(u, kgMean, 1e-10), cv.ShouldEqual, true)
			cv.So(EpsSliceEqual(sd, kgSd, 1e-10), cv.ShouldEqual, true)
		})

	})
}

func TestCapnpLsqSerialization1(t *testing.T) {

	nvar := 7
	preDisk := NewMillerLSQ(nvar, 1)

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("in fuelconst.dat, df.Rows[0] = %v\n", df.Rows[0])
	last := df.Ncol - 1
	for i := range df.Rows {
		preDisk.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:], NAN_TO_ZERO)
	}

	// now save and restore

	fn := "fuel_fit_capn.gob"
	preDisk.GobEncode(fn)

	m := ReadLSQGobFile(fn)

	mat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)

	emat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormatE)
	//fmt.Printf("R is \n%v\n", mat)

	cv.Convey("Given an LSQ restored from disk file (with the fuelcons.dat data loaded, after 48 Includ() calls but before any other calls)", t, func() {

		cv.Convey("Then CompareLSQ should find no differences before / after the disk serialization", func() {
			// CompareLSQ panics if there is a difference, o/w returns true.
			match := CompareLSQ(preDisk, m)
			cv.So(match, cv.ShouldEqual, true)
		})
		cv.Convey("Then the previously seen lsq.Rows should be 48", func() {
			cv.So(m.Nobs, cv.ShouldEqual, 48)
		})

		cv.Convey("Then the R matrix should match the lsq.f90 version output", func() {
			cv.So(mat, cv.ShouldEqual, `
[ 1.000000 4296.916667 7.668333 2362.083333 4.241833 5.565417 2253.041667 57.033333 
 0.000000 1.000000 -0.000031 0.532401 0.000053 0.000508 0.463974 -0.000456 
 0.000000 0.000000 1.000000 -87.904537 0.044881 -1.603833 -215.545385 -2.036955 
 0.000000 0.000000 0.000000 1.000000 -0.000000 -0.000707 1.270032 0.007452 
 0.000000 0.000000 0.000000 0.000000 1.000000 -1.347978 -260.251651 3.888110 
 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 54.974483 0.547643 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 -0.006200 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000]`)

			// check more precisely for those small entries like row 4, column 5
			cv.So(emat, cv.ShouldEqual, `
[ 1.000000e+00 4.296917e+03 7.668333e+00 2.362083e+03 4.241833e+00 5.565417e+00 2.253042e+03 5.703333e+01 
 0.000000e+00 1.000000e+00 -3.139515e-05 5.324008e-01 5.295497e-05 5.081340e-04 4.639739e-01 -4.564454e-04 
 0.000000e+00 0.000000e+00 1.000000e+00 -8.790454e+01 4.488092e-02 -1.603833e+00 -2.155454e+02 -2.036955e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -4.463739e-09 -7.071586e-04 1.270032e+00 7.451966e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -1.347978e+00 -2.602517e+02 3.888110e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 5.497448e+01 5.476433e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -6.199548e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00]`)
		})

		u := m.XStats.Mean()
		sd := m.XStats.Sd()

		//fmt.Printf("u = %v\n", u)
		//fmt.Printf("sd = %v\n", sd)

		kgMean := []float64{4296.916666666666, 7.668333333333333, 2362.0833333333335, 4.241833333333332, 5.565416666666667, 2253.0416666666665, 57.033333333333324}
		kgSd := []float64{4441.10870206372, 0.9507697516051801, 2384.9691182561382, 0.5736237677697418, 3.491507166078876, 2124.270577772204, 5.547026549972453}
		cv.Convey("The mean and stddev should be what we expect", func() {
			//fmt.Printf("u = %v\n", u)
			//fmt.Printf("kgMean = %v\n", kgMean)
			cv.So(EpsSliceEqual(u, kgMean, 1e-10), cv.ShouldEqual, true)
			cv.So(EpsSliceEqual(sd, kgSd, 1e-10), cv.ShouldEqual, true)
		})

		cv.Convey("The preDisk.Rhs and m.Rhs be what we expect", func() {
			cv.So(len(m.Rhs), cv.ShouldEqual, len(preDisk.Rhs))
			for i := range preDisk.Rhs {
				cv.So(len(m.Rhs[i]), cv.ShouldEqual, len(preDisk.Rhs[i]))
				for j := range preDisk.Rhs[i] {
					cv.So(preDisk.Rhs[i][j], cv.ShouldEqual, m.Rhs[i][j])
				}
			}
		})

		cv.Convey("The preDisk.Sserr and m.Sserr be what we expect", func() {
			cv.So(len(m.Sserr), cv.ShouldEqual, len(preDisk.Sserr))
			for i := range preDisk.Sserr {
				cv.So(preDisk.Sserr[i], cv.ShouldEqual, m.Sserr[i])
			}
		})
	})
}
