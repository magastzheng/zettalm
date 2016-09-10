package lsq

// convert from a QR decomposition in MillerLSQ to a Covariance matrix

import (
	"fmt"
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

func TestQR2Cov(t *testing.T) {
	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}
	x1y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1y.Rows[i] = xf.Rows[i]
	}

	nxvar := xf.Ncol - 1
	nyvar := 1
	qr1 := NewMillerLSQ(nxvar, nyvar)
	for i := range x1y.Rows {
		qr1.Includ(1.0, x1y.Rows[i][:nxvar], x1y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	sup := SupplementR(&qr1.UpperTriMatrix, qr1.Rhs[0], qr1.Sserr[0], NoRecycle)

	qr1str := sup.DecodedUpperTriToString("\n", StdFormat6dec)

	// qr1.UpperTriMatrix is only 4x4 now, we want the 5x5 that includes Rhs[0]:
	// qr1str := qr1.DecodedUpperTriToString("\n", StdFormat6dec)

	cv.Convey("Given we've computed the online QR decomposition from smallfuel.dat (first 5 rows) using lsq, our QR should match R's", t, func() {
		cv.So(qr1str, cv.ShouldEqual, `matrix(ncol=5, byrow=TRUE, data=c(
 2.236068, 3.396140, 1794.668159, 122.983739, 1095.673309, 
 0.000000, 1.470261, 1173.645253, -1.787234, 29.060151, 
 0.000000, 0.000000, 1356.721792, -2.136710, -138.149405, 
 0.000000, 0.000000, 0.000000, 4.135247, 30.143798, 
 0.000000, 0.000000, 0.000000, 0.000000, 11.644788))`)
	})

	fmt.Printf("qr1str = %v\n", qr1str)

	covFromQR1, mean1qr, N1qr := qr1.QR2Cov(0)
	fmt.Printf("covFromQR1 = %v\n", covFromQR1)
	fmt.Printf("mean1qr = %v\n", mean1qr)
	fmt.Printf("N1qr = %v\n", N1qr)

	//cov1xx, mean1xx, N1xx := CovAddConstKeepMeanColumn(x1y)
	x1y.AddFirstOnesColumn()
	//cov1xx, mean1xx, N1xx := Cov(x1y)
	cov1xx, mean1xx, N1xx := Cov(x1y)
	fmt.Printf("cov1xx.Colnames = %#v\n", cov1xx.Colnames)
	fmt.Printf("cov1xx = %v\n", cov1xx)
	fmt.Printf("mean1xx = %v\n", mean1xx)
	fmt.Printf("N1xx = %v\n", N1xx)
	fmt.Printf("x1y after AddFirstOnesColumn() = %v\n", x1y.String())

	fmt.Printf("cov1xx = %v\n", cov1xx)
	fmt.Printf("len(cov1xx.A) = %v\n", len(cov1xx.A))

	cv.Convey("Given we've computed the online QR decomposition from smallfuel.dat (first 5 row) using lsq", t, func() {
		cv.Convey("then the QR2Cov() function should return a Cov matrix that matches that computed from X'X directly.", func() {
			cv.So(covFromQR1.String(), cv.ShouldEqual, `matrix(ncol=5, nrow=5, byrow=TRUE, data=c(
 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, 
 -0.000000, 0.540417, 431.391150, -0.656925, 10.681500, 
 0.000000, 431.391150, 804534.300000, -1249.125000, -38331.000000, 
 0.000000, -0.656925, -1249.125000, 6.215000, 91.975000, 
 0.000000, 10.681500, -38331.000000, 91.975000, 5243.500000))`)

			// mismatch only by -0 vs 0
			//cv.So(covFromQR1.String(), cv.ShouldEqual, cov1xx.String())
			cv.So(cov1xx.String(), cv.ShouldEqual, `matrix(ncol=5, nrow=5, byrow=TRUE, data=c(
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
 0.000000, 0.540417, 431.391150, -0.656925, 10.681500, 
 0.000000, 431.391150, 804534.300000, -1249.125000, -38331.000000, 
 0.000000, -0.656925, -1249.125000, 6.215000, 91.975000, 
 0.000000, 10.681500, -38331.000000, 91.975000, 5243.500000))`)

			cv.So(len(covFromQR1.A), cv.ShouldEqual, len(cov1xx.A))
			for i := range covFromQR1.A {
				cv.So(covFromQR1.A[i], cv.ShouldAlmostEqual, cov1xx.A[i], 1e-6)
			}
			for i := range covFromQR1.D {
				cv.So(covFromQR1.D[i], cv.ShouldAlmostEqual, cov1xx.D[i], 1e-6)
			}
			for i := range mean1qr {
				cv.So(mean1qr[i], cv.ShouldAlmostEqual, mean1xx[i], 1e-6)
			}
			cv.So(N1qr, cv.ShouldEqual, N1xx)
		})
	})
}

func TestCov2QR(t *testing.T) {
	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}
	x1y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1y.Rows[i] = xf.Rows[i]
	}

	nxvar := xf.Ncol - 1
	nyvar := 1
	qr1 := NewMillerLSQ(nxvar, nyvar)
	for i := range x1y.Rows {
		qr1.Includ(1.0, x1y.Rows[i][:nxvar], x1y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	//qr1str := qr1.DecodedUpperTriToString("\n", StdFormat6dec)
	sup := SupplementR(&qr1.UpperTriMatrix, qr1.Rhs[0], qr1.Sserr[0], NoRecycle)
	qr1str := sup.DecodedUpperTriToString("\n", StdFormat6dec)

	cv.Convey("Given we've computed the online QR decomposition from smallfuel.dat (first 5 rows) using lsq, our QR should match R's", t, func() {
		cv.So(qr1str, cv.ShouldEqual, `matrix(ncol=5, byrow=TRUE, data=c(
 2.236068, 3.396140, 1794.668159, 122.983739, 1095.673309, 
 0.000000, 1.470261, 1173.645253, -1.787234, 29.060151, 
 0.000000, 0.000000, 1356.721792, -2.136710, -138.149405, 
 0.000000, 0.000000, 0.000000, 4.135247, 30.143798, 
 0.000000, 0.000000, 0.000000, 0.000000, 11.644788))`)
	})

	fmt.Printf("qr1.R = %#v\n", qr1.R)
	fmt.Printf("qr1.D = %#v\n", qr1.D)
	notdecoded_qr1str := sup.UpperTriToString("\n", StdFormat6dec)
	fmt.Printf("notdecoded_qr1str = %s\n", notdecoded_qr1str)

	fmt.Printf("qr1str = %v\n", qr1str)

	covFromQR1, mean1qr, N1qr := qr1.QR2Cov(0)
	fmt.Printf("covFromQR1 = %v\n", covFromQR1)
	fmt.Printf("mean1qr = %v\n", mean1qr)
	fmt.Printf("N1qr = %v\n", N1qr)

	nc := nxvar + 1 + nyvar
	m := NewUpperTriMatrix(nc)
	m.Cov2QR(covFromQR1, N1qr, mean1qr[:(nxvar+1)], mean1qr[(nxvar+1)], !SetQRCompressed)
	restr := m.String6()

	fmt.Printf("expected: \n%s\n", qr1str)
	fmt.Printf("got: restr from recreateQR1 using Cov2QR() is:\n%s\n", restr)
	if qr1str == restr {
		fmt.Printf(" and qr1str == restr, so we're in good shape wth !SetQRCompressed.")
	}

	cv.Convey("Given smallfuel.dat (first 5 rows), starting from QR2Cov() output,", t, func() {
		cv.Convey("then the inverse Cov2QR2() should re-create the MillerLSQ UpperTriMatrix exactly.", func() {
			cv.So(restr, cv.ShouldEqual, qr1str)
		})
	})

	// compressed_qr1str shows it just as it is stored.
	compressed_qr1str := qr1.UpperTriToString("\n", StdFormat6dec)
	if compressed_qr1str == qr1str {
		panic("compressed_qr1str ought to be different from qr1str")
	}
	//fmt.Printf("dump qr1 = \n")
	//goon.Dump(qr1)
	//fmt.Printf("compressed_qr1str = %s\n", compressed_qr1str)

	m2 := NewUpperTriMatrix(nc)
	m2.Cov2QR(covFromQR1, N1qr, mean1qr[:(nxvar+1)], mean1qr[(nxvar+1)], SetQRCompressed)
	restr2 := m2.String6()

	desup, rhs, sserr := m2.DeSupplementR(NoRecycle)
	desupStr := desup.UpperTriToString("\n", StdFormat6dec)

	cv.Convey("Given smallfuel.dat (first 5 rows), starting from QR2Cov() output,", t, func() {
		cv.Convey("when we ask for SetQRCompressed, the inverse Cov2QR2() should re-create the (compressed format) MillerLSQ UpperTriMatrix exactly.", func() {
			cv.So(restr2, cv.ShouldEqual, `matrix(ncol=5, byrow=TRUE, data=c(
 5.000000, 1.518800, 802.600000, 55.000000, 490.000000, 
 0.000000, 2.161667, 798.256512, -1.215590, 19.765303, 
 0.000000, 0.000000, 1840694.021073, -0.001575, -0.101826, 
 0.000000, 0.000000, 0.000000, 17.100264, 7.289480, 
 0.000000, 0.000000, 0.000000, 0.000000, 135.601090))`)
			cv.So(rhs, cv.ShouldResemble, []float64{490, 19.76530332981944, -0.10182589044590314, 7.289480300957488})
			cv.So(sserr, cv.ShouldEqual, 135.60108986194243)

			cv.So(desupStr, cv.ShouldEqual, compressed_qr1str)
			cv.So(desupStr, cv.ShouldEqual, `matrix(ncol=4, byrow=TRUE, data=c(
 5.000000, 1.518800, 802.600000, 55.000000, 
 0.000000, 2.161667, 798.256512, -1.215590, 
 0.000000, 0.000000, 1840694.021073, -0.001575, 
 0.000000, 0.000000, 0.000000, 17.100264))`)
		})
	})

}

func TestMergeOfTwoQR(t *testing.T) {
	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	// x1y
	x1y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1y.Rows[i] = xf.Rows[i]
	}

	nxvar := xf.Ncol - 1
	nyvar := 1
	qr1 := NewMillerLSQ(nxvar, nyvar)
	for i := range x1y.Rows {
		qr1.Includ(1.0, x1y.Rows[i][:nxvar], x1y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// covFromQR1
	covFromQR1, mean1qr, N1qr := qr1.QR2Cov(0)
	fmt.Printf("covFromQR1 = %v\n", covFromQR1)
	fmt.Printf("mean1qr = %v\n", mean1qr)
	fmt.Printf("N1qr = %v\n", N1qr)

	// x2y
	x2y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x2y.Rows[i] = xf.Rows[i+5]
	}

	qr2 := NewMillerLSQ(nxvar, nyvar)
	for i := range x2y.Rows {
		qr2.Includ(1.0, x2y.Rows[i][:nxvar], x2y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// covFromQR2
	covFromQR2, mean2qr, N2qr := qr2.QR2Cov(0)
	fmt.Printf("covFromQR2 = %v\n", covFromQR2)
	fmt.Printf("mean2qr = %v\n", mean2qr)
	fmt.Printf("N2qr = %v\n", N2qr)

	// merge covariances
	covMerged, meanc, combinedN := CombineCov(covFromQR1, covFromQR2, mean1qr, mean2qr, N1qr, N2qr)
	strCovMerged := covMerged.String()

	u := NewUpperTriMatrix(covFromQR1.Ncol)
	u.Cov2QR(covMerged, combinedN, meanc[:(nxvar+1)], meanc[(nxvar+1)], SetQRCompressed)

	desupMerged, rhs2Merged, sserrMerged := u.DeSupplementR(NoRecycle)
	desupMergedStr := desupMerged.String6()

	mergedRDStr := u.String6()
	decodedMergedRDStr := u.DecodedUpperTriToString("\n", StdFormat6dec)

	// full, sequential, the gold-standard to match against
	qrf := NewMillerLSQ(nxvar, nyvar)
	for i := range xf.Rows {
		qrf.Includ(1.0, xf.Rows[i][:nxvar], xf.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	cv.Convey("Given smallfuel.dat (first 5 rows = x1y, second 5 rows = x2y), starting from QR2Cov() from x1y and x2y,", t, func() {
		cv.Convey("when we merge the two Cov matrixes, and then backtransform into a cholesky decomposition, we should get what we expect in the UpperTriMatrix", func() {
			cv.So(strCovMerged, cv.ShouldEqual, `matrix(ncol=5, nrow=5, byrow=TRUE, data=c(
 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, 
 -0.000000, 16.298971, 8972.802722, -10.597169, -117.526333, 
 0.000000, 8972.802722, 5727594.055556, -5502.233333, -80835.000000, 
 0.000000, -10.597169, -5502.233333, 13.664889, 176.277778, 
 0.000000, -117.526333, -80835.000000, 176.277778, 4423.111111))`)
			cv.So(mergedRDStr, cv.ShouldEqual, `matrix(ncol=5, byrow=TRUE, data=c(
 10.000000, 4.001700, 2607.500000, 54.060000, 468.000000, 
 0.000000, 146.690740, 550.513444, -0.650174, -7.210660, 
 0.000000, 0.000000, 7091509.758978, 0.000421, -0.020478, 
 0.000000, 0.000000, 0.000000, 59.717615, 16.074139, 
 0.000000, 0.000000, 0.000000, 0.000000, 13777.628216))`)
			cv.So(decodedMergedRDStr, cv.ShouldEqual, `matrix(ncol=5, byrow=TRUE, data=c(
 3.162278, 12.654487, 8245.638999, 170.952730, 1479.945945, 
 0.000000, 12.111595, 6667.596024, -7.874646, -87.332591, 
 0.000000, 0.000000, 2662.988877, 1.120867, -54.531419, 
 0.000000, 0.000000, 0.000000, 7.727717, 124.216403, 
 0.000000, 0.000000, 0.000000, 0.000000, 117.378142))`)
		})

		cv.Convey("and these should match against the cov and qr generated by a full-sequential processing of all 10 rows in smallfuel.dat", func() {

			for i := range qrf.R {
				cv.So(desupMerged.R[i], cv.ShouldAlmostEqual, qrf.R[i], 1e-6)
			}
			for i := range qrf.D {
				cv.So(desupMerged.D[i], cv.ShouldAlmostEqual, qrf.D[i], 1e-6)
			}

			for i := range qrf.Rhs[0] {
				cv.So(rhs2Merged[i], cv.ShouldAlmostEqual, qrf.Rhs[0][i], 1e-6)
			}
			cv.So(sserrMerged, cv.ShouldAlmostEqual, qrf.Sserr[0], 1e-6)

			cv.So(desupMergedStr, cv.ShouldEqual, `matrix(ncol=4, byrow=TRUE, data=c(
 10.000000, 4.001700, 2607.500000, 54.060000, 
 0.000000, 146.690740, 550.513444, -0.650174, 
 0.000000, 0.000000, 7091509.758978, 0.000421, 
 0.000000, 0.000000, 0.000000, 59.717615))`)

			// covFromQRf
			covFromQRf, meanfqr, Nfqr := qrf.QR2Cov(0)
			if false {
				fmt.Printf("covFromQRf = %v\n", covFromQRf)
				fmt.Printf("meanfqr = %v\n", meanfqr)
				fmt.Printf("Nfqr = %v\n", Nfqr)
			}

			for i := range covFromQRf.A {
				cv.So(covFromQRf.A[i], cv.ShouldAlmostEqual, covMerged.A[i], 1e-6)
			}
			for i := range covFromQRf.D {
				cv.So(covFromQRf.D[i], cv.ShouldAlmostEqual, covMerged.D[i], 1e-6)
			}

		})
	})

}

func TestCombineOfTwoLSQ(t *testing.T) {
	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	// x1y
	x1y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1y.Rows[i] = xf.Rows[i]
	}

	nxvar := xf.Ncol - 1
	nyvar := 1
	qr1 := NewMillerLSQ(nxvar, nyvar)
	for i := range x1y.Rows {
		qr1.Includ(1.0, x1y.Rows[i][:nxvar], x1y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// x2y
	x2y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x2y.Rows[i] = xf.Rows[i+5]
	}

	qr2 := NewMillerLSQ(nxvar, nyvar)
	for i := range x2y.Rows {
		qr2.Includ(1.0, x2y.Rows[i][:nxvar], x2y.Rows[i][nxvar:], NAN_OMIT_ROW)

	}

	// full, sequential, the gold-standard to match against
	qrf := NewMillerLSQ(nxvar, nyvar)
	for i := range xf.Rows {
		qrf.Includ(1.0, xf.Rows[i][:nxvar], xf.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// the main event
	lsqRe := LsqCombineAllRhs(qr1, qr2)

	// now check it out
	verified := CompareLSQ(lsqRe, qrf)

	cv.Convey("We should be able re-create a MillerLSQ from merged data that matches the sequentially computed state (including Rhs[] and Sserr[])", t, func() {
		cv.So(verified, cv.ShouldEqual, true)
	})
}

func TestCombineOfTwoLSQWithTwoY(t *testing.T) {

	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	// x1y
	x1y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1y.Rows[i] = xf.Rows[i]
	}

	nxvar := xf.Ncol - 2
	nyvar := 2
	qr1 := NewMillerLSQ(nxvar, nyvar)
	for i := range x1y.Rows {
		qr1.Includ(1.0, x1y.Rows[i][:nxvar], x1y.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// x2y
	x2y := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x2y.Rows[i] = xf.Rows[i+5]
	}

	qr2 := NewMillerLSQ(nxvar, nyvar)
	for i := range x2y.Rows {
		qr2.Includ(1.0, x2y.Rows[i][:nxvar], x2y.Rows[i][nxvar:], NAN_OMIT_ROW)

	}

	// full, sequential, the gold-standard to match against
	qrf := NewMillerLSQ(nxvar, nyvar)
	for i := range xf.Rows {
		qrf.Includ(1.0, xf.Rows[i][:nxvar], xf.Rows[i][nxvar:], NAN_OMIT_ROW)
	}

	// the main event
	lsqRe := LsqCombineAllRhs(qr1, qr2)

	// now check it out
	verified := CompareLSQ(lsqRe, qrf)

	cv.Convey("We should be able re-create a MillerLSQ from merged data that matches the sequentially computed state (including Rhs[] and Sserr[])", t, func() {
		cv.So(verified, cv.ShouldEqual, true)
	})

}
