package lsq

import (
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

/*
$ cat smallfuel.dat
RoadMls FuelCon   DLic  Fuel_Pop
1.976    557.    52.5    541.
1.250    404.    57.2    524.
1.586    259.    58.0    561.
2.351   2396.    52.9    414.
.431    397.    54.4    410.
1.333   1408.    57.1    457.
11.868   6312.    45.1    344.
2.138   3439.    55.3    467.
8.577   5528.    52.9    464.
8.507   5375.    55.2    498.

# in R
> a=read.table("smallfuel.dat",header=T)
> x=cbind(rep(1,10),as.matrix(a))
> x
        RoadMls FuelCon DLic Fuel_Pop
 [1,] 1   1.976     557 52.5      541
 [2,] 1   1.250     404 57.2      524
 [3,] 1   1.586     259 58.0      561
 [4,] 1   2.351    2396 52.9      414
 [5,] 1   0.431     397 54.4      410
 [6,] 1   1.333    1408 57.1      457
 [7,] 1  11.868    6312 45.1      344
 [8,] 1   2.138    3439 55.3      467
 [9,] 1   8.577    5528 52.9      464
[10,] 1   8.507    5375 55.2      498
> t(x) %*% x
                       RoadMls     FuelCon        DLic    Fuel_Pop
            10.000     40.0170     26075.0     540.600     4680.00
RoadMls     40.017    306.8268    185099.6    2067.945    17670.22
FuelCon  26075.000 185099.5520 119538909.0 1360094.400 11475585.00
DLic       540.600   2067.9445   1360094.4   29347.820   254587.30
Fuel_Pop  4680.000  17670.2190  11475585.0  254587.300  2230048.00
>
*/

func TestXtX(t *testing.T) {
	x, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	xtx, _, _ := XtX(x, !RmMean, !AddIntercept)

	cv.Convey("When smallfuel.dat is loaded into X, t(X) %%*%% X, or XtX(x), should match R (don't remove mean, don't add intercept)", t, func() {
		cv.So(xtx.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 306.826769, 185099.552000, 2067.944500, 17670.219000, 
 185099.552000, 119538909.000000, 1360094.400000, 11475585.000000, 
 2067.944500, 1360094.400000, 29347.820000, 254587.300000, 
 17670.219000, 11475585.000000, 254587.300000, 2230048.000000))`)
	})

	xtxi, _, _ := XtX(x, !RmMean, AddIntercept)
	cv.Convey("When smallfuel.dat is loaded into X, t(X) %%*%% X, or XtX(x), should match R (add intercept)", t, func() {
		cv.So(xtxi.String(), cv.ShouldEqual, `matrix(ncol=5, nrow=5, byrow=TRUE, data=c(
 10.000000, 40.017000, 26075.000000, 540.600000, 4680.000000, 
 40.017000, 306.826769, 185099.552000, 2067.944500, 17670.219000, 
 26075.000000, 185099.552000, 119538909.000000, 1360094.400000, 11475585.000000, 
 540.600000, 2067.944500, 1360094.400000, 29347.820000, 254587.300000, 
 4680.000000, 17670.219000, 11475585.000000, 254587.300000, 2230048.000000))`)
	})

	xtxmi, _, _ := XtX(x, RmMean, AddIntercept)
	cv.Convey("When smallfuel.dat is loaded into X, t(X) %%*%% X, or XtX(x), should match R (remove mean, add intercept)", t, func() {
		cv.So(xtxmi.String(), cv.ShouldEqual, `matrix(ncol=5, nrow=5, byrow=TRUE, data=c(
 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
 0.000000, 146.690740, 80755.224500, -95.374520, -1057.737000, 
 0.000000, 80755.224500, 51548346.500000, -49520.100000, -727515.000000, 
 0.000000, -95.374520, -49520.100000, 122.984000, 1586.500000, 
 0.000000, -1057.737000, -727515.000000, 1586.500000, 39808.000000))`)
	})

	cov, _, _ := Cov(x)
	cv.Convey("When smallfuel.dat is loaded into X, t(X) %%*%% X, or XtX(x), should match R (remove mean, don't add intercept, divide by nr-1)", t, func() {
		cv.So(cov.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 16.298971, 8972.802722, -10.597169, -117.526333, 
 8972.802722, 5727594.055556, -5502.233333, -80835.000000, 
 -10.597169, -5502.233333, 13.664889, 176.277778, 
 -117.526333, -80835.000000, 176.277778, 4423.111111))`)
	})

	L, err, _ := cov.FactorToCholeskyLower(CHOL_RETURN_FRESH_MATRIX)
	if err != nil {
		panic(err)
	}
	cv.Convey("When smallfuel.dat is loaded into X, after Cov(x) if we FactorToCholeskyLower()", t, func() {
		cv.Convey("then the upper tri should be the same, and only the lower tri should be written to", func() {
			cv.So(L.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 4.037198, 0.000000, 0.000000, 0.000000, 
 2222.532008, 887.662959, 0.000000, 0.000000, 
 -2.624882, 0.373622, 2.575906, 0.000000, 
 -29.110864, -18.177140, 41.405468, 39.126047))`)
		})
	})

	L, err, _ = cov.FactorToCholeskyLower(CHOL_OVERWRITE_LOWER_AND_D)
	if err != nil {
		panic(err)
	}
	cv.Convey("When smallfuel.dat is loaded into X, after Cov(x) if we FactorToCholeskyLower()", t, func() {
		cv.Convey("then the upper tri should be the same, and only the lower tri should be written to", func() {
			cv.So(L.String(), cv.ShouldResemble, cov.String())

			// Show the cholesky diagonal, by default after CHOL_OVERWRITE_LOWER_AND_D.
			cv.So(L.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 4.037198, 8972.802722, -10.597169, -117.526333, 
 2222.532008, 887.662959, -5502.233333, -80835.000000, 
 -2.624882, 0.373622, 2.575906, 176.277778, 
 -29.110864, -18.177140, 41.405468, 39.126047))`)

			// now flip ShowD off and just show the variances on the diagonal
			L.ShowD = false
			cv.So(L.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 16.298971, 8972.802722, -10.597169, -117.526333, 
 2222.532008, 5727594.055556, -5502.233333, -80835.000000, 
 -2.624882, 0.373622, 13.664889, 176.277778, 
 -29.110864, -18.177140, 41.405468, 4423.111111))`)
			// the square roots from the cholesky are still there in D, just not shown by String() when ShowD is false.
			cv.So(L.D, cv.ShouldResemble, []float64{4.03719842492566, 887.6629589970097, 2.5759057775456307, 39.1260473287411})
		})
	})

}

func TestCholeskyWith_CHOL_OVERWRITE_VARIANCES(t *testing.T) {
	x, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	cov, _, _ := Cov(x)
	L, err, _ := cov.FactorToCholeskyLower(CHOL_OVERWRITE_VARIANCES)
	if err != nil {
		panic(err)
	}
	cv.Convey("When smallfuel.dat, after Cov(x) if we FactorToCholeskyLower(CHOL_OVERWRITE_VARIANCES)", t, func() {
		cv.Convey("then the diagonal should have the sqrt(Var) instead of Var, reflecting the cholesky results.", func() {
			cv.So(L.String(), cv.ShouldEqual, `matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
 4.037198, 8972.802722, -10.597169, -117.526333, 
 2222.532008, 887.662959, -5502.233333, -80835.000000, 
 -2.624882, 0.373622, 2.575906, 176.277778, 
 -29.110864, -18.177140, 41.405468, 39.126047))`)
			cv.So(L.D, cv.ShouldResemble, []float64{16.298971122222227, 5.727594055555556e+06, 13.664888888888889, 4423.111111111111})
			cv.So(L.String(), cv.ShouldResemble, cov.String())
		})
	})
}

func TestCovMerge(t *testing.T) {
	xf, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}
	x1 := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x1.Rows[i] = xf.Rows[i]
	}

	x2 := &DataFrame{
		Colnames: xf.Colnames,
		Ncol:     xf.Ncol,
		Nrow:     5,
		Rows:     make([][]float64, 5),
	}
	for i := 0; i < 5; i++ {
		x2.Rows[i] = xf.Rows[i+5]
	}

	cov1, mean1, N1 := Cov(x1)
	cov2, mean2, N2 := Cov(x2)

	// gold standard
	covf, meanf, Nf := Cov(xf)

	// what we want to compute to match covf, but in parallel
	covc, meanc, Nc := CombineCov(cov1, cov2, mean1, mean2, float64(N1), float64(N2))

	cv.Convey("Given covariance matrices computed on two halfs of the smallfuel.dat", t, func() {
		cv.Convey("then the CombinedCov() should match the sequentially computed full dataset variance.", func() {
			for i := range covc.A {
				cv.So(covc.A[i], cv.ShouldAlmostEqual, covf.A[i], 1e-6)
			}
			cv.So(Nc, cv.ShouldEqual, int(Nf))
			for i := range meanc {
				cv.So(meanc[i], cv.ShouldAlmostEqual, meanf[i], 1e-6)
			}
		})
	})
}
