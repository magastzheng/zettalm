package lsq

import (
	"fmt"
	"math"
	"testing"

	cv "github.com/glycerine/goconvey/convey"
)

// NIST regression standard, Norris data
//   should give:
/*                                   Standard Deviation
Parameter          Estimate             of Estimate

   B0        -0.262323073774029     0.232818234301152
   B1         1.00211681802045      0.429796848199937E-03

Residual
Standard Deviation   0.884796396144373

R-Squared            0.999993745883712
*/
// http://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Norris.dat
const expectedB0 = -0.262323073774029
const expectedB0_se = 0.232818234301152
const expectedB1 = 1.00211681802045
const expectedB1_se = 0.429796848199937E-03
const expectedResidSd = 0.884796396144373
const expectedR2 = 0.999993745883712

// ordering: y, then x
var Norris = []float64{
	0.1, 0.2,
	338.8, 337.4,
	118.1, 118.2,
	888.0, 884.6,
	9.2, 10.1,
	228.1, 226.5,
	668.5, 666.3,
	998.5, 996.3,
	449.1, 448.6,
	778.9, 777.0,
	559.2, 558.2,
	0.3, 0.4,
	0.1, 0.6,
	778.1, 775.5,
	668.8, 666.9,
	339.3, 338.0,
	448.9, 447.5,
	10.8, 11.6,
	557.7, 556.0,
	228.3, 228.1,
	998.0, 995.8,
	888.8, 887.6,
	119.6, 120.2,
	0.3, 0.3,
	0.6, 0.3,
	557.6, 556.8,
	339.3, 339.1,
	888.0, 887.2,
	998.5, 999.0,
	778.9, 779.0,
	10.2, 11.1,
	117.6, 118.3,
	228.9, 229.2,
	668.4, 669.1,
	449.2, 448.9,
	0.2, 0.5}

// given y, x, fit the model y = B0 + B1*x + B2*(x**2) + B3*(x**3)+ B4*(x**4) + B5*(x**5)
// see http://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler5.dat
/*
               Certified Regression Statistics

                                          Standard Deviation
     Parameter          Estimate             of Estimate

        B0        1.00000000000000         21523262.4678170
        B1        1.00000000000000         23635517.3469681
        B2        1.00000000000000         7793435.24331583
        B3        1.00000000000000         1014755.07550350
        B4        1.00000000000000         56456.6512170752
        B5        1.00000000000000         1123.24854679312

     Residual
     Standard Deviation   23601450.2379268

     R-Squared            0.224668921574940E-02


               Certified Analysis of Variance Table

Source of Degrees of    Sums of               Mean
Variation  Freedom      Squares              Squares           F Statistic

Regression    5    18814317208116.7      3762863441623.33 6.7552445824012241E-03
Residual     15    0.835542680000000E+16 557028453333333.


*/

var wampler5 = []float64{
	7590001, 0,
	-20479994, 1,
	20480063, 2,
	-20479636, 3,
	25231365, 4,
	-20476094, 5,
	20489331, 6,
	-20460392, 7,
	18417449, 8,
	-20413570, 9,
	20591111, 10,
	-20302844, 11,
	18651453, 12,
	-20077766, 13,
	21059195, 14,
	-19666384, 15,
	26348481, 16,
	-18971402, 17,
	22480719, 18,
	-17866340, 19,
	10958421, 20}

// obsoleted by the qr2cov_test.go code
func TestMerge(t *testing.T) {

	N := len(Norris) / 2
	//N := 4

	nxvar := 1
	nyvar := 1
	m1 := NewMillerLSQ(nxvar, nyvar)
	m2 := NewMillerLSQ(nxvar, nyvar)
	mtot := NewMillerLSQ(nxvar, nyvar)

	half := N / 2

	fmt.Printf(" =========== computing mtot\n")
	for i := 0; i < N; i++ {
		k := i * 2
		mtot.Includ(1.0, Norris[k+1:k+2], Norris[k:k+1], NAN_OMIT_ROW)
	}

	fmt.Printf(" =========== computing m1\n")
	for i := 0; i < N; i++ {
		k := i * 2
		if i < half {
			m1.Includ(1.0, Norris[k+1:k+2], Norris[k:k+1], NAN_OMIT_ROW)
		}
	}

	fmt.Printf(" =========== computing m2\n")
	for i := 0; i < N; i++ {
		k := i * 2
		if i >= half {
			m2.Includ(1.0, Norris[k+1:k+2], Norris[k:k+1], NAN_OMIT_ROW)
		}
	}

	fmt.Printf("%d rows in m1\n", m1.Nobs)
	fmt.Printf("%d rows in m2\n", m2.Nobs)
	fmt.Printf("%d rows in mtot\n", mtot.Nobs)

	fmt.Printf("m1.AccumWeightSum = %v\n", m1.AccumWeightSum)
	fmt.Printf("m2.AccumWeightSum = %v\n", m2.AccumWeightSum)

	mmerged := LsqCombineAllRhs(m1, m2)

	fmt.Printf("%d rows in mmerged after merge\n", mmerged.Nobs)
	fmt.Printf("mmerged.AccumWeightSum = %v\n", mmerged.AccumWeightSum)

	totR := UpperTriToString(mtot.R, mtot.Ncol, ";", StdFormat6dec)
	mergedR := UpperTriToString(mmerged.R, mmerged.Ncol, ";", StdFormat6dec)
	m1R := UpperTriToString(m1.R, m1.Ncol, ";", StdFormat6dec)
	m2R := UpperTriToString(m2.R, m2.Ncol, ";", StdFormat6dec)

	fmt.Printf("mtot.R = %v\n", totR)
	fmt.Printf("m1.R = %v\n", m1R)
	fmt.Printf("m2.R = %v\n", m2R)
	fmt.Printf("merged.R = %v\n", mergedR)

	fmt.Printf("mtot.D = %#v\n", mtot.D)
	fmt.Printf("m1.D = %#v\n", m1.D)
	fmt.Printf("m2.D = %#v\n", m2.D)
	fmt.Printf("mmerged.D = %#v\n", mmerged.D)

	//origRmatrix := UpperTriToString(m.R, m.Ncol, ";", StdFormat2dec)

	cv.Convey("Given that we compute each half of the Norris data in separate MillerLSQ", t, func() {
		cv.Convey("when we Merge them, we should get the same R, D, Rhs, Regcf as when we compute them all at once.", func() {

			cv.So(CompareLSQ(mtot, mmerged), cv.ShouldBeTrue)
			err, betaTot := mtot.Regcf([]int{0}, 0)
			if err != nil {
				panic(err)
			}
			mmerged.SS(0)
			err, betaMer := mmerged.Regcf([]int{0}, 0)
			if err != nil {
				panic(err)
			}
			floatcheck(betaTot, betaMer, "regression betas from mtot, mmerged")

			fmt.Printf("betaMer = %#v\n", betaMer)
			cv.So(fEquals(betaMer[0], expectedB0), cv.ShouldBeTrue)
			cv.So(fEquals(betaMer[1], expectedB1), cv.ShouldBeTrue)

			nreq := 2 // CONST, x
			covmat := make([]float64, nreq*(nreq+1)/2)
			sterr := make([]float64, nreq)
			wycol := 0
			err, Var2 := mmerged.Cov(nreq, covmat, sterr, wycol)
			if err != nil {
				panic(err)
			}
			cv.So(fEquals(sterr[0], expectedB0_se), cv.ShouldBeTrue)
			cv.So(fEquals(sterr[1], expectedB1_se), cv.ShouldBeTrue)
			fmt.Printf("covmat is %v\n", covmat)
			fmt.Printf("sterr is %v\n", sterr)
			fmt.Printf("Var2 is %v\n", Var2)
			cv.So(fEquals(math.Sqrt(Var2), expectedResidSd), cv.ShouldBeTrue)
			r2 := 1.0 - mmerged.Rss[wycol][nreq-1]/mmerged.Rss[wycol][0]
			fmt.Printf("R2 (r-squared) is %v\n", r2)
			cv.So(fEquals(r2, expectedR2), cv.ShouldBeTrue)

		})
	})

}

func TestShowRhsOnSmallFuelData(t *testing.T) {

	nvar := 2
	nyvar := 1
	m := NewMillerLSQ(nvar, nyvar)

	df, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("in smallfuel.dat, df.Rows[0] = %v\n", df.Rows[0])
	iniy := df.Ncol - 2
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][0:iniy], df.Rows[i][iniy:(iniy+1)], NAN_OMIT_ROW)

		fmt.Printf("\n after Includ()-ing row %d of smallfuel.dat (x=%v  y=%v), R is:\n", i, df.Rows[i][0:iniy], df.Rows[i][iniy:(iniy+1)])
		mat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)
		fmt.Printf("R is \n%v\n", mat)

		fmt.Printf("len(D)=%d   and D is \n%v\n", len(m.D), m.D)
		fmt.Printf("sqrt(D) is [")
		for k := range m.D {
			fmt.Printf("%v ", math.Sqrt(m.D[k]))
		}
		fmt.Printf("]\n")
		u := m.XStats.Mean()
		fmt.Printf("u = %v\n", u)
	}

	u := m.XStats.Mean()
	sd := m.XStats.Sd()
	m.SS(0)

	fmt.Printf("u = %v\n", u)
	fmt.Printf("sd = %v\n", sd)

	fmt.Printf("after column 3 only as y, Rhs = %v\n", m.Rhs)
	fmt.Printf("after column 3 only as y, Rss = %v\n", m.Rss)
	fmt.Printf("after column 3 only as y, D   = %v\n", m.D)
	fmt.Printf("after column 3 only as y, Sserr = %v\n", m.Sserr)

}

func TestN1N2vsMergeSmallFuelData(t *testing.T) {
	fmt.Printf("\n\nTestShowRhsOn3XSmallFuelData:\n\n")

	df, err := readData("smallfuel.dat")
	if err != nil {
		panic(err)
	}

	fmt.Printf("in smallfuel.dat, using 1st half, rows 0..4\n")
	QR(0, 4, df)

	fmt.Printf("in smallfuel.dat, using 2nd half, rows 5..10\n")
	QR(5, 9, df)

	fmt.Printf("in smallfuel.dat, using all rows 0..10\n")
	QR(0, 9, df)

}

func QR(start int, end int, df *DataFrame) {
	nvar := 3
	nyvar := 1
	m := NewMillerLSQ(nvar, nyvar)

	iniy := df.Ncol - 1
	for i := start; i <= end; i++ {
		m.Includ(1.0, df.Rows[i][0:iniy], df.Rows[i][iniy:(iniy+1)], NAN_OMIT_ROW)
		fmt.Printf("Includ()ing %v %v\n", df.Rows[i][0:iniy], df.Rows[i][iniy:(iniy+1)])
	}

	fmt.Printf("\nNobs = %d\n", m.Nobs)
	mat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)
	fmt.Printf("R is \n%v\n", mat)

	u := m.XStats.Mean()
	sd := m.XStats.Sd()
	m.SS(0)

	fmt.Printf("u = %v\n", u)
	fmt.Printf("sd = %v\n", sd)

	fmt.Printf("Rhs = %v\n", m.Rhs)
	fmt.Printf("Rss = %v\n", m.Rss)
	fmt.Printf("D   = %v\n", m.D)
	fmt.Printf("Sserr = %v\n", m.Sserr)

	fmt.Printf("\nAssembled R with sqrt(D) on diagonal is\n%s\n\n", m.StringRDaigonalSqrtD())

	fmt.Printf("\nRaw R is: %v\n\n", m.R)

	wycol := 0
	err, beta := m.Regcf([]int{1, 2, 3}, wycol)
	if err != nil {
		panic(err)
	}
	fmt.Printf("\n beta = %v\n", beta)

	nreq := 3
	covmat := make([]float64, nreq*(nreq+1)/2)
	sterr := make([]float64, nreq)
	err, Var2 := m.Cov(nreq, covmat, sterr, wycol)

	fmt.Printf("\n covmat =\n%v\n", covmat)

	fmt.Printf("\n sterr  = %v\n", sterr)

	fmt.Printf("\n Var2 = %v\n", Var2)
}
