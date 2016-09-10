package lsq

import "fmt"

// combine two predictions into one, allowing parallel learning

func LsqCombineAllRhs(lsq1 *MillerLSQ, lsq2 *MillerLSQ) (merged *MillerLSQ) {

	nyvar := lsq1.Nyvar
	nxvar := lsq1.Nxvar
	combonc := nxvar + 2 // + 2 because one intercept, one y

	comboSserr := make([]float64, nyvar)
	comboRhs := make([][]float64, nyvar)

	newRbig := NewUpperTriMatrix(combonc)       // overwritten in place each wycol
	newRsmall := NewUpperTriMatrix(combonc - 1) // overwritten in place each wycol

	merged = PrepNewMergedLsq(lsq1, lsq2)

	//fmt.Printf("\n before wycol loop, lsq1 = %#v\n", lsq1)
	//fmt.Printf("\n before wycol loop, lsq2 = %#v\n", lsq2)

	for wycol := 0; wycol < nyvar; wycol++ {
		crhs, csserr := LsqCombine1Rhs(lsq1, lsq2, wycol, newRbig, newRsmall)
		comboSserr[wycol] = csserr
		comboRhs[wycol] = crhs
	}

	merged.UpperTriMatrix = *newRsmall
	merged.Rhs = comboRhs
	merged.Sserr = comboSserr

	//fmt.Printf("\n after wycol loop, merged = %#vn", merged)

	return merged
}

func LsqCombine1Rhs(lsq1 *MillerLSQ, lsq2 *MillerLSQ, wycol int, upbig *UpperTriMatrix, writetoSmall *UpperTriMatrix) (newRhs []float64, newSserr float64) {

	if lsq1.Nyvar != lsq2.Nyvar {
		panic(fmt.Sprintf("LsqCombine1Rhs() must be called on models over the same data set. Nyvar(%v) != lsq2.Nyvar(%v)", lsq1.Nyvar, lsq2.Nyvar))
	}
	if lsq1.Nxvar != lsq2.Nxvar {
		panic(fmt.Sprintf("LsqCombine1Rhs() must be called on models over the same data set. Nxvar(%v) != lsq2.Nxvar(%v)", lsq1.Nxvar, lsq2.Nxvar))
	}

	covFromQR1, mean1qr, N1qr := lsq1.QR2Cov(wycol) // allocates new SquareMatrix covFromQR1
	covFromQR2, mean2qr, N2qr := lsq2.QR2Cov(wycol) // allocates new SquareMatrix covFromQR2

	//fmt.Printf("\n lsq1.R = %v\n", lsq1.R)
	//fmt.Printf("\n lsq1.D = %v\n", lsq1.D)
	//fmt.Printf("\n lsq1.Rhs = %v\n", lsq1.Rhs)

	//fmt.Printf("\n mean1qr = %#v,  N1qr=%v\n", mean1qr, N1qr)
	//fmt.Printf("\n mean2qr = %#v,  N2qr=%v\n", mean2qr, N2qr)

	//fmt.Printf("\n covFromQR1 = %v\n", covFromQR1)
	//fmt.Printf("\n covFromQR2 = %v\n", covFromQR2)

	// merge covariances, compensating for the fact that they aren't mean-subtracted.
	// covMerged is a new SquareMatrix
	covMerged, meanc, combinedN := CombineCov(covFromQR1, covFromQR2, mean1qr, mean2qr, N1qr, N2qr)

	//fmt.Printf("\n covMerged = %v\n", covMerged)
	//fmt.Printf("\n meanc = %v\n", meanc)
	//fmt.Printf("\n combinedN = %v\n", combinedN)

	nc := len(meanc) - 1
	ymean := meanc[nc]
	//fmt.Printf("ymean = %v\n", ymean)
	upbig.Cov2QR(covMerged, combinedN, meanc[:nc], ymean, SetQRCompressed)

	uR, uRhs, uSserr := upbig.DeSupplementR(writetoSmall)

	if uR != writetoSmall {
		panic("AARG! writetoSmall didn't get recycled like it should have, we won't return the results we should have in it!!")
	}

	return uRhs, uSserr
}

// take care of the easy initial combination/merge stuff here.
func PrepNewMergedLsq(lsq1 *MillerLSQ, lsq2 *MillerLSQ) (merged *MillerLSQ) {

	merged = NewMillerLSQ(lsq1.Nxvar, lsq1.Nyvar)

	merged.RowsSeen = lsq1.RowsSeen + lsq2.RowsSeen
	merged.Nobs = lsq1.Nobs + lsq2.Nobs
	merged.CountNaNRowsSkipped = lsq1.CountNaNRowsSkipped + lsq2.CountNaNRowsSkipped

	merged.AccumWeightSum = lsq1.AccumWeightSum + lsq2.AccumWeightSum
	merged.Nxvar = lsq1.Nxvar
	merged.Nyvar = lsq1.Nyvar
	merged.NanApproach = lsq1.NanApproach
	merged.UseMeanSd = lsq1.UseMeanSd

	merged.XStats = lsq1.XStats
	merged.XStats.Merge(&lsq2.XStats)
	merged.YStats = lsq1.YStats
	merged.YStats.Merge(&lsq2.YStats)

	// Rss // filled in by SS() from D, Rhs, and Sserr

	// stuff I'm not using and therefore not worrying about merging at the moment.
	//	Xmean
	//	Xsd
	//	Ymean
	//	Ysd
	//  Vorder

	return merged
}

/*
	// covFromQR1
	covFromQR1, mean1qr, N1qr := qr1.QR2Cov(0)
	fmt.Printf("covFromQR1 = %v\n", covFromQR1)
	fmt.Printf("mean1qr = %v\n", mean1qr)
	fmt.Printf("N1qr = %v\n", N1qr)

	// covFromQR2
	covFromQR2, mean2qr, N2qr := qr2.QR2Cov(0)
	fmt.Printf("covFromQR2 = %v\n", covFromQR2)
	fmt.Printf("mean2qr = %v\n", mean2qr)
	fmt.Printf("N2qr = %v\n", N2qr)

	// merge covariances
	covMerged, meanc, combinedN := CombineCov(covFromQR1, covFromQR2, mean1qr, mean2qr, N1qr, N2qr)

	u := NewUpperTriMatrix(covFromQR1.Ncol)
	u.Cov2QR(covMerged, combinedN, nxvar, meanc, []float64{}, SetQRCompressed)
*/
