package lsq

import "fmt"

// convert from a QR decomposition in MillerLSQ to a Covariance matrix
/*
R2cov=function(R, X) {

  nc = ncol(X)
  nr = nrow(X)
  umat=matrix(ncol=nc, nrow=nc,data=0)
  u = colMean(X)

  for (i in 1:nc) {
    for (j in 1:nc) {
      umat[i,j] = u[i]*u[j]*nr
    }
  }

  cv=(t(R) %*% R - umat)/(nr-1)

  cv
}
*/

//
// can't call MillerLSQ.Cov(), since that is the covariance of the betas,
//  not the covariance of the design matrix X.
//
// always expects a ymeans value now
//
func (qr *MillerLSQ) QR2Cov(wycol int) (covFromQR *CovMatrix, means []float64, N float64) {

	if wycol < 0 {
		panic("wycol must be >= 0")
	}
	if wycol >= qr.Nyvar {
		panic("wycol must be < qr.Nyvar")
	}

	nc := qr.Ncol + 1

	sup := SupplementR(&qr.UpperTriMatrix, qr.Rhs[wycol], qr.Sserr[wycol], NoRecycle)

	nr := qr.AccumWeightSum
	N = nr

	covFromQR = NewCovMatrix(nc)
	xmeans := qr.XStats.Mean()
	ymeans := qr.YStats.Mean()

	// fill the mean of the 1s column
	means = append([]float64{1}, xmeans...)
	// just do one y-target column, wycol, to avoid colinearity issues
	// and to try and preserve the accuracy of the y-predictions as
	// much as possible.

	means = append(means, ymeans[wycol])

	//fmt.Printf("\n QR2Cov(): nc = %d  len(xmeans)=%d  len(ymeans)=%d  len(means)=%d\n", nc, len(xmeans), len(ymeans), len(means))
	//fmt.Printf("\n QR2Cov(): means = %#v\n", means)

	if len(means) != nc {
		panic(fmt.Sprintf("len(means)==%d must be equal to nc==%d\n", len(means), nc))
	}

	fNr := float64(nr)
	fNrMinus1 := fNr - 1
	var u2, raw float64

	//R: umat[i,j] = u[i]*u[j]*nr
	//R: cv=(t(R) %*% R - umat)/(nr-1)
	for i := 0; i < nc; i++ {
		for j := 0; j < nc; j++ {
			u2 = means[i] * means[j] * fNr

			raw = sup.DotColumns(i, j) - u2

			covFromQR.A[nc*i+j] = raw / fNrMinus1
		}
	}

	return
}

// and the inverse
/*
# R code:
# there are two square roots, set neg.root
# to true if you want the negative one
#
xcov2xR=function(cv, nc, nr, xcolmean, ymean, neg.root=FALSE) {

  if (length(ymean) != 1) {
    stop("ymean must be of length 1")
  }

  nc = nc+1
  umat=matrix(ncol=nc, nrow=nc,data=0)
  u = c(xcolmean, ymean)

  for (i in 1:nc) {
    for (j in 1:nc) {
      umat[i,j] = u[i]*u[j]*nr
    }
  }

  rr = ((nr-1)*cv) + umat

  ch = chol(rr)
  if (ch[1,1]>0 && neg.root) {
    return(-ch)
  }
  ch
}
*/
const SetQRCompressed = true

//
// *Always* assume exactly one y now, and thus one ymean, now.
//
// Doesn't do the DeSupplementR() for you, do that after Cov2QR().
//
func (m *UpperTriMatrix) Cov2QR(cov *CovMatrix, nr float64, xmean []float64, ymean float64, bCompress bool) {

	Numxvar := len(xmean) - 1
	if xmean[0] != 1 {
		panic("xmean[0] should be 1, for the const column")
	}

	nc := Numxvar + 2 // + 1 is for the constant, +1 for y

	if nc != m.Ncol {
		panic(fmt.Sprintf("m.Ncol is the wrong size(%d). We expected %d since len(xmean) was %d\n", m.Ncol, nc, len(xmean)))
	}

	means := append(xmean, ymean)

	S := NewSquareMatrix(nc)
	fNr := nr
	fNrMinus1 := fNr - 1

	var u2, raw float64

	for i := 0; i < nc; i++ {
		for j := i; j < nc; j++ {
			u2 = means[i] * means[j] * fNr
			raw = cov.At(i, j)*fNrMinus1 + u2
			S.Set(i, j, raw)
		}
	}

	// with CHOL_OVERWRITE_VARIANCES, chol will be the same as S.
	// This is fine because S is a temp (square for convenience)
	// matrix anyway, and we copy into m at the end.
	//
	chol, err, singList := S.FactorToCholeskyLower(CHOL_OVERWRITE_VARIANCES)
	//chol, err, singList := S.FactorToCholeskyLower(CHOL_RETURN_FRESH_MATRIX)

	if err != nil {
		panic(err)
	}
	if len(singList) != 0 {
		fmt.Printf("singList = %v\n", singList)
	}

	// pick the positive root for the diagonal
	for i := 0; i < chol.Ncol; i++ {
		if chol.At(i, i) < 0 {
			fmt.Printf("\n after chol, found negative root on the diagonal, flipping signs on column %d\n", i)
			chol.ColMultiply(i, -1)
		}
	}

	// put into m
	//fmt.Printf("m.Ncol = %d\n", m.Ncol)
	//fmt.Printf("chol.Ncol = %d\n", chol.Ncol)

	// do diagonal first, then the rest, to be sure that the multipliers
	// on the diagonal are set before everything else when bCompress is true.
	for i := 0; i < nc; i++ {
		for j := i; j < nc; j++ {
			// chol result is in the lower triangle, so j,i not i,j
			raw = chol.At(j, i)
			//fmt.Printf("raw = %v  i=%d  j=%d\n", raw, i, j)
			if bCompress {
				m.SetCompressed(i, j, raw)
			} else {
				m.Set(i, j, raw)
			}
		}
	}
}
