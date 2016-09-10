package lsq

import (
	"errors"
	"fmt"
	"math"
)

// merge 2 covariance matrices

type SquareMatrix struct {
	Ncol     int
	A        []float64
	Colnames []string

	// D is the diagonal of the cholesky factorization, for when we write the
	// cholesky to the lower left triangle of A. D is left unallocated
	// until a cholesky is generated.
	D []float64

	// ShowD is set as side-effect of CHOL_OVERWRITE_LOWER_AND_D. This makes
	// String()-ification of the matrix will show D, but the underlying
	// diagonal in A is still there. Flip ShowD to false to see the A
	// contents (the Variances) again.
	//
	ShowD bool // when String() displays, put the D over the diagonal
}

type CovMatrix struct {
	SquareMatrix
}
type CholeskyMatrix struct {
	SquareMatrix
}

func NewCovMatrix(ncol int) *CovMatrix {
	return &CovMatrix{*NewSquareMatrix(ncol)}
}

func NewSquareMatrix(ncol int) *SquareMatrix {
	return &SquareMatrix{
		Ncol: ncol,
		A:    make([]float64, ncol*ncol),
	}
}

//
// Compute a Covariance matrix from a DataFrame (but without remove the means).
// Essentially this is a Matrix multiply: return t(X) %*% X,
// adding an implicit constant 1 column as the first column on the left.
//
// if addIntercept, compute XtX as if there was a column of 1's before
//  all other columns in x.
//
const RmMean = true
const AddIntercept = true

func XtX(x *DataFrame, rmMean bool, addIntercept bool) (s *SquareMatrix, meanvec []float64, nrow int) {

	off := 0 // offset for intercept
	if addIntercept {
		off = 1
	}

	scol := x.Ncol + off
	s = NewSquareMatrix(scol)
	if addIntercept {
		s.Colnames = append([]string{"const"}, x.Colnames...)
	} else {
		s.Colnames = x.Colnames
	}

	tracker := NewSdTracker(x.Ncol)
	for i := range x.Rows {
		tracker.AddObs(x.Rows[i], 1)
	}

	// mean stays zero if !rmMean
	mean := tracker.Mean()
	//sd := tracker.Sd()

	nr := x.Nrow
	nc := x.Ncol
	var sum float64

	if addIntercept {
		if rmMean {
			// leave zeros on first column and first row
		} else {
			// compute the products that touch the 0-th column
			s.A[0] = float64(nr)
			for j := 0; j < nc; j++ {
				sum = 0
				for k := 0; k < nr; k++ {
					sum += x.Rows[k][j] // * 1, the 1 being the implicit first column we add to x.
				}
				s.A[j+off] = sum        // s.A[0,j+1]
				s.A[scol*(j+off)] = sum // s.A[j+1, 0]
			}
		}
	} else {
		off = 0
	}

	// and the rest, those dot-products that don't touch the implicit 0-th column of 1's.
	var meani, meanj float64
	for i := 0; i < nc; i++ {
		for j := i; j < nc; j++ {

			// dot product between columns of x
			sum = 0
			meani = mean[i]
			meanj = mean[j]
			if !rmMean {
				meani = 0
				meanj = 0
			}
			for k := 0; k < nr; k++ {
				sum += (x.Rows[k][i] - meani) * (x.Rows[k][j] - meanj)
			}
			s.A[scol*(i+off)+(j+off)] = sum
			s.A[scol*(j+off)+(i+off)] = sum
		}
	}
	return s, mean, nr
}

func (m *SquareMatrix) String() string {

	//formatFunc := StdFormatE
	formatFunc := StdFormat6dec
	linesep := "\n"
	ncol := m.Ncol

	var i, j int
	// R-code formatted, to make it easy to bring into R.
	s := fmt.Sprintf("matrix(ncol=%d, nrow=%d, byrow=TRUE, data=c(%s", ncol, ncol, linesep)
	for i = 0; i < ncol; i++ {

		for j = 0; j < ncol; j++ {
			if j < i {
				// for formatFunc puts a space before it adds the number.
				s += formatFunc(m.A[ncol*i+j]) + ","
			} else if j == i {
				if m.ShowD {
					s += formatFunc(m.D[i])
				} else {
					s += formatFunc(m.A[ncol*i+j])
				}

				if i < ncol-1 {
					s += ","
				}
			} else {
				// j > i, upper triangle
				s += formatFunc(m.A[ncol*i+j]) + ","
			}
		}

		if i < ncol-1 {
			s += " " + linesep
		} else {
			s += "))"
		}
	}

	return s
}

//
// Cov : compute covariance matrix of X.
//
func Cov(x *DataFrame) (s *CovMatrix, meanvec []float64, nrow int) {
	xtx, mv, nr := XtX(x, RmMean, !AddIntercept)
	xtx.DivideBy(float64(x.Nrow - 1))
	return &CovMatrix{*xtx}, mv, nr
}

func CovAddConstKeepMean(x *DataFrame) (s *CovMatrix, meanvec []float64, nrow int) {
	xtx, mv, nr := XtX(x, !RmMean, AddIntercept)
	xtx.DivideBy(float64(x.Nrow - 1))
	return &CovMatrix{*xtx}, mv, nr
}

func CovKeepMean(x *DataFrame) (s *CovMatrix, meanvec []float64, nrow int) {
	xtx, mv, nr := XtX(x, !RmMean, !AddIntercept)
	xtx.DivideBy(float64(x.Nrow - 1))
	return &CovMatrix{*xtx}, mv, nr
}

func (m *SquareMatrix) DivideBy(d float64) {
	for i := range m.A {
		m.A[i] /= d
	}
}

type CholOpt int

const (
	CHOL_OVERWRITE_LOWER_AND_D = 0
	CHOL_RETURN_FRESH_MATRIX   = 1
	CHOL_OVERWRITE_VARIANCES   = 2
)

// if opt == CHOL_OVERWRITE_LOWER_AND_D, then we return a pointer to ourselves, m,
//  and the diagonal of the cholesky decomposition is written to m.D in this case,
//  so as to not distrub the variances on the covariance matrix.
//
// CHOL_OVERWRITE_VARIANCES means we ignore D and write everything to A, destroying
//  the variances there. If you want to keep them, you could copy them to D first.
//
func (m *SquareMatrix) FactorToCholeskyLower(opt CholOpt) (L *SquareMatrix, err error, singularList []int) {

	n := m.Ncol
	if opt == CHOL_RETURN_FRESH_MATRIX {
		L = NewSquareMatrix(n)
	} else {
		// even if CHOL_OVERWRITE_VARIANCES, still write to D
		// at first, so we can read the Variances as needed first.
		L = m
		if opt == CHOL_OVERWRITE_VARIANCES {
			L.ShowD = false
		} else {
			L.ShowD = true
		}

		// only allocate D if we have to.
		if len(m.D) == n {
			// already allocated. Reuse, by zero it first.
			zero_out(m.D)
		} else {
			m.D = make([]float64, n)
		}
	}

	var s, root float64
	for i := 0; i < n; i++ {
		for j := 0; j < (i + 1); j++ {
			s = 0
			for k := 0; k < j; k++ {
				s += L.A[i*n+k] * L.A[j*n+k] // k < j <= i, so these read only the lower triangle
			}
			if i == j {
				tmp := m.A[i*n+i] - s
				if tmp < 0.0 {
					// the singularList lets us know if there is more than one
					singularList = append(singularList, i)
					err = errors.New(fmt.Sprintf("Cholesky failed: input matrix A had negative on diagnonal; not semi-positive-definite, at i = %v", singularList))
					if opt == CHOL_RETURN_FRESH_MATRIX {
						L.A[i*n+j] = 0
					} else {
						m.D[i] = 0
					}
				} else {
					root = math.Sqrt(tmp)
					if opt == CHOL_RETURN_FRESH_MATRIX {
						L.A[i*n+j] = root
					} else {
						m.D[i] = root
					}
				}
			} else {
				if opt == CHOL_RETURN_FRESH_MATRIX {
					L.A[i*n+j] = ((m.A[j*n+i] - s) / L.A[j*n+j])
				} else {
					L.A[i*n+j] = ((m.A[j*n+i] - s) / m.D[j])
				}
			}
		}
	}

	if opt == CHOL_OVERWRITE_VARIANCES {
		// swap Variances and the Cholesky diagonal (their square roots)
		for i := 0; i < n; i++ {
			L.A[i*n+i], m.D[i] = m.D[i], L.A[i*n+i]
		}
	}

	return
}

func CombineCov(cov1 *CovMatrix, cov2 *CovMatrix, mean1 []float64, mean2 []float64, N1 float64, N2 float64) (merged *CovMatrix, meanc []float64, combinedN float64) {

	n := cov1.Ncol

	if cov1.Ncol != cov2.Ncol {
		panic("cov1 and cov2 must match dimensions")
	}
	if len(mean1) != n {
		panic(fmt.Sprintf("len(mean1)==%d must be the same as cov1.Ncol(%d)", len(mean1), cov1.Ncol))
	}
	if len(mean2) != n {
		panic(fmt.Sprintf("len(mean2)==%d must be the same as cov1.Ncol(%d)", len(mean2), cov1.Ncol))
	}

	merged = NewCovMatrix(n)
	combinedN = N1 + N2

	fN1 := float64(N1)
	fN2 := float64(N2)
	fcombinedN := float64(combinedN)

	// meanc
	meanc = make([]float64, n)
	for i := range mean1 {
		meanc[i] = (fN1*mean1[i] + fN2*mean2[i]) / fcombinedN
	}

	// merged
	var ans float64
	for u := 0; u < n; u++ {
		for v := 0; v < n; v++ {
			ans = (fN1-1)*cov1.A[n*u+v] + (fN2-1)*cov2.A[n*u+v] + (mean2[u]-mean1[u])*(mean2[v]-mean1[v])*fN1*fN2/fcombinedN
			merged.A[n*u+v] = ans / (fcombinedN - 1)
		}
	}

	return
}

func (m *SquareMatrix) ColumnDotProduct(i, j int) float64 {
	var sum float64
	n := m.Ncol
	for k := 0; k < n; k++ {
		sum += m.A[n*k+i] * m.A[n*k+j]
	}
	return sum
}

func (m *SquareMatrix) Transpose() {
	n := m.Ncol
	var tmp float64
	for i := 0; i < n; i++ {
		// only do one half, or you'll swap *back* everything on the 2nd half.
		for j := i + 1; j < n; j++ {
			// swap i,j and j,i
			tmp = m.A[n*i+j]
			m.A[n*i+j] = m.A[n*j+i]
			m.A[n*j+i] = tmp
		}
	}
}

func (m *SquareMatrix) At(i, j int) float64 {
	return m.A[m.Ncol*i+j]
}

func (m *SquareMatrix) Set(i, j int, val float64) {
	m.A[m.Ncol*i+j] = val
}

func (m *SquareMatrix) RowMultiply(irow int, scalar float64) {
	if scalar == 1.0 {
		return
	}
	n := m.Ncol
	for j := 0; j < n; j++ {
		m.A[n*irow+j] *= scalar
	}
}

func (m *SquareMatrix) ColMultiply(jcol int, scalar float64) {
	if scalar == 1.0 {
		return
	}
	n := m.Ncol
	for i := 0; i < n; i++ {
		m.A[n*i+jcol] *= scalar
	}
}
