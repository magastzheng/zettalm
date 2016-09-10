package lsq

import (
	"fmt"
	"math"
)

// UpperTriMatrix:
// factor out the upper-right triangular
// matrix construct from MillerLSQ into its own
// reusable type, UpperTriMatrix.
type UpperTriMatrix struct {
	Ncol  int
	R_dim int       // len(R)
	R     []float64 // row-major, upper triangular without the diagonal
	D     []float64 // the diagonal elements

	Row_ptr []int // easy to get at each row in turn

	// 0-based Rowlen in R[] is just (Ncol - i) -1, not counting
	// the diagonal entry. e.g. so on the last row, where i=Ncol-1,
	// we have 0 entries in R for that row.

}

// At: i and j are 0-based request for matrix element m[i,j]
func (m *UpperTriMatrix) At(i, j int) float64 {
	if i == j {
		return m.D[i]
	}
	if j < i {
		return 0
	}
	return m.R[m.Row_ptr[i]+j-i-1]
}

// Set assigns val to m[i,j]
func (m *UpperTriMatrix) Set(i, j int, val float64) {
	if i == j {
		m.D[i] = val
		return
	}
	if j < i {
		panic("trying to Set() on lower part of UpperTrimatrix, which must be all 0s.")
	}
	m.R[m.Row_ptr[i]+j-i-1] = val
}

func (m *UpperTriMatrix) SetCompressed(i, j int, val float64) {
	if i == j {
		m.D[i] = val * val
		return
	}

	//fmt.Printf("\n in SetCompressed(i=%d, j=%d, val=%v), m.D[i]=%v, val/m.D[i]=%v\n", i, j, val, m.D[i], val/m.D[i])

	if j < i {
		panic("trying to SetCompressed() on lower part of UpperTrimatrix, which must be all 0s.")
	}
	m.R[m.Row_ptr[i]+j-i-1] = val / math.Sqrt(m.D[i])
}

func (m *UpperTriMatrix) RowMultiply(irow int, scalar float64) {
	if scalar == 1.0 {
		return
	}
	n := m.Ncol
	m.D[irow] *= scalar
	numj := n - irow - 1
	begin := m.Row_ptr[irow]
	for j := 0; j < numj; j++ {
		m.R[begin+j] *= scalar
	}
}

func (m *UpperTriMatrix) AtDecompressed(i, j int) float64 {
	di := math.Sqrt(m.D[i])
	if i == j {
		return di
	}
	if j < i {
		return 0
	}
	return di * m.R[m.Row_ptr[i]+j-i-1]
}

func NewUpperTriMatrix(ncol int) *UpperTriMatrix {

	m := &UpperTriMatrix{
		Ncol: ncol,
	}
	m.R_dim = m.Ncol * (m.Ncol - 1) / 2

	// the upper triangle, without diagonal. For precision preservation,
	// columns are multiplied up by their D[column] entry
	// when displayed, but otherwise kept factored, so the value
	// stored in m.R[pos for column j] is actually 1/D[j] of its
	// 'actual' value. This is a classic technique in numerical
	// code to preserve accuracy.
	m.R = make([]float64, m.R_dim)

	// the Diagonal
	m.D = make([]float64, m.Ncol)
	//	for i := range m.D {
	//		m.D[i] = 1.0
	//	}

	m.Row_ptr = make([]int, m.Ncol)
	// m.Row_ptr[i] is the position of element R(i,i+1) in array m.R.
	// Since the diagonal is implicit 1's at R(i,i), this is the
	// first relevant/interesting entry.

	// Row_ptr entries are accurately 0-based or 1-based depending on Adj
	m.Row_ptr[1-Adj] = 1 - Adj
	for i := 2; i <= m.Ncol-1; i++ {
		m.Row_ptr[i-Adj] = m.Row_ptr[(i-1)-Adj] + m.Ncol - i + 1
	}

	// the m.Row_ptr[m.Ncol-Adj] entry should never be used; it doesn't exist since
	// there is only a 1 on the diagnonal in the R(Ncol, Ncol) 1-based position.
	// so here we make sure this last entry is out-of-bounds for m.R
	m.Row_ptr[m.Ncol-Adj] = m.R_dim + 1

	return m
}

func (m *UpperTriMatrix) String() string {
	return m.UpperTriToString("\n", StdFormatE)
}

func (m *UpperTriMatrix) String6() string {
	return m.UpperTriToString("\n", StdFormat6dec)
}

// linesep defaults to "\n" unless supplied differently
func (m *UpperTriMatrix) UpperTriToString(linesep string, formatFunc func(float64) string) string {
	ncol := m.Ncol

	if linesep == "" {
		linesep = "\n"
	}
	var i, j, pos int
	// R-code formatted, to make it easy to bring into R.
	s := fmt.Sprintf("matrix(ncol=%d, byrow=TRUE, data=c(%s", ncol, linesep)
	pos = 0
	for i = 1; i <= ncol; i++ {

		for j = 1; j <= ncol; j++ {
			if j < i {
				// for formatFunc puts a space before it adds the number.
				s += formatFunc(0.0) + ","
			} else if j == i {
				s += formatFunc(m.D[i-Adj])
				if i < ncol {
					s += ","
				}
			} else {
				// j > i, upper triangle
				s += formatFunc(m.R[pos]) + ","
				pos++
			}
		}

		if i < ncol {
			s += " " + linesep
		} else {
			s += "))"
		}
	}

	return s
}

func (m *UpperTriMatrix) DecodedUpperTriToString(linesep string, formatFunc func(float64) string) string {
	ncol := m.Ncol

	if linesep == "" {
		linesep = "\n"
	}
	var i, j, pos int
	// R-code formatted, to make it easy to bring into R.
	s := fmt.Sprintf("matrix(ncol=%d, byrow=TRUE, data=c(%s", ncol, linesep)
	pos = 0
	for i = 1; i <= ncol; i++ {
		di := math.Sqrt(m.D[i-Adj])
		for j = 1; j <= ncol; j++ {
			if j < i {
				// for formatFunc puts a space before it adds the number.
				s += formatFunc(0.0) + ","
			} else if j == i {
				s += formatFunc(di)
				if i < ncol {
					s += ","
				}
			} else {
				// j > i, upper triangle
				s += formatFunc(di*m.R[pos]) + ","
				pos++
			}
		}

		if i < ncol {
			s += " " + linesep
		} else {
			s += "))"
		}
	}

	return s
}

func (m *UpperTriMatrix) ZeroOut() {
	zero_out(m.D)
	zero_out(m.R)
}

// CholeksyFactor:
// PRE: matrix A is a symmetric n x n matrix.
// POST: we return in m the square matrix containing
// the *upper* triangular cholesky decomposition
//
// This routine is used for test purposes; obviously
// the MillerLSQ routines themselves provide
// the same functionality; and handle infinite data.
//
func (m *UpperTriMatrix) CholeskyFactor(A []float64) (err error) {
	var n int = m.Ncol
	if len(A) != n*n {
		panic(fmt.Sprintf("A input needs to be n x n, where n = %d, but len(A) was %d", n, len(A)))
	}

	U := make([]float64, n*n)

	var s float64
	for i := 0; i < n; i++ {
		for j := 0; j < (i + 1); j++ {
			s = 0
			for k := 0; k < j; k++ {
				//fmt.Printf("U[k=%d *n+ i=%d]=%v\n", k, i, U[k*n+i])
				//fmt.Printf("U[k=%d *n+ j=%d]=%v\n", k, j, U[k*n+j])

				s += U[k*n+i] * U[k*n+j]
				//fmt.Printf("at i,j,k = %d, %d, %d,  s = %v\n", i, j, k, s)
			}
			if i == j {
				tmp := A[i*n+i] - s
				if tmp < 0.0 {
					err = fmt.Errorf("Cholesky failed: input matrix A had negative on diagnonal; not semi-positive-definite, at i = %d", i)
					return err
				}
				//fmt.Printf("CholeksyFactor() diagnonal[%d] as sqrt of (A[%d*n+%d]=%v - %v) = %v\n", i, i, i, A[i*n+i], s, tmp)
				U[j*n+i] = math.Sqrt(tmp)
				//fmt.Printf("CholeskyFactor() j=%d, i=%d, U[j*n+i] = %v\n", j, i, U[j*n+i])

			} else {
				U[j*n+i] = (1.0 / U[j*n+j] * (A[j*n+i] - s))
			}
			//fmt.Printf("chol.go: at end of i=%d  j=%d  we have U = :\n%v\n", i, j, U)
		}
	}

	m.ImportFromSymmetricSquare(U)
	return nil
}

func (m *UpperTriMatrix) ImportFromSymmetricSquare(A []float64) {
	n := m.Ncol
	if len(A) != n*n {
		panic(fmt.Sprintf("non-conforming shape error: A input needs to be m.Ncol x m.Ncol == %d, where m.Ncol = %d; but len(A) was %d", n*n, n, len(A)))
	}

	nextr := 0
	for i := 0; i < n; i++ {
		m.D[i] = A[i*n+i]
		for j := i + 1; j < n; j++ {
			m.R[nextr] = A[i*n+j]
			nextr++
		}
	}
}

// decompress undoes the row-multiplier D that lsq uses.
func (m *UpperTriMatrix) ExtractSquareR(decompress bool) (s *SquareMatrix) {
	s = NewSquareMatrix(m.Ncol)

	ncol := m.Ncol

	var di, denorm float64
	var i, j, pos int
	for i = 0; i < ncol; i++ {
		if decompress {
			di = math.Sqrt(m.D[i])
			denorm = di
		} else {
			di = m.D[i]
			denorm = 1
		}
		for j = 0; j < ncol; j++ {
			if j < i {
				// lower triangle, all zeros
				s.A[i*ncol+j] = 0
			} else if j == i {
				s.A[i*ncol+i] = di
			} else {
				// j > i, upper triangle
				s.A[i*ncol+j] = denorm * m.R[pos]
				pos++
			}
		}
	}

	return
}
