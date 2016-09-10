package lsq

// dot product of two columns in the QR of lsq.

func (u *UpperTriMatrix) DotColumns(i int, j int) float64 {
	var sum float64
	n := u.Ncol
	for k := 0; k < n; k++ {
		sum += u.AtDecompressed(k, i) * u.AtDecompressed(k, j)
	}
	return sum
}
