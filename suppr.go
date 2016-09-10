package lsq

// R + [Rhs, Sserr]' --SupplementR()--> SupplementedR -> xtraCov matrix.
// then xtraCov -> SupplmentedR ---DeSupplementR()---> R + [Rhs, Sserr]'.

// generate a new UpperTriMatrix without the last column and row of m.
// The former last column is returned in Rhs and Sserr, with Sserr holding
// m's lower-right corner element, and Rhs having all the rest. Rhs is the
// last column of m without the bottom entry.
//
// keep it all compressed. Generates all new data, no sharing of pointers.
// If recycle is provided and the right size (recycle == m.Ncol -1) then
// we will right to recycle and return a pointer to recycle in desup.

var NoRecycle *UpperTriMatrix = nil

func (m *UpperTriMatrix) DeSupplementR(recycle *UpperTriMatrix) (desup *UpperTriMatrix, Rhs []float64, Sserr float64) {
	ncol := m.Ncol
	m1 := ncol - 1

	if recycle != nil && recycle.Ncol == m1 {
		recycle.ZeroOut()
		desup = recycle
	} else {
		//fmt.Printf("warning: no recycling in DeSupplementR\n")
		desup = NewUpperTriMatrix(m1)
	}
	Rhs = make([]float64, m1)
	Sserr = m.At(m1, m1)

	// R and D
	for i := 0; i < m1; i++ {
		for j := i; j < m1; j++ {
			desup.Set(i, j, m.At(i, j))
		}
	}

	// Rhs
	for i := 0; i < m1; i++ {
		Rhs[i] = m.At(i, m1)
	}

	return desup, Rhs, Sserr
}

// get a super-sized, supplemented upper triangular R matrix back.
// will try to recycle and return a pointer to recycle if provided, to
// avoid allocation. recycle can be nil to force allocation.
func SupplementR(r *UpperTriMatrix, Rhs []float64, Sserr float64, recycle *UpperTriMatrix) (sup *UpperTriMatrix) {
	ncol := r.Ncol
	p1 := ncol + 1

	if recycle != nil && recycle.Ncol == p1 {
		recycle.ZeroOut()
		sup = recycle
	} else {
		// recycle not provided, or wrong size. So allocate.
		//fmt.Printf("warning: no recycling in SupplementR\n")
		sup = NewUpperTriMatrix(p1)
	}
	var i, j int
	for i = 0; i < ncol; i++ {
		for j = i; j < ncol; j++ {
			sup.Set(i, j, r.At(i, j))
		}
	}

	j = ncol
	for i := 0; i < ncol; i++ {
		sup.Set(i, j, Rhs[i])
	}
	sup.Set(i, i, Sserr)

	return sup
}
