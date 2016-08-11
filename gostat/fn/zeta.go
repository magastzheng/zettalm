// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

import (
	"math"
)

// Riemann zeta function ζ(s)
func RiemannZeta(s float64) float64 {
	var (
		x float64
	)

	b := []float64{
		1.6666666666666666666666666666666666666666666666667e-1,
		-3.3333333333333333333333333333333333333333333333333e-2,
		2.3809523809523809523809523809523809523809523809524e-2,
		-3.3333333333333333333333333333333333333333333333333e-2,
		7.5757575757575757575757575757575757575757575757576e-2,
		-2.5311355311355311355311355311355311355311355311355e-1,
		1.1666666666666666666666666666666666666666666666667e+0,
		-7.0921568627450980392156862745098039215686274509804e+0,
	}

	switch s {
	case 0:
		x = -0.5
	case 1:
		x = math.NaN()
	case 1.5:
		x = 2.61237534868548834334856756792407163057080065240006340757332824881492776768827286099624386812631195238297
	case 2:
		x = 1.64493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870
	case 3:
		x = 1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915
	case 4:
		x = 1.08232323371113819151600369654116790277475095191872690768297621544412061618696884655690963594169991
	case 5:
		x = 1.03692775514336992633136548645703416805708091950191281197419267790380358978628148456004310655713333
	case 6:
		x = 1.01734306198444913971451792979092052790181749003285356184240866400433218290195789788277397793853517
	case 7:
		x = 1.00834927738192282683979754984979675959986356056523870641728313657160147831735573534609696891385132
	case 8:
		x = 1.00407735619794433937868523850865246525896079064985002032911020265258295257474881439528723037237197
	case 9:
		x = 1.00200839282608221441785276923241206048560585139488875654859661590978505339025839895039306912716958
	case 10:
		x = 1.00099457512781808533714595890031901700601953156447751725778899463629146515191295439704196861038565

	default:
		a := 12
		aa := float64(a)
		k := 8
		x = 0.0

		for i := 1; i < a; i++ {
			ii := float64(i)
			x += 1.0 / math.Pow(ii, s)
		}

		x += 1.0/((s-1)*math.Pow(aa, (s-1))) + 1/(2*math.Pow(aa, s))
		term := (s / 2) / math.Pow(aa, (s+1))
		x += term * b[0]

		for i := 1; i < k; i++ {
			ii := float64(i) + 1
			term *= (s + 2*ii - 2) * (s + 2*ii - 3) / (aa * aa * 2 * ii * (2*ii - 1))
			x += term * b[i]
		}
	}
	return x
}