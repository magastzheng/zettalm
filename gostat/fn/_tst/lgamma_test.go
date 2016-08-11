package fn

import (
	"fmt"
	"testing"
)

// Test of ln(Gamma) against R:lgamma
func Test_LnGamma(t *testing.T) {
	const delta = 1e-4

	fmt.Println("Testing LnGamma #1")
	a := 3.15
	x := LnGamma(a)
	y := 0.8359236

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

	fmt.Println("Testing LnGamma #2")
	a = 432.123
	x = LnGamma(a)
	y = 2188.191

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

	fmt.Println("Testing LnGamma #3")
	a = 0.0000000675
	x = LnGamma(a)
	y = 16.51114

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

	fmt.Println("Testing LnGamma #4")
	a = 3.15e-19
	x = LnGamma(a)
	y = 42.60171

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

	fmt.Println("Testing LnGamma #5")
	a = 3.15e-99
	x = LnGamma(a)
	y = 226.8085

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

	fmt.Println("Testing LnGamma #6")
	a = -12.33
	x = LnGamma(a)
	y = -19.53042

	if abs(x-y) > delta {
		fmt.Println("failed: ", x, y)
		t.Error()
	}

}
