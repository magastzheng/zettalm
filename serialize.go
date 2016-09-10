package lsq

import (
	"bytes"
	"encoding/gob"
	"fmt"
	"log"
	"os"
)

func (lsq *MillerLSQ) GobEncode(outFileName string) {

	gob.Register(lsq)
	gob.Register(SdTracker{})

	var network bytes.Buffer        // Stand-in for a network connection
	enc := gob.NewEncoder(&network) // Will write to network.

	// Encode (send) some values.
	err := enc.Encode(lsq)
	if err != nil {
		log.Fatal("encode error:", err)
	}

	file, err := os.Create(outFileName)
	if err != nil {
		panic(fmt.Sprintf("problem creating lsq outfile '%s': %s", outFileName, err))
	}
	defer file.Close()
	drainable := network
	drainable.WriteTo(file)
}

func ReadLSQGobFile(fname string) *MillerLSQ {

	f, err := os.Open(fname)
	if err != nil {
		panic(err)
	}

	// Decode (receive) and print the values.
	dec := gob.NewDecoder(f) // Will read from f

	var lsq MillerLSQ
	err = dec.Decode(&lsq)
	if err != nil {
		log.Fatal("decode error 1:", err)
	}

	return &lsq
}
