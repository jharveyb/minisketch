# Math improvements

## Multipliers

https://en.wikipedia.org/wiki/Toom%E2%80%93Cook_multiplication

## FFT

https://cr.yp.to/f2mult.html

LCH additive FFT

https://blog.lambdaclass.com/additive-fft-background/


# Existing code as examples

## FFT

additive DFT in GF(2), with SSE2 and PCLMULQDQ, both Cantor and Mateer-Gao
Maybe also AVX2, BMI1, and POPCNT?

https://github.com/kunzjacq/Additive_DFTs

INRIA library for arithmetic in GF(2)

https://gitlab.inria.fr/gf2x/gf2x

Sage + C++ implementation

https://github.com/mtbadakhshan/additive-fft

FFT over finite fields

https://github.com/lambdaclass/lambdaworks/blob/main/crates/math/src/fft/README.md

Circuit generator for Frobenius Additive FFT

https://github.com/fast-crypto-lab/Frobenius_AFFT
