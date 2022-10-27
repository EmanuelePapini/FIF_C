#include <fftw3.h>

double* realFFT(double *, int);

fftwl_complex* fft_dir(double *, int);

double* fft_inv(fftwl_complex *, int);
