/* 
    
    Collection of FFT functions needed by Fif.c
    
    This file is part of the C implementation of Fast Iterative Filtering (FIF)
    

    Authors: Igor Bertello, Emanuele Papini
    Affiliation(s): IAPS - INAF

    Dependencies: FFTW3

*/

#include "FFT.h"
#include <stdlib.h>
//#include <fftw3.h>


double* realFFT(double * f, int N){

  fftwl_complex *in=NULL, *out=NULL;
  fftwl_plan p;
 
  double *fout;
  bool foutB;

  
  fout=(double *)malloc(sizeof(double)*N);
  foutB=1;
  
  
  in=(fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex)*N);
  out=(fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * N);
  p   = fftwl_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fflush(NULL);

  // compila il vettore in
  for(int i=0; i<N; i++){
     in[i][0]=f[i];
     in[i][1]=0;
  }

  fftwl_execute(p); 

  for(int i=0; i<N; i++){
     fout[i]=out[i][0];
  }
  
  fftwl_destroy_plan(p);
  fftwl_free(in); fftwl_free(out); // */

return fout;
}




fftwl_complex* fft_dir(double *f, int N){

  double *fout;
     
  fftwl_complex *in, *out;
  fftwl_plan p;

  in  = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N);
  out = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N);
  p   = fftwl_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  // compila il vettore i
  for(int i=0; i<N; i++){
     in[i][0]=f[i];
     in[i][1]=0;
  }
  
  fftwl_execute(p); 

  fftwl_destroy_plan(p);
  fftwl_free(in); 

return out;

}


double* fft_inv(fftwl_complex *f, int N){

  double *fout;
     
  fout=(double *)malloc(sizeof(double)*N);


  fftwl_complex  *out;
  fftwl_plan p;
  
  out = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N);

  p   = fftwl_plan_dft_1d(N, f, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwl_execute(p); 
  // compila il vettore out
  for(int i=0; i<N; i++){
     fout[i]=out[i][0]/N;
  } 

  fftwl_destroy_plan(p);

  fftwl_free(out); 

return fout;

}
