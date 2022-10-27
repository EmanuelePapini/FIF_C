/* 
    Main file for testing the C implementation of Fast Iterative Filtering (FIF)
    
    written based on  https://github.com/Acicone/FIF/blob/master/Example_v8.m

    Authors: Igor Bertello, Emanuele Papini
    Affiliation(s): IAPS - INAF, University of L'Aquila (Italy)

*/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include "Fif.h"
//#define pi 3.145926
#define lung 1000
#define maxIMF 10 

int main(int varn, char *varc[]){

  Fif_t IMFs;
  Fif_t *IMFappo;
  double f[lung];
  double r=1;
  double t,dt;
  int* prova;
  clock_t begin_T, end_T;
  double dT;
  double pi;
  int numIMF = maxIMF;
  char fname[13];
  dt = 0.001; 
  pi = acos(-1);
// create the signal to analise
  srand(time_t(NULL));
  for(int i=0; i<lung; i++){
    t = i*dt;
    f[i] = (2 *pow(t-0.5,2) +0.2)*sin(20*pi*t + 0.2*cos(40*pi*t)) +
           4*pow(t-0.5,2) + 1;
    //f[i] =5*sin(i*dt*pi/500)+5*sin(i*pi/200);//+r;
  }
sprintf(fname, "signal.dat");
FILE *fil = fopen(fname,"w");
for (int i=0; i<lung;i++){
    fprintf(fil,"%f\n",f[i]);
}
fclose(fil);

  begin_T=clock();
// RUN FIF
  IMFs=FIF_v2_1(f,lung,&numIMF);
  end_T=clock();
  dT=(float)(end_T-begin_T)/CLOCKS_PER_SEC;
  printf("dt=%8.6f\n", dT);
  printf("# IMFs extracted: %d (max %d)\n",numIMF,maxIMF);
  
  IMFappo = &IMFs;

  for (int j=0; j<numIMF; j++)
  {
      sprintf(fname, "IMF%d.dat",j+1);
      FILE *fil = fopen(fname,"w");
      for (int i=0; i<lung;i++){
          fprintf(fil,"%f\n",IMFappo->dati[i]);
      }
      fclose(fil);
      IMFappo = IMFappo->next;
  }
  {
      sprintf(fname, "residual.dat");
      FILE *fil = fopen(fname,"w");
      for (int i=0; i<lung;i++){
          fprintf(fil,"%f\n",IMFappo->dati[i]);
      }
      fclose(fil);
  }
return 0;
}



