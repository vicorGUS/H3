#include "DMC.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

void DMC(int I, int i0, int Nmax, double* x[Nmax], double* ET, double* E_avg, double dtau, double min, double max, double gamma){

    int m = 0;
    
    int* N = malloc((I+1) * sizeof(int));
    double* V = malloc(Nmax * sizeof(double));
    double* W = malloc(Nmax * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));
    
    for (int i=0; i<(I+1); i++) {
        N[i] = 0;
    }
    
    N[0] = 200;
    ET[0] = 0.5;

    for (int j=0; j<N[0]; j++){
        x[j][0] = -5. + (0.05 * j);
    }
    
    for ( int j=N[0]; j<Nmax; j++ ) {
          x[j][0] = 0.0;
    }
    
    printf("START: N: %d, E: %f\n", N[0], ET[0]);
    
    for (int i=0; i<I; i++){
        ET_avg[i] = 0;
        if (N[i]<Nmax) {
            for (int j=0; j<N[i]; j++){
                
                x[j][i] += sqrt(dtau) * randn(0., 1.);
                
                if (x[j][i]>5.0) {
                    x[j][i] = 5.0;
                }
                if (x[j][i]<-5.0) {
                    x[j][i] = -5.0;
                }
                    
                if (fabs(x[j][i])<5.0) {
                        
                    V[j] = (1./2) * pow(1 - exp(- x[j][i]), 2);
                        
                    W[j] = exp(- (V[j] - ET[i]) * dtau);
                    m = (int)(W[j] + ((double) rand() / (RAND_MAX)));
                        
                    N[i+1] += m;
                        
                    if (N[i+1] < Nmax){
                            
                        for (int k=(N[i+1]-m); k<N[i+1]; k++){
                            x[k][i+1] = x[j][i];
                            }
                        }
                    }
                }

            for (int i´=(i0+1); i´<i; i´++) {
                ET_avg[i] += 1. / (i-i0) * ET[i´];
            }
            
            ET[i+1] = ET_avg[i] - gamma * log((double)N[i] / N[0]);
            
            if ((i-i0+1) != 0) {
                ET_avg[i+1] = 1 / (i - i0 + 1.) * ET[i+1] + (i - i0) / (i - i0 + 1.) * ET_avg[i];
            }
           
            printf("i: %d N: %d E: %f\n",i, N[i+1], ET_avg[i+1]);
        }
    }
    
}
