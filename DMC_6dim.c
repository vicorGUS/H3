#include "DMC.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int pos_idx(int i, int j, int k, int I, int N){
    if(i>=I || j> N || k>6){
        return(-1);
    }
    return (k * (I+1) * N) + (j * I) + i;
}

void DMC_6dim(int I, int i0, int* N, int Nmax, double* pos, double* ET, double* ET_avg, double* ET_avg_acc, double dtau, double gamma){

    int m = 0;
    double* V = malloc(Nmax * sizeof(double));
    double* W = malloc(Nmax * sizeof(double));
    double r1 = 0.7 + rand() / RAND_MAX;
    double r2 = 0.7 + rand() / RAND_MAX;
    double t1 = acos(2. * rand() / RAND_MAX - 1);
    double t2 = acos(2. * rand() / RAND_MAX - 1);
    double p1 = 2. * 3.14 * rand() / RAND_MAX;
    double p2 = 2. * 3.14 * rand() / RAND_MAX;
    
    for (int i=0; i<(I+1); i++) {
        N[i] = 0;
    }
    
    N[0] = 1000;
    ET[0] = -3.0;
    

    for (int j=0; j<N[0]; j++){
        pos[pos_idx(0,j,0,I,N[I])] = (double) r1 * sin(t1) * cos(p1); //x1
        pos[pos_idx(0,j,1,I,N[I])] = (double) r2 * sin(t1) * sin(p1); //y1
        pos[pos_idx(0,j,2,I,N[I])] = (double) r1 * cos(t1);           //z1
        pos[pos_idx(0,j,3,I,N[I])] = (double) r2 * sin(t2) * cos(p2); //x2
        pos[pos_idx(0,j,4,I,N[I])] = (double) r2 * sin(t2) * sin(p2); //y2
        pos[pos_idx(0,j,5,I,N[I])] = (double) r2 * cos(t2);           //z2 
    }
    
    for (int j=N[0]; j<Nmax; j++ ) {
        for (int k=0; k<Nmax; k++ ){
            pos[pos_idx(0,j,k,I,N[I])] = 0.0;
        }
    }
    
    printf("START: N: %d, E: %f\n", N[0], ET[0]);
    
    for (int i=0; i<I; i++){
        ET_avg[i] = 0;
        if (N[i]<Nmax) {
            for (int j=0; j<N[i]; j++){
                for (int k = 0; k<6; k++){
                    pos[pos_idx(i,j,k,I,N[I])] += sqrt(dtau) * randn(0., 1.);
                    V[j] = (1./2) * pow(1 - exp(- pos[pos_idx(i,j,k,I,N[I])]), 2);
                        
                    W[j] = exp(- (V[j] - ET[i]) * dtau);
                    m = (int)(W[j] + ((double) rand() / (RAND_MAX)));
                        
                    N[i+1] += m/6;
                        
                    if (N[i+1] < Nmax){
                        for (int k=(N[i+1]-m); k<N[i+1]; k++){
                            pos[pos_idx(i+1,j,k,I,N[I])] = pos[pos_idx(i,j,k,I,N[I])];
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
            ET_avg_acc[i] = ET_avg[i+1];
           
            printf("i: %d N: %d E: %f\n",i, N[i+1], ET_avg[i+1]);
        }
    }
}