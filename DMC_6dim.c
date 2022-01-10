#include "DMC.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int pos_idx(int i, int j, int k, int I, int N){
    return (i * 6 * N) + (j * 6) + k;
}

void DMC_6dim(int I, int i0, int* N, int Nmax, double* pos, double* ET, double* ET_avg, double dtau, double gamma){

    int m;
    double* V = malloc(Nmax * sizeof(double));
    double* W = malloc(Nmax * sizeof(double));
    double r1 = 0.7 + (double) rand() / RAND_MAX;
    double r2 = 0.7 + (double) rand() / RAND_MAX;
    double r12;
    double t1 = acos(2. * (double) rand() / RAND_MAX - 1);
    double t2 = acos(2. * (double) rand() / RAND_MAX - 1);
    double p1 = 2. * 3.14 * (double) rand() / RAND_MAX;
    double p2 = 2. * 3.14 *  (double) rand() / RAND_MAX;
    
    for (int i=1; i<(I+1); i++) {
        N[i] = 0;
    }
    
    N[0] = 1000;
    ET[0] = -3.0;

    for (int j=0; j<N[0]; j++){
        pos[pos_idx(0,j,0,I,Nmax)] = (double) r1 * sin(t1) * cos(p1); //x1
        pos[pos_idx(0,j,1,I,Nmax)] = (double) r1 * sin(t1) * sin(p1); //y1
        pos[pos_idx(0,j,2,I,Nmax)] = (double) r1 * cos(t1);           //z1
        pos[pos_idx(0,j,3,I,Nmax)] = (double) r2 * sin(t2) * cos(p2); //x2
        pos[pos_idx(0,j,4,I,Nmax)] = (double) r2 * sin(t2) * sin(p2); //y2
        pos[pos_idx(0,j,5,I,Nmax)] = (double) r2 * cos(t2);           //z2
    }
    
    for (int j=N[0]; j<Nmax; j++ ) {
        for (int k=0; k<6; k++ ){
            pos[pos_idx(0,j,k,I,Nmax)] = 0.0;
        }
    }
    
    printf("START: N: %d, E: %f\n", N[0], ET[0]);
    
    for (int i=0; i<I; i++){
        ET_avg[i] = 0;
        for (int j=0; j<N[i]; j++){
            for (int k = 0; k<6; k++){
                pos[pos_idx(i,j,k,I,Nmax)] += sqrt(dtau) * randn(0., 1.);
            }
                    
            r1 = sqrt(pow(pos[pos_idx(i,j,0,I,Nmax)], 2) + pow(pos[pos_idx(i,j,1,I,Nmax)], 2) + pow(pos[pos_idx(i,j,2,I,Nmax)], 2));
                    
            r2 = sqrt(pow(pos[pos_idx(i,j,3,I,Nmax)], 2) + pow(pos[pos_idx(i,j,4,I,Nmax)], 2) + pow(pos[pos_idx(i,j,5,I,Nmax)], 2));
            
            r12 = sqrt(pow(pos[pos_idx(i,j,3,I,Nmax)] - pos[pos_idx(i,j,0,I,Nmax)], 2) + pow(pos[pos_idx(i,j,4,I,Nmax)] - pos[pos_idx(i,j,1,I,Nmax)], 2) + pow(pos[pos_idx(i,j,5,I,Nmax)] - pos[pos_idx(i,j,2,I,Nmax)], 2));
                    
            V[j] = (-2./r1 - 2./r2 + 1./r12);
                    
            W[j] = exp( -(V[j] - ET[i]) * dtau);

            m = (int)(W[j] + ((double) rand() / (RAND_MAX)));

            N[i+1] += m;
            
            for (int k = 0; k<6; k++){
                for (int l=(N[i+1]-m); l<N[i+1]; l++){
                    pos[pos_idx(i+1,l,k,I,Nmax)] = pos[pos_idx(i,j,k,I,Nmax)];
                }
            }
            
        }
        
        if (N[i+1]>Nmax) {
            N[i+1] = (int) floor(N[i+1] / 2.);
            printf("%d\n", N[i+1]);
        }

        for (int i´=0; i´<i; i´++) {
            ET_avg[i] += (1. / i) * ET[i´];
        }
        
        ET_avg[0] = -3.0;
            
        ET[i+1] = ET_avg[i] - gamma * log((double)N[i] / N[0]);
        
        if ((i-i0+1) != 0) {
            ET_avg[i+1] = 1 / (i - i0 + 1.) * ET[i+1] + (i - i0) / (i - i0 + 1.) * ET_avg[i];
            }
        
//        printf("i: %d N: %d E: %f\n",i, N[i], ET_avg[i]);
    }
    printf("N: %d E: %f\n",N[I], ET_avg[I]);
}
