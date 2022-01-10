#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "DMC.h"

int main(){
    
    int I = 85000;
    int i0 = 0;
    int Nmax = 4000;
    
//    double* x[Nmax];
//
//    for (int j=0; j<Nmax; j++) {
//        x[j] = (double*)malloc(I * sizeof(double));
//    }
//
//    int count = 0;
//    for (int j = 0; j < Nmax; j++)
//            for (int i = 0; i < (I+1); i++)
//                x[j][i] = ++count;
    
    double* ET = malloc((I+1) * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));

    double* pos = malloc(6 * (I+1) * Nmax * sizeof(double));
    int* N = malloc((I+1) * sizeof(int));
    
    double dtau = 0.01;
    double gamma = 0.5;
    
    srand((unsigned int)time(NULL));
    
    //DMC(I, i0, N, Nmax, x, ET, ET_avg, dtau, gamma);
    
    DMC_6dim(I, i0, N, Nmax, pos, ET, ET_avg, dtau, gamma);
    
//    FILE *fp;
//    fp = fopen("DMC.csv", "wb");
//    for (int i=0; i<I; i++) {
//        fprintf(fp, "%f, %f, %d\n", ET_avg[i], ET[i], N[i]);
//    }
//    fclose(fp);
//
//    FILE *f;
//    f = fopen("x.csv", "wb");
//    for (int j=0; j<N[I]; j++) {
//        fprintf(f, "%f\n", x[j][I]);
//    }
//    fclose(f);
    
    FILE * f;
    f = fopen("DMC_6dim.csv", "wb");
    for (int i=0; i<I; i++) {
        fprintf(f, "%f, %d\n", ET_avg[i], N[i]);
    }
    fclose(f);
    return 0;
}
