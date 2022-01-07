#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "DMC.h"
#include "DMC.c"
#include "DMC_6dim.c"


int main(){
    /*
    int I = 20000;
    int i0 = 10000;
    int Nmax = 2000;
    int* N = malloc((I+1) * sizeof(int));
    
    double* x[Nmax];
    
    for (int j=0; j<Nmax; j++) {
        x[j] = (double*)malloc(I * sizeof(double));
    }*/

    int I = 200;
    int i0 = 10000;
    int Nmax = 10000;
    int* N = malloc((I+1) * sizeof(int));
    
    double* pos = malloc(6 * (I+1) * (I+1) * sizeof(double));
    double* ET = malloc((I+1) * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));
    double* ET_avg_acc = malloc((I+1) * sizeof(double));
    
    double dtau = 0.02;
    double min = -5.;
    double max = 5.;
    double gamma = 0.5;
    
    srand((unsigned int)time(NULL));
    
    // DMC(I, i0, N, Nmax, x, ET, ET_avg, ET_avg_acc, dtau, min, max, gamma);
    DMC_6dim(I, i0, N, Nmax, pos, ET, ET_avg, ET_avg_acc, dtau, min, max, gamma);
    
    FILE *fp;
    fp = fopen("DMC.csv", "wb");
    for (int i=0; i<I; i++) {
        fprintf(fp, "%f,%f,%f\n", ET[i], ET_avg[i], ET_avg_acc[i]);
    }
    fclose(fp);
    FILE *f2;
    f2 = fopen("x.csv", "wb");
    for (int j = 0; j < N[I]; j++){
        // fprintf(f2, "%f\n", pos[j][I]);
        fprintf(f2, "%f,%f,%f,%f,%f,%f\n", pos[pos_idx(I,j,0,I,N[I])],pos[pos_idx(I,j,1,I,N[I])], pos[pos_idx(I,j,2,I,N[I])], pos[pos_idx(I,j,3,I,N[I])], pos[pos_idx(I,j,4,I,N[I])], pos[pos_idx(I,j,5,I,N[I])]);
    }
    fclose(f2);
    return 0;
}