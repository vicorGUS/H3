#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "DMC.h"
#include "DMC_6dim.c"


int dim6(){
    
    int I = 20000;
    int i0 = 10000;
    int Nmax = 2000;
    int* N = malloc((I+1) * sizeof(int));
    
    double* pos[Nmax][6];
    
    for (int j=0; j<Nmax; j++) {
        pos[j] = (double*)malloc(6 * sizeof(double));
    }
    
    int count = 0;
    for (int j = 0; j < Nmax; j++)
            for (int i = 0; i < (I+1); i++)
                pos[j][i] = ++count;
    
    double* ET = malloc((I+1) * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));
    double* ET_avg_acc = malloc((I+1) * sizeof(double));
    
    double dtau = 0.01;
    double min = -5.;
    double max = 5.;
    double gamma = 0.5;
    
    srand((unsigned int)time(NULL));
    
    DMC(I, i0, N, Nmax, pos, ET, ET_avg, ET_avg_acc, dtau, min, max, gamma);
    
    FILE *fp;
    fp = fopen("DMC.csv", "wb");
    for (int i=0; i<I; i++) {
        fprintf(fp, "%f,%f,%f\n", ET[i], ET_avg[i], ET_avg_acc[i]);
    }
    fclose(fp);

    FILE *f2;
    f2 = fopen("x.csv", "wb");
    for (int j = 0; j < N[I]; j++){
        fprintf(f2, "%f,%f,%f,%f,%f,%f\n", pos[j][0],pos[j][1], pos[j][2], pos[j][3], pos[j][4], pos[j][5]);
    }
    fclose(f2);
    return 0;
}