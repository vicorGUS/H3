#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "DMC.h"


int main(){
    
    int I = 200000;
    int i0 = 10000;
    int Nmax = 2000;
    
    double* x[Nmax];
    
    for (int j=0; j<Nmax; j++) {
        x[j] = (double*)malloc(I * sizeof(double));
    }
    
    int count = 0;
    for (int j = 0; j < Nmax; j++)
            for (int i = 0; i < (I+1); i++)
                x[j][i] = ++count;
    
    double* ET = malloc((I+1) * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));
    
    double dtau = 0.02;
    double min = -5.;
    double max = 5.;
    double gamma = 0.5;
    
    srand((unsigned int)time(NULL));
    
    DMC(I, i0, Nmax, x, ET, ET_avg, dtau, min, max, gamma);
    
    FILE *fp;
    fp = fopen("DMC.csv", "wb");
    for (int i=0; i<I; i++) {
        fprintf(fp, "%f\n", ET_avg[i]);
    }
    fclose(fp);
    return 0;
}
