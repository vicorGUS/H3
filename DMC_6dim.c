#include "DMC.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int pos_idx(int i, int j, int k, int I, int N){
    return (i * 6 * N) + (j * 6) + k;
}

double r(double x, double y, double z){
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

double V(double R[6]){
    double r1 = r(R[0], R[1], R[2]);
    double r2 = r(R[3], R[4], R[5]);
    double r12 = r(R[3] - R[0], R[4] - R[1], R[5] - R[2]);

    return -2./r1 - 2./ r2 + 1./r12;
}

double * diff(double R[6], double dtau){
    static double R_R[6];
    for (int k = 0; k<6; k++){
        R_R[k] = sqrt(dtau) * randn(0., 1.) + R[k];
    }
    return R_R;
}

double * drift(double R[6], double R2[6], double a, double dtau){
    static double R_F[6]; 
    double r1 = r(R[0], R[1], R[2]);
    double r2 = r(R[3], R[4], R[5]);
    double r12 = r(R[3] - R[0], R[4] - R[1], R[5] - R[2]);

    for (int k = 0; k<3; k++){
        R_F[k] = (-2. * R[k] / r1 - 1./(2. * pow(1.+a*r12,2)) * (R[k+3] - R[k])/r12) *dtau + R2[k];
        R_F[k+3] = (-2. * R[k+3] / r2 - 1./(2. * pow(1.+a*r12,2)) * (R[k+3] - R[k])/r12) * dtau + R2[k+3];
    }

    return R_F;
}

double EL(double R[6], double a){
    double r1 = r(R[0], R[1], R[2]);
    double r2 = r(R[3], R[4], R[5]);
    double r12 = r(R[3] - R[0], R[4] - R[1], R[5] - R[2]);

    double EL = -4. - 1./ (r12 * pow(1.+a*r12,3)) - 1./ (4. * pow(1.+a*r12,4)) + 1./r12;
    for(int k=0; k<3; k++){
        EL += (R[k]/r1 - R[k+3]/r2) * (R[k] - R[k+3]) / (r12*pow(1.+a*r12,2));
    }
    return EL;
}

void DMC_6dim(int I, int i0, int Nmax, double* pos, double dtau, double gamma, double a, double ET_0, int N_0){

    int m;
    double R1 = 0.7 + (double) rand() / RAND_MAX;
    double R2 = 0.7 + (double) rand() / RAND_MAX;
    double t1 = acos(2. * (double) rand() / RAND_MAX - 1);
    double t2 = acos(2. * (double) rand() / RAND_MAX - 1);
    double p1 = 2. * 3.149265 * (double) rand() / RAND_MAX;
    double p2 = 2. * 3.149265 * (double) rand() / RAND_MAX;
    double* R = malloc(6* sizeof(double));
    //double* R_2 = malloc(6* sizeof(double));

    double* ET = malloc((I+1) * sizeof(double));
    double* ET_avg = malloc((I+1) * sizeof(double));
    int* N = malloc((I+1) * sizeof(int));
    
    for (int i=1; i<(I+1); i++) {
        N[i] = 0;
    }

    N[0] = N_0;
    ET[0] = ET_0;

    for (int j=0; j<N[0]; j++){
        pos[pos_idx(0,j,0,I,Nmax)] = (double) R1 * sin(t1) * cos(p1); //x1
        pos[pos_idx(0,j,1,I,Nmax)] = (double) R1 * sin(t1) * sin(p1); //y1
        pos[pos_idx(0,j,2,I,Nmax)] = (double) R1 * cos(t1);           //z1
        pos[pos_idx(0,j,3,I,Nmax)] = (double) R2 * sin(t2) * cos(p2); //x2
        pos[pos_idx(0,j,4,I,Nmax)] = (double) R2 * sin(t2) * sin(p2); //y2
        pos[pos_idx(0,j,5,I,Nmax)] = (double) R2 * cos(t2);           //z2
    }
    
    for (int j=N[0]; j<Nmax; j++) { 
        for (int k=0; k<6; k++ ){
            pos[pos_idx(0,j,k,I,Nmax)] = 0.0;
        }
    }
    
    printf("START: N: %d, E: %f\n", N[0], ET[0]);
    
    for (int i=0; i<I; i++){    
        ET_avg[i] = 0;
        for (int j=0; j<N[i]; j++){

            for(int k = 0; k<6; k++){
                R[k] = pos[pos_idx(i,j,k,I,Nmax)];
            }

            // Using the spherical potential
            /*
            R = diff(R, dtau);
            double W = exp( -(V(R) - ET[i]) * dtau);
            */

            // With importance sampling
            /* 
            R = drift(R, R, a, dtau);
            R = diff(R, dtau);
            double W = exp( -(EL(R, a) - ET[i]) * dtau);
            */

            // With importance sampling more exact
            //* 
            R = drift(R, R, a, dtau);
            R = diff(R, dtau);
            double W = exp( -(EL(R, a) - ET[i]) * dtau);
            //*/

            // Add or delete walkers
            m = (int)(W + ((double) rand() / (RAND_MAX)));

            if (abs(m)>1000){
                m = 1;
            }
            else{
                N[i+1] += m;
            }

            // importance sampling more accurate
            R = diff(R, dtau/2.);
            R = drift(R, R, a, dtau/2.);

            // Update the position values
            for(int k = 0; k<6; k++){
                pos[pos_idx(i,j,k,I,Nmax)] = R[k];
            }

            // Duplicate the walkers
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

        for (int i´=0; i´<=i; i´++) {
            ET_avg[i] += (1. / i) * ET[i´];
        }
        
        ET_avg[0] = -3.0;
        ET[i+1] = ET_avg[i] - gamma * log((double)N[i] / N[0]);
        
        if ((i-i0+1) != 0) {
            ET_avg[i+1] = 1 / (i - i0 + 1.) * ET[i+1] + (i - i0) / (i - i0 + 1.) * ET_avg[i];
        }
        
        //printf("i: %d N: %d E: %f\n",i, N[i], ET_avg[i]);
    }
    printf("END:   N: %d E: %f\n",N[I], ET_avg[I]);

    FILE *fp;
    fp = fopen("DMC.csv", "wb");
    for (int i=0; i<I; i++) {
        fprintf(fp, "%i,%f,%f\n",N[i], ET[i], ET_avg[i]);
    }
    fclose(fp);
}