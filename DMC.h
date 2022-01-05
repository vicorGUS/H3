//
//  DMC.h
//  H3
//
//  Created by Victor Gustafsson on 2021-12-21.
//

#ifndef DMC_h
#define DMC_h

#include <stdio.h>

void DMC(int I, int i0, int Nmax, double* x[Nmax], double* ET, double* E_avg, double dtau, double min, double max, double gamma);
double randn (double mu, double sigma);

#endif /* DMC_h */
