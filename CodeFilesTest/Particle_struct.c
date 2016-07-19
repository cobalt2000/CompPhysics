//
//  Particle_struct.c
//  
//
//  Created by Valerie Plaus on 7/19/16.
//
//

#include "Particle_struct.h"

datatype struct {
    double  mass; //The mass eigenvalue of this particle state.
    double* eigenvec; //The vector composition of this particle state.
    double Y_ud; //The effective Yukawa couplings for up with this particle.
    double Y_d; //The effective Yukawa couplings for down with this particle.
    double Y_c; //The effective Yukawa couplings for charm with this particle.
    double Y_s; //The effective Yukawa couplings for strange with this particle.
    double Y_t; //The effective Yukawa couplings for top with this particle.
    double Y_b; //The effective Yukawa couplings for bottom with this particle.
} particle;

datatype struct { // These would either be input by the user or randomized.  Each pointer points to an array of real numbers, and the number of elements is equal to the number of Higgs fields, if the model is holomorphic, or effectively so.  If not, the number of elements would be twice the number of Higgs fields.
    double* Y_ud; //The Yukawa couplings for the first quark family, up and down, for each Higgs field.
    double* Y_cs; //The Yukawa couplings for the second quark family, the charm and strange quarks, for each Higgs field.
    double* Y_tb;  //The  Yukawa couplings for the third quark family, the top and bottom quarks, for each Higgs field.
} Yukawa;

void Calc_Yukawa (double *evec, double *Y, double *Y1, double *Y2){
    double n=sizeof(evec)/sizeof(evec[0]);
    int i;
    for (i=0;i<n;i+=2){
        Y1 += evec[i]*Y[i];
        Y2 += evec[i+1]*Y[i+1];
    }
    return;
}
