//
//  Particle_struct.c
//  
//
//  Created by Valerie Plaus on 7/19/16.
//
//

#include "Particle_struct.h"


//Thing to remember: there will need to be a loop creating the array of particles after we call the eigenvector/eigenvalue function.


/*

struct { // These would either be input by the user or randomized.  Each pointer points to an array of real numbers, and the number of elements is equal to the number of Higgs fields, if the model is holomorphic, or effectively so.  If not, the number of elements would be twice the number of Higgs fields.
    double* Y_ud; //The Yukawa couplings for the first quark family, up and down, for each Higgs field.
    double* Y_cs; //The Yukawa couplings for the second quark family, the charm and strange quarks, for each Higgs field.
    double* Y_tb;  //The  Yukawa couplings for the third quark family, the top and bottom quarks, for each Higgs field.
} Yukawa;
*/
 
 
void calc_yukawa (particle *Higgs, double *Y1, double *Y2, double *Y3, double *YL){
    /*
     Higgs   input/output   Pointer to the struct for the particle.
     Y1      input   Pointer to the array for the yukawas for 1st family, down/up.
     Y2      input   Pointer to the array for the yukawas for 2nd family, strange/charm.
     Y3      input   Pointer to the array for the yukawas for 3rd family, bottom/top.
     */
    const double n=sizeof(Higgs[0].evec)/sizeof(Higgs[0].evec[0]);
    int i;
    for (i=0;i<n;i+=2){
        Higgs[0].Y_d += Higgs[0].evec[i]*Y1[i];
        //Contracts the down-type Higgs field(s) with corresponding Yukawa couplings for down quark.
        Higgs[0].Y_u += Higgs[0].evec[i+1]*Y1[i+1];
        //Contracts the up-type Higgs field(s) with corresponding Yukawa couplings for up quark.
        Higgs[0].Y_s += Higgs[0].evec[i]*Y2[i];
        Higgs[0].Y_c += Higgs[0].evec[i+1]*Y2[i+1];
        Higgs[0].Y_b += Higgs[0].evec[i]*Y3[i];
        Higgs[0].Y_t += Higgs[0].evec[i+1]*Y3[i+1];
        Higgs[0].Y_tau += Higgs[0].evec[i]*YL[i];
        Higgs[0].Y_mu += Higgs[0].evec[i+1]*YL[i+1];
    }
    return;
}

void fill_struct (const double m, double *A, double *Y1, double *Y2, double *Y3, double *YL, particle *Higgs){
    /*
     m      input   The eigenvalue of the mass state.
     A      input   Pointer to the eigenvector of the mass state.
     Y1     input   Pointer to the array holding the first family Yukawa couplings.
     Y2     input   Pointer to the array holding the second family Yukawa couplings.
     Y3     input   Pointer to the array holding the third family Yukawa couplings.
     Higgs  output  Pointer to the struct for the particle in question.
     */
    Higgs[0].mass=m;
    Higgs[0].evec=A;
    Higgs[0].evec_size=4;//sizeof(A)/sizeof(A[0]);
    calc_yukawa(Higgs,Y1,Y2,Y3,YL);
// Moving these parts to the main, eventually.
    Higgs[0].decay=(double*)malloc(sizeof(double)*Higgs[0].decay_size);
    Higgs[0].branching_frac=(double*)malloc(sizeof(double)*Higgs[0].bf_size);
    memset(Higgs[0].decay, 0, sizeof(double)*Higgs[0].decay_size);
    memset(Higgs[0].branching_frac, 0, sizeof(double)*Higgs[0].bf_size);
    //Branching fraction function could to be called here. (BFF)
}

double dot_prod (double *A, double *B){
    double Ans;
    int i;
    int n=sizeof(A)/sizeof(A[0]);
    int m=sizeof(B)/sizeof(B[0]);
    if(m!=n){
        printf ("The arrays are not the same length.");
        return 0;
    }
    else {
        for (i=0;i<=n;i++){
            Ans+=A[i]*B[i];
        }
    }
    return Ans;
}

double W_H_rotation (double *A, double *B){
    double Ans;
    int i;
    int n=sizeof(A)/sizeof(A[0]);
    int m=sizeof(B)/sizeof(B[0]);
    if(m!=n){
        printf ("The arrays are not the same length.");
        return 0;
    }
    else {
        for (i=0;i<=n;i+=2){
            Ans+=A[i]*B[i];
            Ans-=A[i+1]*B[i+1];
        }
    }
    return Ans;
}
