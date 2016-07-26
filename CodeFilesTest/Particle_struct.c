//
//  Particle_struct.c
//  
//
//  Created by Valerie Plaus on 7/19/16.
//
//

#include "Particle_struct.h"

typedef struct {
    double  mass; //The mass eigenvalue of this particle state.
    double Y_d; //The effective Yukawa couplings for down with this particle.
    double Y_u; //The effective Yukawa couplings for up with this particle.
    double Y_s; //The effective Yukawa couplings for strange with this particle.
    double Y_c; //The effective Yukawa couplings for charm with this particle.
    double Y_b; //The effective Yukawa couplings for bottom with this particle.
    double Y_t; //The effective Yukawa couplings for top with this particle.
    int evec_size;  //This tells us the number of elements in the eigenvector.
    double* eigenvec; //The vector composition of this particle state.
    int bf_size; //This tells us the number of decay modes of the given particle.
    double* branching_frac;  //The array containing the various decay modes of the particle.
} particle;

//Thing to remember: there will need to be a loop creating the array of particles after we cal the eigenvector/eigenvalue function.


/*

struct { // These would either be input by the user or randomized.  Each pointer points to an array of real numbers, and the number of elements is equal to the number of Higgs fields, if the model is holomorphic, or effectively so.  If not, the number of elements would be twice the number of Higgs fields.
    double* Y_ud; //The Yukawa couplings for the first quark family, up and down, for each Higgs field.
    double* Y_cs; //The Yukawa couplings for the second quark family, the charm and strange quarks, for each Higgs field.
    double* Y_tb;  //The  Yukawa couplings for the third quark family, the top and bottom quarks, for each Higgs field.
} Yukawa;
*/
 
 
void calc_yukawa (Particle *Higgs, double *Y1, double *Y2, double *Y3){
    /*
     evec   input   Pointer to the array for the eigenvector composition of the particle.
     Y      input   Pointer to the array for the yukawas for whichever family, up/down, charm/strange, top/bottom.
     Y1     output  Calculates the effective Yukawa coupling for the down-type quark. Needs the location in the struct for the particle.
     Y2     output  Calculates the effective Yukawa coupling for the up-type quark.  Needs the location in struct for the particle.
     */
    const double n=sizeof(evec)/sizeof(evec[0]);
    int i;
    for (i=0;i<n;i+=2){
        Higgs.Y_d += Higgs.evec[i]*Y1[i];
        Higgs.Y_u += Higgs.evec[i+1]*Y1[i+1];
        Higgs.Y_s += Higgs.evec[i]*Y2[i];
        Higgs.Y_c += Higgs.evec[i+1]*Y2[i+1];
        Higgs.Y_b += Higgs.evec[i]*Y3[i];
        Higgs.Y_t += Higgs.evec[i+1]*Y3[i+1];
    }
    return;
}

void fill_struct (const double m, double *A, double *Y1, double *Y2, double *Y3, Particle *Higgs){
    /*
     m      input   The eigenvalue of the mass state.
     A      input   Pointer to the eigenvector of the mass state.
     Y1     input   Pointer to the array holding the first family Yukawa couplings.
     Y2     input   Pointer to the array holding the second family Yukawa couplings.
     Y3     input   Pointer to the array holding the third family Yukawa couplings.
     Higgs  output  Pointer to the struct for the particle in question.
     */
    Higgs.mass=m;
    Higgs.evec=A;
    Higgs.evec_size=sizeof(A)/sizeof(A[0]);
    calc_yukawa(Higgs,Y1,Y2,Y3);
    //Branching fraction function would need to be called here. (BFF)
    //Higgs.bf_size=sizeof(branching_frac)/sizeof(branching_frac[0]);
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
