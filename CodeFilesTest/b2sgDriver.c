//
//  b2sgDriver.c
//  
//
//  Created by Valerie Plaus on 7/11/16.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "Particle_struct.h"
#include "b_to_sgamma.h"



/*****************************************
 A first try at writing a driver for b->s gamma.
 *****************************************/

int main(){
    
//    int i,j, k;  // Standard indexs
//    int info;    // This is a flag used to tell if the function worked correctly.
    
    const double PIVALUE = 4.0*atan(1.0);
    const double vev = 246.0;
    const double tanb = 2;
    const double v1 = vev/sqrt(1.0+tanb*tanb);
    const double v2 = (tanb)*vev/sqrt(1.0+tanb*tanb);
    const double mt = 173; //Mass of the top particle, Particle Data Group 2015 Summary Tables
    const double mb = 4.66; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
    const double mc = 1.41; //Mass of the charm particle, Particle Data Group 2015 Summary Tables
    const double ms = 0.096; //Mass of the strange particle, Particle Data Group 2016 Summary Tables
    const double mu = 0.0022; //Mass of the up particle, Particle Data Group 2016 Summary Tables
    const double md = 0.0047; //Mass of the down particle, Particle Data Group 2016 Summary Tables
    const double mtau = 1.78; //Mass of the up particle, Particle Data Group 2016 Summary Tables
    const double mmu = 0.105; //Mass of the down particle, Particle Data Group 2016 Summary Tables
//    const double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
//    const double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
//    const double xw = 0.2312; //weak mixing angle from Particle Data Group 2015 physical constants
//    printf("No chhiggs yet.\n");
    particle *chhiggs;
//    printf("No malloc yet.\n");

    chhiggs=(particle*)malloc(sizeof(particle)*1);
//    printf("Malloc works?\n");

    chhiggs[0].mass=730; //The mass eigenvalue of this particle state.
//    printf("Mass win.\n");

    chhiggs[0].Y_d=sqrt(2.0)*md/v1; //The effective Yukawa couplings for down with this particle.
    chhiggs[0].Y_u=sqrt(2.0)*mu/v2; //The effective Yukawa couplings for up with this particle.
    chhiggs[0].Y_s=sqrt(2.0)*ms/v1; //The effective Yukawa couplings for strange with this particle.
    chhiggs[0].Y_c=sqrt(2.0)*mc/v2; //The effective Yukawa couplings for charm with this particle.
    chhiggs[0].Y_b=sqrt(2.0)*mb/v1; //The effective Yukawa couplings for bottom with this particle.
    chhiggs[0].Y_t=sqrt(2.0)*mt/v2; //The effective Yukawa couplings for top with this particle.
    chhiggs[0].Y_tau=sqrt(2.0)*mtau/v1; //The effective Yukawa couplings for tau with this particle.
    chhiggs[0].Y_mu=sqrt(2.0)*mmu/v1; //The effective Yukawa couplings for muon with this particle.
    chhiggs[0].evec_size=1;  //This tells us the number of elements in the eigenvector.
 //   printf("No eigenvector yet.\n");

    chhiggs[0].evec =(double*)malloc(sizeof(double)*10); //The vector composition of this particle state.
    
//    int bf_size; //This tells us the number of decay modes of the given particle.
//    int decay_size; //This tells us the number of decay modes of the given particle.
//    double* decay;  //The array containing the various decay modes of the particle.
//    double* branching_frac;  //The array containing the various decay modes of the particle.
//    printf("Malloc is happy after initialization.\n");
    double BF;
    double pull;
    
//    printf("The SM and effective top couplings are: %lg and %lg . \n", ytsm, Yteff);

    
    
    // Call solver
    bsg_nlo( chhiggs, &BF, &pull );
    printf("Malloc is happy after call.\n");

    // Print solution
    printf("The pull for the loop order correction for the charge higgs loop is %lg .\n",pull);
    
    
    return 0;
    
}
