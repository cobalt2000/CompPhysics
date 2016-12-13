//
//  DriverTest.c
//  
//
//  Created by Valerie Plaus on 7/11/16.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "b_to_sgamma.h"



/*****************************************
 A first try at writing a driver.
 *****************************************/

void main(){
    
//    int i,j, k;  // Standard indexes
//    int info;    // This is a flag used to tell if the function worked correctly.
    
    const double PIVALUE = 4.0*atan(1.0);
    const double vev = 246.0;
    const double xw = 0.2312; //weak mixing angle from Particle Data Group 2015 physical constants
    const double alphaw = 1.0/128; //fine structure constant Particle Data Group 2015 physical constants
    const double g2sm = sqrt(4.0*PIVALUE*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
    const double gwsm = vev*g2sm; //Fermi coupling??
    const double mt = 173; //Mass of the top particle, Particle Data Group 2015 Summary Tables
    const double mb = 4.66; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
    const double mc = 1.41; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
    const double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
    const double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    const double ycsm = sqrt(2.0)*mc/vev; //Standard Model value for the bottom Yukawa coupling
    const double Yteff=ytsm*2.0; // Effective top coupling in the new model
    const double Ybeff=ybsm*2.0; // Effective bottom coupling in the new model
    const double Yceff=ycsm*2.0; // Effective bottom coupling in the new model
    const double gweff=2.0*gwsm;
    const double higgsmass=130; // Test mass top coupling in the new model
    double BF;
    double pull;
    int i;
    
    particle *chhiggs;
    chhiggs = (particle*)malloc( sizeof(particle));
    double *vector;
    vector = (double*)malloc(sizeof(double)*4);
        
    for(i=0; i<4; i++){
        vector[i]=1.0/2.0;
    }
        // initialize each element in the array to an even-weighted, normalized eigenvector
        
    chhiggs[0].mass=higgsmass;
    chhiggs[0].evec=vector;
    chhiggs[0].Y_d=0;
    chhiggs[0].Y_u=0;
    chhiggs[0].Y_s=0;
    chhiggs[0].Y_c=Yceff;
    chhiggs[0].Y_b=Ybeff;
    chhiggs[0].Y_t=Yteff;
    
        
    printf("The charged higgs mass is: %lg . \n", chhiggs[0].mass);
    printf("The effective yukawa couplings are: %lg and %lg . \n", ybsm, Ybeff);
    
    // Call solver
    double *ptr_BF=&BF;
    double *ptr_pull=&pull;
    bsg_nlo(chhiggs, ptr_BF, ptr_pull);
    
    // Print solution
    //printf("The brancing fractions are\n")
    //printf("The ratio of the Higgs production from gamma-gamma with respect to the SM is %lg .\n",xsgaga);
    //printf("The ratio of the Higgs production from gluon-gluon with respect to the SM is %lg .\n",xsglgl);
    
    
    return;
    
}
