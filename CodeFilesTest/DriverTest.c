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
#include "HiggsLoops.h"



/*****************************************
 A first try at writing a driver.
 *****************************************/

int main(){
    
//    int i,j, k;  // Standard indexs
//    int info;    // This is a flag used to tell if the function worked correctly.
    
    const double PIVALUE = 4.0*atan(1.0);
    const double vev = 246.0;
    const double mt = 173; //Mass of the top particle, Particle Data Group 2015 Summary Tables
    const double mb = 4.66; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
    const double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
    const double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    const double xw = 0.2312; //weak mixing angle from Particle Data Group 2015 physical constants
    const double alphaw = 1.0/128; //fine structure constant Particle Data Group 2015 physical constants
    const double g2sm = sqrt(4.0*PIVALUE*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
    const double gwsm = vev*g2sm; //Fermi coupling??
    const double Yteff=ytsm*2.0; // Effective top coupling in the new model
    const double Ybeff=ybsm*2.0; // Effective bottom coupling in the new model
    const double gweff=2.0*gwsm;
    const double higgsmass=130; // Test mass top coupling in the new model
    double xsgaga;
    double xsglgl;
    
    printf("The SM and effective top couplings are: %lg and %lg . \n", ytsm, Yteff);
    printf("The SM and effective bottom couplings are: %lg and %lg . \n", ybsm, Ybeff);
    printf("The sample Higgs mass is: %lg . \n", higgsmass);
    
    
    // Call solver
    xsgaga=h2gaga( higgsmass, Yteff, Ybeff, gweff );
    xsglgl=h2glgl( higgsmass, Yteff, Ybeff );
    
    // Print solution
    printf("The ratio of the Higgs production from gamma-gamma with respect to the SM is %lg .\n",xsgaga);
    printf("The ratio of the Higgs production from gluon-gluon with respect to the SM is %lg .\n",xsglgl);
    
    
    return 0;
    
}
