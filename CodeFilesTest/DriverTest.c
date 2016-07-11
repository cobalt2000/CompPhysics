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
    
    double vev = 246.0;
    double mt = 173; //Mass of the top particle, Particle Data Group 2015 Summary Tables
    double mb = 4.66; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
    double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
    double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    double Yteff=ytsm*2.0; // Effective top coupling in the new model
    double Ybeff=ybsm*2.0; // Effective bottom coupling in the new model
    double gw=1.0;
    double higgsmass=130; // Test mass top coupling in the new model
    double xsgaga;
    double xsglgl;
    
    printf("The SM and effective top couplings are: %lg and %lg", ytsm, Yteff);
    printf("The SM and effective bottom couplings are: %lg and %lg", ybsm, Ybeff);
    printf("The sample Higgs mass is: %lg", higgsmass);
    
    
    // Call solver
    xsgaga=h2gaga( higgsmass, Yteff, Ybeff, gw );
    xsglgl=h2glgl( higgsmass, Yteff, Ybeff );
    
    // Print solution
    printf("The ratio of the Higgs production from gamma-gamma with respect to the SM is %lg",xsgaga);
    printf("The ratio of the Higgs production from gluon-gluon with respect to the SM is %lg",xsglgl);
    
    
    return 0;
    
}
