//
//  HiggsLoops.h
//  
//
//  Created by Valerie Plaus on 7/5/16.
//
// I'm declaring and defining a lot of values here because then I can just copy and paste them into different header files.

#ifndef HiggsLoops_h
#define HiggsLoops_h

#include <stdio.h>
#include <math.h>

#endif /* HiggsLoops_h */

double pi = 4.0*atan(1.0);
double vev = 246.0;
double mw = 80.38; //Mass of the W particle, Particle Data Group 2015 Summary Tables
double mz = 91.187; //Mass of the Z particle, Particle Data Group 2015 Summary Tables
double mt = 173; //Mass of the top particle, Particle Data Group 2015 Summary Tables
double mb = 4.66; //Mass of the bottom particle, Particle Data Group 2015 Summary Tables
double mc = 1.27; //Mass of the charm particle, Particle Data Group 2015 Summary Tables
double ms = 0.095; //Mass of the strange particle, Particle Data Group 2015 Summary Tables
double xw = 0.2312; //weak mixing angle from Particle Data Group 2015 physical constants
double alphaw = 1.0/128; //fine structure constant Particle Data Group 2015 physical constants
double alphas = 0.118; //strong coupling constant Particle Data Group 2015 physical constants
double g2sm = sqrt(4.0*pi*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling

double h2glgl(double mhin,double ytin,double ybin)
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 gluons, based on effective top (ytin) and bottom (ybin)
//c yukawa couplings.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
// Original code by Gabe Shaughnessy, UW Madison, 2011,
//    adapted for C by Valerie Plaus, Wittenberg University, 2016

double h2gaga(double mhin,double ytin,double ybin,double gw)
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 photons, based on effective top and bottom yukawa couplings,
//c and scaled weak coupling.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
// Original code by Gabe Shaughnessy, UW Madison, 2011,
//    adapted for C by Valerie Plaus, Wittenberg University, 2016

double complex fvector(double tau)

double complex ffermion(double tau)

double complex fscalar(double tau)

double complex func(double complex tau) 
