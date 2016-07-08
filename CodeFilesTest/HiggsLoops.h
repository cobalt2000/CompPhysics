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

const double PIVALUE 4.0*ATAN(1.0);

double h2glgl(double mhin,double ytin,double ybin);
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 gluons, based on effective top (ytin) and bottom (ybin)
//c yukawa couplings.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
// Original code by Gabe Shaughnessy, UW Madison, 2011,
//    adapted for C by Valerie Plaus, Wittenberg University, 2016

double h2gaga(double mhin,double ytin,double ybin,double gw);
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 photons, based on effective top and bottom yukawa couplings,
//c and scaled weak coupling.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
// Original code by Gabe Shaughnessy, UW Madison, 2011,
//    adapted for C by Valerie Plaus, Wittenberg University, 2016

double complex fvector(double tau);

double complex ffermion(double tau);

double complex fscalar(double tau);

double complex func(double complex tau);

#endif /* HiggsLoops_h */

