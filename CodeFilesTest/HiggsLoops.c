//
//  HiggsLoops.c
//  
//
//  Created by Valerie Plaus on 7/5/16.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "HiggsLoops.h"



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double h2glgl(double mhin,double ytin,double ybin){
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 gluons, based on effective top (ytin) and bottom (ybin)
//c yukawa couplings.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
// Original code by Gabe Shaughnessy, UW Madison, 2011,
//    adapted for C by Valerie Plaus, Wittenberg University, 2016
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    


//    double ytin   input   Effective top coupling for the Higgs mass eigenstate
//    double ybin   input   Effective bottom coupling for the Higgs mass eigenstate
//    double mhin   input   Mass of the Higgs mass eigenstate
double PIVALUE = 4.0*ATAN(1.0);
    double tau; //This variable calculates ability for a particle to split into particle/antiparticle pairs.  If >1, the interaction amplitude is real; if <1, the interaction aplitude is complex.
    double gwsm; //Not sure yet
    double complex amp;
    double complex ffermion,fvector;
    double ampsq,ampsqSM;
    double prefactor;

//cccc	SM Reference
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
    double g2sm = sqrt(4.0*PIVALUE*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
    double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
    double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    gwsm = vev*g2sm; //Fermi coupling??
    prefactor = (alphas*alphas)*(g2sm*g2sm)*(mhin*mhin*mhin)*(1.0/(512.0*(PIVALUE*PIVALUE*PIVALUE)*(mw*mw)))*2.0;
    

    amp = 0.0+0.0*I;
//c	top contribution to the interaction matrix element in the Standard Model
    tau = 4.0*(mt*mt)*(1.0/(mhin*mhin));
    amp = amp + 3.0*ffermion(tau);
    //3 for the 3 colours of quarks, 2/3 for the charge of up-type quarks, and the fermionic form factor for the loop.

//c	bottom contribution to the interaction matrix element in the Standard Model
    tau = 4.0*(mb*mb)*(1.0/(mhin*mhin));
    amp = amp + 3.0*ffermion(tau);
    //3 for the 3 colours of quarks, -1/3 for the charge of down-type quarks, and the fermionic form factor for the loop.

    ampsqSM = amp*conj(amp)*prefactor;


    amp = 0.0+0.0*I;
//c	top contribution to the interaction matrix element in the new model
    tau = 4.0*(mt*mt)*(1.0/(mhin*mhin));
    amp = amp + 3.0*(ytin/ytsm)*ffermion(tau);
    //3 for the 3 colours of quarks, 2/3 for the charge of up-type quarks, and the fermionic form factor for the loop.


//c	bottom contribution to the interaction matrix element in the new model
    tau = 4.0*(mb*mb)*1.0/(mhin*mhin);
    amp = amp + 3.0*(ybin*(1.0/ybsm))*ffermion(tau);
    //3 for the 3 colours of quarks, -1/3 for the charge of down-type quarks, and the fermionic form factor for the loop.

    ampsq = amp*conj(amp)*prefactor;

//Ratio of the gluon fusion to higgs with respect to the Standard Model
    h2glgl = ampsq/ampsqSM;

    return h2glgl;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double h2gaga(double mhin,double ytin,double ybin,double gw) {
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 photons, based on effective top and bottom yukawa couplings,
//c and scaled weak coupling.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double ytin   input   Effective top coupling for the Higgs mass eigenstate
//    double ybin   input   Effective bottom coupling for the Higgs mass eigenstate
//    double mhin   input   Mass of the Higgs mass eigenstate
    
double PIVALUE = 4.0*ATAN(1.0);
    double gwsm;  //not sure
    double tau; //This variable calculates ability for a particle to split into particle/antiparticle pairs.  If >1, the particles and the interaction amplitude are real; if <1, the particles are virtual and the interaction amplitude is complex.
    double complex amp;
//    double complex ffermion,fvector; These are functions, not variables
    double ampsq,ampsqSM;
    double prefactor;


//cccc	SM Reference
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
    double g2sm = sqrt(4.0*PIVALUE*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
    double ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
    double ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    gwsm = vev*g2sm; //Fermi coupling??
    prefactor = (alphaw*alphaw)*(g2sm*g2sm)*(mhin*mhin*mhin) *1.0/(1024.0*(PIVALUE*PIVALUE*PIVALUE)*(mw*mw))*2.0;
    // This calculation includes all the couplings for each vertex

    amp = 0.0+0.0*I; //Initialize the amplitude and allow each contribution to be accumulate for each virtual/real particle in the loop.

//c	top contribution to the interaction matrix element in the Standard Model
    tau = 4.0*(mt*mt)/(mhin*mhin);
    amp = amp + 3.0*pow((2.0/3.0),2)*ffermion(tau);
    //3 for the 3 colours of quarks, 2/3 for the charge of up-type quarks, and the fermionic form factor for the loop.

//c	bottom contribution to the interaction matrix element in the Standard Model
    tau = 4.0*(mb*mb)/(mhin*mhin);
    amp = amp + 3.0*pow((-1.0/3.0),2)*ffermion(tau);
    //3 for the 3 colours of quarks, -1/3 for the charge of down-type quarks, and the fermionic form factor for the loop.

//c	W boson contribution to the interaction matrix element in the Standard Model
    tau = 4.0*(mw*mw)/(mhin*mhin);
    amp = amp + (1.0*1.0)*fvector(tau);
    //1 for the colourless-ness of the W, 1 for the charge of the W, and the vector form factor for the loop.

    ampsqSM = amp*conj(amp)*prefactor;

    amp = 0.0+0.0*I;
//c	top contribution to the interaction matrix element in the new model
    tau = 4.0*(mt*mt)/(mhin*mhin);
    amp = amp + 3.0*(ytin/ytsm)*pow((2.0/3.0),2)*ffermion(tau);
    //3 for the 3 colours of quarks, 2/3 for the charge of up-type quarks, and the fermionic form factor for the loop.


//c	bottom contribution to the interaction matrix element in the new model
    tau = 4.0*(mb*mb)/(mhin*mhin);
    amp = amp + 3.0*(ybin/ybsm)*pow((-1.0/3.0),2)*ffermion(tau);
    //3 for the 3 colours of quarks, -1/3 for the charge of down-type quarks, and the fermionic form factor for the loop.


//c	W boson contribution to the interaction matrix element in the new model
    tau = 4.0*(mw*mw)/(mhin*mhin);
    amp = amp + (gw/gwsm)*(1.0*1.0)*fvector(tau);
    //1 for the colourless-ness of the W, 1 for the charge of the W, and the vector form factor for the loop.

    ampsq = amp*conj(amp)*prefactor;


    h2gaga = ampsq/ampsqSM;

    return h2gaga;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex fvector(double complex tau) {
// Calculates the form factor for a vector-like particle loop.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double complex func;

    double complex formvector;
    formvector = 2.0+3.0*tau+3.0*tau*(2.0-tau)*func(tau);
    return formvector;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex ffermion(double complex tau){
// Calculates the form factor for a fermionic particle loop.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double complex func;

    double complex formfermion;
    formfermion = -2.0*tau*(1.0+(1.0-tau)*func(tau));
    return formfermion;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex fscalar(double complex tau) {
// Calculates the form factor for a salar particle loop.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double complex func;

    double complex formscalar;
    formscalar = tau*(1.0-tau*func(tau));
//	fscalar = (1.0,0.0)
    return formscalar;

}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex func(double tau) {
// This function calculates what happens when the particles flowing in the loop are virtual (initial state mass<4*intermediate state mass) or real.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double complex what;
    double PIVALUE = 4.0*ATAN(1.0);
  
    if(tau>1.0) {
        what = pow(asin(sqrt(1.0/tau)),2)+ 0.0*I;
    }
    else {
        what = -.250*(log((1.0+sqrt(1.0-tau))/(1.0-sqrt(1.0-tau)))-(PIVALUE*I));
    }

    return what;
}

