//
//  HiggsLoops.c
//  
//
//  Created by Valerie Plaus on 7/5/16.
//
//



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
    
#include "HiggsLoops.h"


    double ytin,ybin,mhin; //input from the user for the effective top coupling, bottom coupling, and mass of the higgs
    double tau,gwsm;
    double complex amp;
    double complex ffermion,fvector;
    double ampsq,ampsqSM;
    double prefactor;

//cccc	SM Reference: some declared and defined in the header file
//    xw = 0.2312; //weak mixing angle from Particle Data Group 2015 physical constants
//    alphaw = 1.0/128; //fine structure constant Particle Data Group 2015 physical constants
//    alphas = 0.118; //strong coupling constant Particle Data Group 2015 physical constants
//    g2sm = sqrt(4.0*pi*alphaw)/sqrt(xw); //Standard Model value for isospin coupling
//    ytsm = sqrt(2.0)*mt/vev;  //Standard Model value for the top Yukawa coupling
//    ybsm = sqrt(2.0)*mb/vev; //Standard Model value for the bottom Yukawa coupling
    gwsm = vev*g2sm; //Fermi coupling??
    prefactor = (alphas*alphas)*(g2sm*g2sm)*(mhin*mhin*mhin)*(1.0/(512.0*(pi*pi*pi)*(mw*mw)))*2.0;
    

    amp = 0.0+0.0*I;
//c	top contribution in the Standard Model
    tau = 4.0*(mt*mt)*(1.0/(mhin*mhin));
    amp = amp + 3.0*ffermion(tau);

//c	bottom contribution in the Standard Model
    tau = 4.0*(mb*mb)*(1.0/(mhin*mhin));
    amp = amp + 3.0*ffermion(tau);

    ampsqSM = amp*conj(amp)*prefactor;


    amp = 0.0+0.0*I;
//c	top contribution in the new model
    tau = 4.0*(mt*mt)*(1.0/(mhin*mhin));
    amp = amp + 3.0*(ytin/ytsm)*ffermion(tau);


//c	bottom contribution in the new model
    tau = 4.0*(mb*mb)*1.0/(mhin*mhin);
    amp = amp + 3.0*(ybin*(1.0/ybsm))*ffermion(tau);

    ampsq = amp*conj(amp)*prefactor;

//Ratio of the gluon fusion to higgs with respect to the Standard Model
    h2glgl = ampsq/ampsqSM;

return
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double h2gaga(double mhin,double ytin,double ybin,double gw) {
//c These calculate the production ratio of CP even Higgs bosons
//c from 2 photons, based on effective top and bottom yukawa couplings,
//c and scaled weak coupling.
//c Coupling is based on the composition of the CP even Higgs.
//c (This is not the ratio of the couplings)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include 'HiggsLoops.h'

    double ytin,ybin,gwsm,mhin,gw;
    double tau;
    double complex amp;
    double complex ffermion,fvector;
    double ampsq,ampsqSM;
    double alphaw,prefactor;


//cccc	SM Reference
xw = 0.2312;
alphaw = 1.0/128.0;
g2sm = sqrt(4.0*pi*alphaw)*(1.0/sqrt(xw));
ytsm = sqrt(2.0)*mt/vev;
ybsm = sqrt(2.0)*mb/vev;
gwsm = vev*g2sm;
prefactor = (alphaw*alphaw)*(g2sm*g2sm)*(mhin*mhin*mhin) *1.0/(1024.0*(pi*pi*pi)*(mw*mw))*2.0;

    amp = 0.0+0.0*I;
//c	top contribution in the Standard Model
    tau = 4.0*(mt*mt)/(mhin*mhin);
    amp = amp + 3.0*pow((2.0/3.0),2)*ffermion(tau);

//c	bottom contribution in the Standard Model
    tau = 4.0*(mb*mb)/(mhin*mhin);
    amp = amp + 3.0*pow((-1.0/3.0),2)*ffermion(tau);

//c	W boson contribution in the Standard Model
    tau = 4.0*(mw*mw)/(mhin*mhin);
    amp = amp + (1.0*1.0)*fvector(tau);

    ampsqSM = amp*conj(amp)*prefactor;

    amp = 0.0+0.0*I;
//c	top contribution in the new model
    tau = 4.0*(mt*mt)/(mhin*mhin);
    amp = amp + 3.0*(ytin/ytsm)*pow((2.0/3.0),2)*ffermion(tau);


//c	bottom contribution in the new model
    tau = 4.0*(mb*mb)/(mhin*mhin);
    amp = amp + 3.0*(ybin/ybsm)*pow((-1.0/3.0),2)*ffermion(tau);

//c	W boson contribution in the new model
    tau = 4.0*(mw*mw)/(mhin*mhin);
    amp = amp + (gw/gwsm)*(1.0*1.0)*fvector(tau);

    ampsq = amp*conj(amp)*prefactor;


    h2gaga = ampsq/ampsqSM;

return
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex fvector(double tau) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include 'HiggsLoops.h'
    double tau;
    double complex func;

    fvector = 2.0+3.0*tau+3.0*tau*(2.0-tau)*func(tau);
return fvector
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex ffermion(double tau){
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include 'HiggsLoops.h'
    double tau;
    double complex func;

    ffermion = -2.0*tau*(1.0+(1.0-tau)*func(tau));
return ffermion
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex fscalar(double tau) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include 'HiggsLoops.h'
    double tau;
    double complex func;

    fscalar = tau*(1.0-tau*func(tau));
//	fscalar = (1.0,0.0)
return fscalar

}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex func(double complex tau) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include <math.h>
#include 'HiggsLoops.h'
    double tau;

    if(tau>1.0) {
        func = pow(asin(sqrt(1.0/tau)),2)+ 0.0*I;}
    else {
        func = -.250*(log((1.0+sqrt(1.0-tau))/(1.0-sqrt(1.0-tau)))-(pi*I));}
endif

return func
}

