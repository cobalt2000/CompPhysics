//
//  b_to_sgamma.c
//  
//
//  Created by Valerie Plaus on 7/15/16.
//
//

#include "b_to_sgamma.h"
#include <math.h>
#include "Particle_struct.h"
#define DEL .99

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
void bsg_nlo(   particle *chhiggs,double *BF,double *pull){
//    double bsg_nlo(double *xAu,double *xAd,double *xMH,double xnhiggs,double *BF,double *pull){

    /*c
c	Calculates the B -> X_S + gamma BF from the Standard model and 2HDM
c
c	M. Ciuchini, G. Degrassi, P. Gambino, and G. F. Giudice, Nucl. Phys. B527, 21 (1998).
c
c
c	Higgs - Input  - array of structs of type "particle", containing the mass, eigenvector, and effective yukawa couplings for that particle.
c
c	BF - Output - the result BF(B -> X_S + gamma + X)
c	pull - Output - the number of std. dev. away from experiment
c		(taking into account experimental uncertainty and 10% theoretical uncertainty)
c
c
c
c	G. Shaughnessy 5/3/2012
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

    
    
//real*8 xAu(maxhiggs),xAd(maxhiggs),xMh(maxhiggs)
    int i;
//Variables for the final step of the calculation.  They are physically observed constants, but are defined below for thought-flow.
    double BFthref,BFthrefunc, BFexp,BFexpunc,BFcenu,relunc;
//Variables for the intermediate steps of the calculations.
    double delNPc,delNPSL,delNPgam,del;
//Variables holding the final branching braction, the whole point of these functions.
//    double BF, pull;
    
    //Physical constants like pi, or the masses of particles.
    const double pibsg = 4*atan(1);
    
    const double mcbsg = 1.41;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mbartpole = 175;
    const double nf = 5;  //The top does not participate.
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;
    const double alphae = 1.0/130.3;
    
    //c	For comparison w/ the Literature
    //c	mb = 4.8d0
    //c	mc = mb - 3.39d0
    //c	MW = 80.33d0
    
    //The energy scale at which calculations will occur
    const double mubarb = mbbsg;
    const double mub = mbbsg;
    //c	mub = 5d0
    
    //Parameters for the calculation from the paper.  I moved them to the correct sub-programs.
/*    const double nf = 5;  //The top does not participate.
    const double beta0 = 11-2.0/3.0*nf;
    const double beta1 = 102-38.0/3.0*nf;
    
    const double gam0m = 8;
    const double gam1m = 404.0/3.0-40.0/9.0*nf;*/
    
    const double lam1 = 0;	//! drops out anyway
    const double lam2 = 0.12; //! GeV**2

    const double eta = alphasbsg(mwbsg)/alphasbsg(mub); //the ratio of the strong couplings at the two renormalization scales, the W and b poles?
    const double z = (mcbsg*mcbsg)/(mbbsg*mbbsg);

    const double CKMproductSq = 0.95;

    //order = 'NLO';
    /*    if(order='LO') {
     beta1 = 0;
     gam1m = 0;
     }
     */
    
//Now we begin the particle-specific pieces.  Except that now that's relegated to the places were the sums over particles are needed.
/*    nhiggs = sizeof(chhiggs)/sizeof(chhiggs[0]);
    for (i=1;i<nhiggs;i++){
        Au[i] = Higgs.Y_u[i];
        Ad[i] = xAd[i];
        MHbsg[i] = xMH[i];
        //c        write(*,*) xMH(i),xAu(i),xAd(i)
        //c        write(*,*) MHbsg(i),Au(i),Ad(i)
    }*/
    

    delNPc = -1/9*lam2/C0b(7,eta)*(C0b(2,eta) - C0b(1,eta)/6);

    delNPSL = lam1/2 + 3*lam2/2*(1-4*pow((1-z),4)/f(z));
    delNPgam = lam1/2 - 9/2*lam2;
    del = DEL;

//c	write(*,*)(C0beff(i),i=1,8)
//c	write(*,*)C1beff(7)
//c	stop


    BFcenu = 0.1049;	//! from PDG

    *BF = BFcenu * CKMproductSq * 6*alphae/(pibsg*f(z)*kap(z))
        *mbrun(mub)*mbrun(mub)/(mbbsg*mbbsg)*(DtermSq(mub, chhiggs) + Aterm(mub, eta, chhiggs))
        *(1- delNPSL / (mbbsg*mbbsg) + delNPgam/(mbbsg*mbbsg) + delNPc / (mcbsg*mcbsg));

//c	Comparison with Experiment
//c	Taken from HFAG world average: 1010.1589
    BFexp = 3.55*0.0001; //d-4 in Fortran
    BFexpunc = 0.256*0.0001; //d-4 in Fortran

//c	Taken as the SM reference
    BFthref = 3.62*0.0001; //d-4 in Fortran
    BFthrefunc = 0.33*0.0001; //d-4 in Fortran

    relunc = BFthrefunc/BFthref;


    *pull = fabs((*BF-BFexp)/sqrt(BFexpunc*BFexpunc + (*BF * relunc)*(*BF * relunc)));




    return;
    }



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double DtermSq(double Q,    particle *chhiggs) {
    /* Now we're actually getting to a calculation.
     Q  input   energy scale of the calculation.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mcbsg = 1.41;
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;
    const double pibsg = 4*atan(1);
    double sumreal,sumimag;
    double r1r,r1c,r2r,r2c,r7,r8r,r8c;
    double logz;
    double zeta3;
    double dtermreal,dtermimag;
    double moose;

    zeta3 = 1.20206;

    //The energy scale at which calculations will occur
    const double mubarb = mbbsg;
    const double mub = mbbsg;
    //c	mub = 5d0
    const double eta = alphasbsg(mwbsg)/alphasbsg(mub); //the ratio of the strong couplings at the two renormalization scales, the W and b poles?
    const double z = (mcbsg*mcbsg)/(mbbsg*mbbsg);


    logz = log(z);


    r2r = 2/243*(-833+144*pibsg*pibsg*pow(z,1.5)
            + (1728-180*pibsg*pibsg-1296*zeta3+(1296-324*pibsg*pibsg)*logz + 108*log(z)*logz+36*logz*logz*logz)*z
                + (648+72*pibsg*pibsg+(432-216*pibsg*pibsg)*logz+36*logz*logz*logz)*z*z
                + (-54-84*pibsg*pibsg+1092*logz-756*logz*logz)*z*z*z);
    
    r2c = 16*pibsg/81*(-5+(45-3*pibsg*pibsg+9*logz+9*logz*logz)*z
                + (-3*pibsg*pibsg+9*logz*logz)*z*z + (28-12*logz)*z*z*z);
    

    r7 = -10/3-8/9*pibsg*pibsg;

    r8r = -4/27*(-33+2*pibsg*pibsg);
    r8c = 24/27*pibsg;

    r1r = -1/6*r2r;
    r1c = -1/6*r2c;




    sumreal = 0;
    sumreal = sumreal + C0beff(1,eta,chhiggs)*(r1r + gam0eff(1,7)*log(mbbsg/mub));
    sumreal = sumreal + C0beff(2,eta,chhiggs)*(r2r + gam0eff(2,7)*log(mbbsg/mub));
    sumreal = sumreal + C0beff(7,eta,chhiggs)*(r7 + gam0eff(7,7)*log(mbbsg/mub));
    sumreal = sumreal + C0beff(8,eta,chhiggs)*(r8r + gam0eff(8,7)*log(mbbsg/mub));


    sumimag = 0;
    sumimag = sumimag + C0beff(1,eta,chhiggs)*r1c;
    sumimag = sumimag + C0beff(2,eta,chhiggs)*r2c;
    sumimag = sumimag + C0beff(8,eta,chhiggs)*r8c;


    dtermreal = C0beff(7,eta,chhiggs) + alphasbsg(mub)/(4*pibsg)*(C1beff(7,eta,chhiggs) + sumreal);
    dtermimag = alphasbsg(mub)/(4*pibsg)*(sumimag);

    moose = (dtermreal*dtermreal + dtermimag*dtermimag);

//c	write(*,*)C0beff(7),alphas(mub)/(4d0*pi)*C1beff(7),sumreal,sumimag

    return moose;
}



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double gam0eff(int i,int j) {
    //Why does he have J as an input here?  He doesn't use it.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    const double pibsg = 4*atan(1);
    
    const double mcbsg = 1.41;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mbartpole = 175; //"m bar" at the t pole (renormalization scheme)
    const double nf = 5;  //The top does not participate.
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;
    const double alphae = 1.0/130.3;
    
    //c	For comparison w/ the Literature
    //c	mb = 4.8d0
    //c	mc = mb - 3.39d0
    //c	MW = 80.33d0
    
    //The energy scale at which calculations will occur
    const double mubarb = mbbsg;
    const double mub = mbbsg;
    //c	mub = 5d0

    double moose;

    switch(i){
        case 1: moose = -208*1.0/243.0;
        case 2: moose = 416*1.0/81.0;
        case 3: moose = -176*1.0/81.0;
        case 4: moose = -152*1.0/243.0;
        case 5: moose = -6272*1.0/81.0;
        case 6: moose = 4624*1.0/243.0;
        case 7: moose = 32*1.0/3.0;
        case 8: moose = -32*1.0/9.0;
    }


    moose *= alphasbsg(mub)*1.0/(4.0*pibsg);

    return moose;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Aterm(double Q, double eta,   particle *chhiggs) {
    /* Correction for low energy bremsstrahlang.
     Q  input   energy of calculation
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    //Error: Gabe defines Q as an input, and never uses it. Instead, he uses the bottom mass, defined in mub=mbbsg above.  I used Q, so that if you change the calculation energy above, it actually changes here :).
    
    const double pibsg = 4*atan(1);
    
    const double mcbsg = 1.41;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mbartpole = 175;
    const double nf = 5;  //The top does not participate.
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;
    const double alphae = 1.0/130.3;
    const double del = DEL;
    
    //c	For comparison w/ the Literature
    //c	mb = 4.8d0
    //c	mc = mb - 3.39d0
    //c	MW = 80.33d0
    


    double sum,moose;
    int i,j;

    sum = 0;
    for (i=1;i<=8;i++){
        for (j=i;i<=8;j++){
            sum = sum + C0beff(i,eta, chhiggs)*C0beff(j,eta,chhiggs)*fij(i,j);
        }
//c	 write(*,*)i,C0beff(i)
    }

//c	write(*,*)sum

    moose = (exp(-alphasbsg(Q)*log(del)*(7+2*log(del))/(3*pibsg))-1) * C0beff(7,eta,chhiggs)*C0beff(7,eta,chhiggs)
         + alphasbsg(Q)/pibsg * sum;


    return moose;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double mbart(double Q){
    /* Correction for low energy bremsstrahlang.
     Q  input   energy of calculation
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;
    const double pibsg = 4*atan(1);
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mbartpole = 175;
    //Parameters for the calculation from the paper
    const double nf = 5;  //The top does not participate.
    const double beta0 = 11-2.0/3.0*nf;
    const double beta1 = 102-38.0/3.0*nf;
    const double gam0m = 8;
    const double gam1m = 404.0/3.0-40.0/9.0*nf;


    what = mbartpole * pow((alphasbsg(Q)/alphasbsg(mtbsg)),(gam0m/(2*beta0)))*
        (1 + alphasbsg(mtbsg)/(4*pibsg)*gam0m/(2*beta0)*(gam1m/gam0m-beta1/beta0)*
     		(alphasbsg(Q)/alphasbsg(mtbsg)-1));

//c	write(*,*)mbart,mbartpole,alphas(Q),gam0m,beta0,gam1m,beta1,nf
    return what;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double fij(int i,int j) {
    /* List of the matrix elements for the possible interactions.  Each element listed below corresponds to a particular particle interaction and Wilson coefficient.
     i  input   column index for matrix element
     j  input   row index for matrix element
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double A;
    const double del = DEL;

//    if(del!=0.99) {
//write(*,*)'Delta not 0.99'
//stop
//    }
/*
if(i.eq.1.and.j.eq.1) A =  0.0009;
if(i.eq.1.and.j.eq.2) A = -0.0113;
if(i.eq.1.and.j.eq.3) A = -0.0035;
if(i.eq.1.and.j.eq.4) A =  0.0006;
if(i.eq.1.and.j.eq.5) A = -0.0459;
if(i.eq.1.and.j.eq.6) A = -0.0600;
if(i.eq.1.and.j.eq.7) A = -0.0030;
if(i.eq.1.and.j.eq.8) A =  0.0010;

if(i.eq.2.and.j.eq.2) A =  0.0340;
if(i.eq.2.and.j.eq.3) A =  0.0210;
if(i.eq.2.and.j.eq.4) A = -0.0035;
if(i.eq.2.and.j.eq.5) A =  0.2754;
if(i.eq.2.and.j.eq.6) A =  0.3599;
if(i.eq.2.and.j.eq.7) A =  0.0182;
if(i.eq.2.and.j.eq.8) A = -0.0061;

if(i.eq.3.and.j.eq.3) A =  0.0140;
if(i.eq.3.and.j.eq.4) A = -0.0047;
if(i.eq.3.and.j.eq.5) A =  0.3277;
if(i.eq.3.and.j.eq.6) A =  0.0666;
if(i.eq.3.and.j.eq.7) A =  0.0421;
if(i.eq.3.and.j.eq.8) A = -0.0140;

if(i.eq.4.and.j.eq.4) A =  0.0088;
if(i.eq.4.and.j.eq.5) A = -0.0546;
if(i.eq.4.and.j.eq.6) A =  0.1570;
if(i.eq.4.and.j.eq.7) A = -0.0070;
if(i.eq.4.and.j.eq.8) A =  0.0023;

if(i.eq.5.and.j.eq.5) A =  1.9369;
if(i.eq.5.and.j.eq.6) A =  0.9506;
if(i.eq.5.and.j.eq.7) A =  0.5926;
if(i.eq.5.and.j.eq.8) A = -0.1975;


if(i.eq.6.and.j.eq.6) A =  0.9305;
if(i.eq.6.and.j.eq.7) A =  0.0002;
if(i.eq.6.and.j.eq.8) A = -0.0001;

if(i.eq.7.and.j.eq.7) A =  3.4211;
if(i.eq.7.and.j.eq.8) A =  0.3897;

if(i.eq.8.and.j.eq.8) A =  3.2162;
*/
    switch(i){
        case 1: switch (j){
            case 1: A =  0.0009; break;
            case 2: A = -0.0113; break;
            case 3: A = -0.0035; break;
            case 4: A =  0.0006; break;
            case 5: A = -0.0459; break;
            case 6: A = -0.0600; break;
            case 7: A = -0.0030; break;
            case 8: A =  0.0010; break;
        }  break;
    
        case 2: switch(j) {
            case 2: A =  0.0340; break;
            case 3: A =  0.0210; break;
            case 4: A = -0.0035; break;
            case 5: A =  0.2754; break;
            case 6: A =  0.3599; break;
            case 7: A =  0.0182; break;
            case 8: A = -0.0061; break;
        } break;
    
        case 3: switch(j){
            case 3: A =  0.0140; break;
            case 4: A = -0.0047; break;
            case 5: A =  0.3277; break;
            case 6: A =  0.0666; break;
            case 7: A =  0.0421; break;
            case 8: A = -0.0140; break;
        } break;
    
        case 4: switch(j) {
            case 4: A =  0.0088; break;
            case 5: A = -0.0546; break;
            case 6: A =  0.1570; break;
            case 7: A = -0.0070; break;
            case 8: A =  0.0023; break;
        } break;
    
        case 5: switch(j) {
            case 5: A =  1.9369; break;
            case 6: A =  0.9506; break;
            case 7: A =  0.5926; break;
            case 8: A = -0.1975; break;
        } break;
            
    
        case 6: switch(j) {
            case 6: A =  0.9305; break;
            case 7: A =  0.0002; break;
            case 8: A = -0.0001; break;
        } break;
                
        case 7: switch(j) {
            case 7: A =  3.4211;  break;
            case 8: A =  0.3897; break;
        } break;
    
        case 8: switch(j) {
            case 8: A =  3.2162; break;
        }
    }

    return A;
    }



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double mbrun(const double Q) {
    /* The running value of the bottom mass. Eqn 42
     Q   input   energy scale at which the calling function needs the bottom mass.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    const double pibsg = 4*atan(1);

    double moose = Q*(1-4/3*alphasbsg(Q)/pibsg);

    return moose;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double alphasbsg(const double Q){
    /* The running value of the strong coupling constant. Eqn 42
     Q   input   energy scale at which the calling function needs the strong coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    const double pibsg = 4*atan(1);
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;

    double moose;
    double v,alphasMZ;
    //Parameters for the calculation from the paper
    const double nf = 5;  //The top does not participate.
    const double beta0 = 11-2.0/3.0*nf;
    const double beta1 = 102-38.0/3.0*nf;


    alphasMZ = 0.118;

    v = 1+beta0*alphasMZ/(2*pibsg)*log(Q/mzbsg);
//c	write(*,*)v,beta0,beta1,alphasMZ,pi,Q,mz
//c	stop

    moose = alphasMZ / v *(1-beta1/beta0*alphasMZ/(4*pibsg)*log(v)/v);

    return moose;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double kap(const double z) {
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    const double pibsg = 4*atan(1);
    const double mcbsg = 1.41;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mbartpole = 175;
    const double nf = 5;  //The top does not participate.
    const double mzbsg = 91.18;
    const double mwbsg = 80.39;
    const double alphae = 1.0/130.3;
    
    //c	For comparison w/ the Literature
    //c	mb = 4.8d0
    //c	mc = mb - 3.39d0
    //c	MW = 80.33d0
    
    //The energy scale at which calculations will occur
    const double mubarb = mbbsg;
    const double mub = mbbsg;
    //c	mub = 5d0
    double thing = 1 - 2*alphasbsg(mubarb)/(3*pibsg)*h(z)/f(z);

//    if(order='LO') kap = 1;

    return thing;
}



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double f(const double z) {
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing = 1-8*z+8*z*z*z-z*z*z*z-12*z*z*log(z);
    return thing;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double h(const double z) {
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    const double pibsg = 4*atan(1);
    double term1,term2,term3,term4,term5,term6,term7,sum;


    term1 = -(1-z*z)*(25/4-239/3*z+25/4*z*z);
    term2 = z*log(z)*(20+90*z-4/3*z*z+17/3*z*z*z);
    term3 = z*z*(log(z)*log(z))*(36+z*z);
    term4 = (1-z*z)*(17/3-64/3*z+17/3*z*z)*log(1-z);
    term5 = -4*(1+30*z*z+z*z*z*z)*log(z)*log(1-z);
    term6 = -(1+16*z*z+z*z*z*z)*(6*Li2(z) - pibsg*pibsg);
    term7 = -32*pow(z,1.5)*(1+z)*(pibsg*pibsg-4*Li2(sqrt(z))+4*Li2(-sqrt(z))-2*log(z)*log((1-sqrt(z))/(1+sqrt(z))));



    sum = term1 + term2 + term3 + term4 + term5 + term6 + term7;

    return sum;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Li2(const double z) {			//!either positive or negative arguments allowed
    //Undefined in the paper.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing=0;

    if(z<0) {
        thing = 0.5*Li2p(z*z) - Li2p(fabs(z));
        return thing; //Why twice?
    }
    else {
        thing = Li2p(z);
    }


    return thing;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Li2p(double z)	{		//! positive arguments only
    //Undefined in the paper.
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double what;
    const double pibsg = 4*atan(1);
    double sum;
    int i;

    sum = 0;
    z = fabs(z);

    if(z<1) {
        for (i=1;i<=50;i++) {
            sum = sum + pow(z,i)/(1*i)*(1*i);
        }
    }
    else {     //(z.ge.1d0) then
        sum = (pibsg*pibsg)/3.0-0.5*(log(z)*log(z));  //AAAHHHHH What did you do Gabe? Must check paper.
        for (i=1;i<=50;i++) {
            sum = sum - 1/(i*i*pow(z,i));
        }
    }
//    what = sum;

//    return what;
    return sum;
}




//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C1beff(const int i, double eta,    particle *chhiggs){
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing;
    int j;
    double sum;
    const double mwbsg = 80.39;
    double* ai;
    ai = (double*)malloc(sizeof(double)*8);
    ai[0]=0.6087;
    ai[1]=0.6957;
    ai[2]=0.2609;
    ai[3]=-0.5217;
    ai[4]=0.4086;
    ai[5]=-0.4230;
    ai[6]=-0.8994;
    ai[7]=0.1456;
    
    double* ei;
    ei = (double*)malloc(sizeof(double)*8);
    ei[0]=5.7064;
    ei[1]=-3.8412;
    ei[2]=0;
    ei[3]=0;
    ei[4]=-1.9043;
    ei[5]=-0.1008;
    ei[6]=0.1216;
    ei[7]=0.0183;
    
    double* fi;
    fi = (double*)malloc(sizeof(double)*8);
    fi[0]= -17.3023;
    fi[1]= 8.5027;
    fi[2]= 4.5508;
    fi[3]= 0.7519;
    fi[4]= 2.0040;
    fi[5]= 0.7476;
    fi[6]= -0.5385;
    fi[7]= 0.0914;
    
    double* gi;
    gi = (double*)malloc(sizeof(double)*8);
    gi[0]= 14.8088;
    gi[1]= -10.8090;
    gi[2]= -0.8740;
    gi[3]= 0.4218;
    gi[4]= -2.9347;
    gi[5]= 0.3971;
    gi[6]= 0.1600;
    gi[7]= 0.0225;
    
    double x = mbart(mwbsg)*mbart(mwbsg)/( mwbsg*mwbsg); //This looks like the run mass of some particle...

    if(i==7) {
        sum = 0;
        for (j=0;j<=7;j++) {
            sum = sum + (ei[j]*eta*E(x)+fi[j]+gi[j]*eta)*pow(eta,ai[j]);
        }

thing = pow(eta,(39/23))*C1Weff(7,chhiggs) + 8/3*( pow(eta,(37/23))-pow(eta,(39/23)) )*C1Weff(8,chhiggs)
            + ( 297664/14283*pow(eta,(16/23)) - 7164416/357075*pow(eta,(14/23))
            +   256868/14283*pow(eta,(37/23)) - 6698884/357075*pow(eta,(39/23)) )*C0W(8,chhiggs)
            + 37208/4761*( pow(eta,(39/23))-pow(eta,(16/23)) )*C0W(7,chhiggs) + sum;

    }
    else {
//  stop 'wtf c1beff'  ??
    }


//c	write(*,*)C1beff
//c	write(*,*)eta**(39d0/23d0)*C1Weff(7)+8d0/3d0*(eta**(37d0/23d0)-eta**(39d0/23d0))*C1Weff(8)
//c	write(*,*)C1beff - (eta**(39d0/23d0)*C1Weff(7)+8d0/3d0*(eta**(37d0/23d0)-eta**(39d0/23d0))*C1Weff(8))
//c	stop

//    if(order='LO') thing = 0;

    free(ai);
    free(ei);
    free(fi);
    free(gi);
    
    return thing;
    }

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0beff(const int i, const double eta,    particle *chhiggs) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    int j;
    double sum, thing;

    double* ai;
    ai = (double*)malloc(sizeof(double)*8);
    ai[0]=0.6087;
    ai[1]=0.6957;
    ai[2]=0.2609;
    ai[3]=-0.5217;
    ai[4]=0.4086;
    ai[5]=-0.4230;
    ai[6]=-0.8994;
    ai[7]=0.1456;
    
    double* hi;
    hi = (double*)malloc(sizeof(double)*8);
    hi[0]= 2.2996;
    hi[1]= -1.0880;
    hi[2]= -0.4286;
    hi[3]= -0.07143;
    hi[4]= -0.6494;
    hi[5]= -0.0380;
    hi[6]= -0.0186;
    hi[7]= -0.0057;
/*
if(i.eq.7) then
sum = 0d0
do j=1,8
sum = sum + eta**ai(j)*hi(j)
enddo

C0beff = eta**(16d0/23d0)*C0W(7)+8d0/3d0*(eta**(14d0/23d0)-eta**(16d0/23d0))*C0W(8) + sum

elseif(i.eq.8) then
C0beff = (C0W(8)+313063d0/363036d0)*eta**(14d0/23d0) -
.	0.9135d0*eta**(0.4086d0) + 0.0873d0*eta**(-0.4230d0)
.	- 0.0571d0*eta**(-0.8994d0) + 0.0209d0*eta**(0.1456d0)
else
C0beff = C0b(i)
endif
 */
    switch(i) {
        case 7: sum = 0;
            for (j=0;j<=7;j++) {
                sum = sum + pow(eta,ai[j])*hi[j]; //This looks like an integral?
            }
            
            thing = pow(eta,(16/23))*C0W(7, chhiggs)+8/3*(pow(eta,(14/23))-pow(eta,(16/23)))*C0W(8, chhiggs)
            + sum;
            break;
        case 8: thing = (C0W(8, chhiggs)+313063/363036)*pow(eta,(14/23)) -
            	0.9135*pow(eta,(0.4086)) + 0.0873*pow(eta,(-0.4230))
            	- 0.0571*pow(eta,(-0.8994)) + 0.0209*pow(eta,(0.1456));
            break;
        default: thing = C0b(i,eta);
            break;
    }
    free(ai);
    free(hi);
    return thing;
}

/*  I need to figure out how to use this information, which is called here.
 c	Arrays for running from mu_W to mu_b
	real*8 ai(8),ei(8),fi(8),gi(8),hi(8)
	data ai/0.6087d0,0.6957d0,0.2609d0,-0.5217d0,0.4086d0,-0.4230d0,-0.8994d0,0.1456d0/
	data ei/5.7064d0,-3.8412d0,0d0,0d0,-1.9043d0,-0.1008d0,0.1216d0,0.0183d0/
	data fi/-17.3023d0,8.5027d0,4.5508d0,0.7519d0,2.0040d0,0.7476d0,-0.5385d0,0.0914d0/
	data gi/14.8088d0,-10.8090d0,-0.8740d0,0.4218d0,-2.9347d0,0.3971d0,0.1600d0,0.0225d0/
	data hi/2.2996d0,-1.0880d0,-0.4286d0,-0.07143d0,-0.6494d0,-0.0380d0,-0.0186d0,-0.0057d0/

 double* ai;
 ai = (double*)malloc(sizeof(double)*8);
 ai[0]=0.6087;
 ai[1]=0.6957;
 ai[2]=0.2609;
 ai[3]=-0.5217;
 ai[4]=0.4086;
 ai[5]=-0.4230;
 ai[6]=-0.8994;
 ai[7]=0.1456;

 double* ei;
 ei = (double*)malloc(sizeof(double)*8);
 ei[0]=5.7064;
 ei[1]=-3.8412;
 ei[2]=0;
 ei[3]=0;
 ei[4]=-1.9043;
 ei[5]=-0.1008;
 ei[6]=0.1216;
 ei[7]=0.0183;
 
 double* fi;
 fi = (double*)malloc(sizeof(double)*8);
 fi[0]= -17.3023;
 fi[1]= 8.5027;
 fi[2]= 4.5508;
 fi[3]= 0.7519;
 fi[4]= 2.0040;
 fi[5]= 0.7476;
 fi[6]= -0.5385;
 fi[7]= 0.0914;
 
 double* gi;
 gi = (double*)malloc(sizeof(double)*8);
 gi[0]= 14.8088;
 gi[1]= -10.8090;
 gi[2]= -0.8740;
 gi[3]= 0.4218;
 gi[4]= -2.9347;
 gi[5]= 0.3971;
 gi[6]= 0.1600;
 gi[7]= 0.0225;
 
 double* hi;
 hi = (double*)malloc(sizeof(double)*8);
 hi[0]= 2.2996;
 hi[1]= -1.0880;
 hi[2]= -0.4286;
 hi[3]= -0.07143;
 hi[4]= -0.6494;
 hi[5]= -0.0380;
 hi[6]= -0.0186;
 hi[7]= -0.0057;
 
 free(ai);
 free(ei);
  free(fi);
  free(gi);
  free(hi);*/

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0b(const int i, double eta) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing;

/*
if(i.eq.1) then
C0b = -eta**(-12d0/23d0) + eta**(6d0/23d0)
elseif(i.eq.2) then
C0b = 1d0/3d0*eta**(-12d0/23d0) + 2d0/3d0*eta**(6d0/23d0)
elseif(i.eq.3) then
C0b = -1d0/27d0*eta**(-12d0/23d0) + 2d0/63d0*eta**(6d0/23d0)
.	- 0.0659d0*eta**0.4086d0 + 0.0595d0*eta**(-0.423d0)
.	- 0.0218d0*eta**(-0.8994d0) + 0.0335d0*eta**(0.1456d0)
elseif(i.eq.4) then
C0b = 1d0/9d0*eta**(-12d0/23d0) + 1d0/21d0*eta**(6d0/23d0)
.	+ 0.0237d0*eta**0.4086d0 - 0.0173d0*eta**(-0.423d0)
.	- 0.1336d0*eta**(-0.8994d0) - 0.0316d0*eta**(0.1456d0)
elseif(i.eq.5) then
C0b = 1d0/108d0*eta**(-12d0/23d0) - 1d0/126d0*eta**(6d0/23d0)
.	+ 0.0094d0*eta**0.4086d0 - 0.0100d0*eta**(-0.423d0)
.	+ 0.0010d0*eta**(-0.8994d0) - 0.0017d0*eta**(0.1456d0)
elseif(i.eq.6) then
C0b = -1d0/36d0*eta**(-12d0/23d0) - 1d0/84d0*eta**(6d0/23d0)
.	+ 0.0108d0*eta**0.4086d0 + 0.0163d0*eta**(-0.423d0)
.	+ 0.0103d0*eta**(-0.8994d0) + 0.0023d0*eta**(0.1456d0)
else
C0b = C0beff(i)

endif
*/

    switch(i){
        case 1: thing = -pow(eta,(-12/23)) + pow(eta,(6/23));
            break;
        case 2: thing = 1/3*pow(eta,(-12/23)) + 2/3*pow(eta,(6/23));
            break;
        case 3: thing = -1/27*pow(eta,(-12/23)) + 2/63*pow(eta,(6/23))
                    - 0.0659*pow(eta,0.4086) + 0.0595*pow(eta,(-0.423))
                - 0.0218*pow(eta,(-0.8994)) + 0.0335*pow(eta,(0.1456));
            break;
        case 4: thing = 1/9*pow(eta,(-12/23)) + 1/21*pow(eta,(6/23))
            	+ 0.0237*pow(eta,0.4086) - 0.0173*pow(eta,(-0.423))
            	- 0.1336*pow(eta,(-0.8994)) - 0.0316*pow(eta,(0.1456));
            break;
        case 5: thing = 1/108*pow(eta,(-12/23)) - 1/126*pow(eta,(6/23))
            	+ 0.0094*pow(eta,0.4086) - 0.0100*pow(eta,(-0.423))
            	+ 0.0010*pow(eta,(-0.8994)) - 0.0017*pow(eta,(0.1456));
            break;
        case 6: thing = -1/36*pow(eta,(-12/23)) - 1/84*pow(eta,(6/23))
            	+ 0.0108*pow(eta,0.4086) + 0.0163*pow(eta,(-0.423))
            	+ 0.0103*pow(eta,(-0.8994)) + 0.0023*pow(eta,(0.1456));
            break;
    }
    
    return thing;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C1Weff(const int i,    particle *chhiggs) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing;
    double Q,x,y;
    const double mwbsg = 80.39;


    Q = mwbsg;
    x = mbart(Q)*mbart(Q)*1.0/(mwbsg*mwbsg);
    y = Q*Q*1.0/(mwbsg*mwbsg);  //Placing this here to shorten the calculation.
//    double x = mbart(mwbsg)*mbart(mwbsg)/( mwbsg*mwbsg); Is this what he meant?

    thing = 0;
/*if(i.eq.1) then
C1Weff = 15d0+6d0*dlog(q**2/mwbsg**2)
return
elseif(i.eq.4) then
C1Weff = E(x) - 2d0/3d0 + 2d0/3d0*dlog(q**2/mwbsg**2) + C1WHeff(i)
return
elseif(i.eq.7) then
C1Weff = G7(x) + Delta7(x)*dlog(q**2/mwbsg**2) + C1WHeff(i)
return
elseif(i.eq.8) then
C1Weff = G8(x) + Delta8(x)*dlog(q**2/mwbsg**2) + C1WHeff(i)
return
endif
*/

    //These look like log(1) when q=mwbsg, as assigned above...

    switch(i){
        case 1: thing = 15+6*log(y);
            break;
        case 4: thing = E(x) - 2/3 + 2/3*log(y) + C1WHeff(i, chhiggs);
            break;
        case 7: thing = G7(x) + Delta7(x)*log(y) + C1WHeff(i, chhiggs);
            break;
        case 8: thing = G8(x) + Delta8(x)*log(y) + C1WHeff(i, chhiggs);
            break;
        default: break;
    }

    return thing;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double E(const double x){
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(-18+11*x+x*x)/(12*(x-1)*(x-1)*(x-1))
            + x*x*(15-16*x+4*x*x)/(6*(x-1)*(x-1)*(x-1)*(x-1))*log(x)
            - 2/3*log(x);


    return what;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G7(const double x) {
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;
    double term1,term2,term3,term4;


    term1 = (-16*x*x*x*x-122*x*x*x+80*x*x-8*x)/(9*pow((x-1),4))*Li2(1-1/x);
    term2 = (6*x*x*x*x+46*x*x*x-28*x*x)/(3*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*pow(log(x),2);
    term3 = (-102*pow(x,5)-588*x*x*x*x-2262*x*x*x+3244*x*x-1364*x+208)/(81*pow((x-1),5))*log(x);
    term4 = (1646*x*x*x*x+12205*x*x*x-10740*x*x+2509*x-436)/(486*(x-1)*(x-1)*(x-1)*(x-1));

    what = term1 + term2 + term3 + term4;


    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G8(const double x) {
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;
    double term1,term2,term3,term4;


    term1 = (-4*x*x*x*x+40*x*x*x+41*x*x+x)/(6*(x-1)*(x-1)*(x-1)*(x-1))*Li2(1-1/x);
    term2 = (-17*x*x*x-31*x*x)/(2*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*pow(log(x),2);
    term3 = (-210*x*x*x*x*x+1086*x*x*x*x+4893*x*x*x+2857*x*x-1994*x+280)/(216*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*log(x);
    term4 = (737*x*x*x*x-14102*x*x*x-28209*x*x+610*x-508)/(1296*(x-1)*(x-1)*(x-1)*(x-1));

    what = term1 + term2 + term3 + term4;

    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta7(const double x) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;

    whatever = (208-1111*x+1086*x*x+383*x*x*x+82*x*x*x*x)/(81*(x-1)*(x-1)*(x-1)*(x-1))
            + (2*x*x*(14-23*x-3*x*x))/(3*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta8(const double x) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;


    whatever = (140-902*x-1509*x*x-398*x*x*x+77*x*x*x*x)/(108*(x-1)*(x-1)*(x-1)*(x-1))
        + (x*x*(31+17*x))/(2*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0Weff(const int i,    particle *chhiggs) {
    //Why does this function exist?
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing = C0W(i, chhiggs);

    return thing;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0W(const int i,    particle *chhiggs) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     i   input   which coefficient is being calculated.
     */

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    double thing;
    double Q,x;
    const double mwbsg = 80.39;

    Q = mwbsg;

    x = mbart(Q)*mbart(Q)*1.0/(mwbsg*mwbsg);

//c	write(*,*)x,mbart(Q),mw,q
//c	stop
 /*   if(i=2) {
        thing = 1;
        return thing;
    }

    if(i=6) {
        thing = 0;
        return thing;
    }

    if(i=7) {
        thing = F17(x);
    }

    if(i=8) {
        thing = F18(x);
    }

    thing += C0WHiggs(i);
*/
    
    switch(i){
        case 2: thing = 1; break;
        case 6: thing = 0; break;
        case 7: thing = F17(x); break;
        case 8: thing = F18(x); break;
    }
  
    thing += C0WHiggs(i,chhiggs);
    
    return thing;
    }

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F17(const double x){
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(7-5*x-8*x*x)/(24*(x-1)*(x-1)*(x-1))
            + x*x*(3*x-2)/(4*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F18(const double x) {
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(2+5*x-x*x)/(8*(x-1)*(x-1)*(x-1))
            - 3*x*x/(4*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0WHiggs(int i,    particle *chhiggs) {
    /* Correction to the SM contributions by the charged Higgs at LO to the effective Wilson coefficients.
     i   input   which coefficient correction is being calculated.
     
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    
    int j;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    const double mwbsg = 80.39;
    double sum;
    double x;
    const int nhiggs = sizeof(chhiggs)/sizeof(chhiggs[0]);
    double Au,Ad;

    sum = 0;
    for (j=0;j<nhiggs;j++){
        x = mbart(mwbsg)*mbart(mwbsg)/(chhiggs[j].mass*chhiggs[j].mass);
        Au = chhiggs[j].Y_t/(sqrt(2.0)*mtbsg/246.0);
        Ad = chhiggs[j].Y_b/(sqrt(2.0)*mbbsg/246.0);
        if(i==7) sum = sum + Au*Au/3*F17(x) - Au*Ad*F27(x);
        if(i==8) sum = sum + Au*Au/3*F18(x) - Au*Ad*F28(x);
    }

//C0WHiggs = sum
    return sum;
    }

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F27(const double x) {
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;

    whatever = x*(3-5*x)/(12*(x-1)*(x-1)) + x*(3*x-2)/(6*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F28(const double x) {
    /* Correction to the SM contributions by the charged Higgs at LO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever = x*(3-x)/(4*(x-1)*(x-1)) - x/(2*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C1WHeff(int i,    particle *chhiggs){
// Summation of each effective Wilson coefficient (i) over all Higgses (j).
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    int j;
    const double mwbsg = 80.39;
    const double mbbsg = 4.75;
    const double mtbsg = 175;
    double x;
    double sum =0;
    const int nhiggs = sizeof(chhiggs)/sizeof(chhiggs[0]);
    double Au,Ad;


//    whatever = 0;
    return sum; //Why is this here?  It returns 0, then a double, if there's a double? Why not wait?
//    sum = 0;
    
    for (j=1;nhiggs;j++){ //I believe that this is where the function starts adding up contributions from multiple Higgs, hence the index "nhiggs".
        x = mbart(mwbsg)*mbart(mwbsg)/(chhiggs[j].mass*chhiggs[j].mass);
        Au = chhiggs[j].Y_t/(sqrt(2.0)*mtbsg/246.0);
        Ad = chhiggs[j].Y_b/(sqrt(2.0)*mbbsg/246.0);

        if(i==4) sum = sum + EH(x,Au,Ad);
        if(i==7) sum = sum + G7H(x,Au,Ad)
            + Delta7H(x,Au,Ad)*log(mwbsg*mwbsg/(chhiggs[j].mass*chhiggs[j].mass))
            - 4/9*EH(x,Au,Ad);
        if(i==8) sum = sum + G8H(x,Au,Ad)
            + Delta8H(x,Au,Ad)*log(mwbsg*mwbsg/(chhiggs[j].mass*chhiggs[j].mass))
            - 1/6*EH(x,Au,Ad);
    }

//    whatever = sum;

    return sum;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double EH(const double x,const double xAu,const double xAd){
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
//real*8 xAu,xAd,x

    whatever = xAu*xAu*(x*(16-29*x+7*x*x)/(36*pow((x-1),3))+x*(3*x-2)/(6*pow((x-1),4))*log(x));


    return whatever;
}



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G7H(const double x,const double xAu,const double xAd){
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
    double term1,term2;


    term1 = xAd*xAu*4/3*x*(
                           4*(-3+7*x-2*x*x)/(3*pow((x-1),3))*Li2(1-1/x)+
                           (8-14*x-3*x*x)/(3*pow((x-1),4))*pow(log(x),2)+
                           (2*(-3-x+12*x*x-2*x*x*x)/(3*pow((x-1),4)))*log(x)+
                           (7-13*x+2*x*x)/pow((x-1),3)
                           );
    term2 = xAu*xAu*2/9*x*(
                          x*(18-37*x+8*x*x)/pow((x-1),4)*Li2(1-1/x) +
                          x*(-14+23*x+3*x*x)/pow((x-1),5)*pow(log(x),2) +
                          (-50+251*x-174*x*x-192*x*x*x+21*x*x*x*x)/(9*pow((x-1),5))*log(x) +
                          (797-5436*x+7569*x*x-1202*x*x*x)/(108*pow((x-1),4))
                           );

    whatever = term1 + term2;

    return whatever;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta7H(const double x,const double xAu,const double xAd){
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
//    double xAu,xAd,x;

    whatever = xAu*xAd*2/9*x*(
                              (21-47*x+8*x*x)/pow((x-1),3)
                              +2*(-8+14*x+3*x*x)/pow((x-1),4)*log(x)
                              )
                +xAu*xAu*2/9*x*(
                                (-31-18*x+135*x*x-14*x*x*x)/(6*pow((x-1),4))
                                +x*(14-23*x-3*x*x)/pow((x-1),5)*log(x)
                                );

    return whatever;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G8H(const double x,const double xAu,const double xAd){
    /* Correction to the SM contributions to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//real*8 xAu,xAd,x
    double whatever;
    double term1,term2;


    term1 = xAd*xAu*1/3*x*( (-36+25*x-17*x*x)/(2*pow((x-1),3))*Li2(1-1/x)
                           +(19+17*x)/(pow((x-1),4))*pow(log(x),2)
                           +(-3-187*x+12*x*x-14*x*x*x)/(4*pow((x-1),4))*log(x)
                           +3*(143-44*x+29*x*x)/(8*pow((x-1),3)) );
    term2 = xAu*xAu*1/6*x*( x*(30-17*x+13*x*x)/pow((x-1),4)*Li2(1-1/x)
                           -x*(31+17*x)/pow((x-1),5)*pow(log(x),2)
                           +(-226+817*x+1353*x*x+318*x*x*x+42*x*x*x*x)/(36*pow((x-1),5))*log(x)
                           +(1130-18153*x+7650*x*x-4451*x*x*x)/(216*pow((x-1),4)) );

    whatever = term1 + term2;



    return whatever;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta8H(const double x,const double xAu,const double xAd){
/* Correction to the SM contributions to the effective Wilson coefficients.
    x   input   running mass of the top quark at the W mass, divided by the W mass.
    xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
    xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
*/
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double xAu,xAd,x;
    double what;


    what = xAu*xAd*1/3*x*( (81-16*x+7*x*x)/(2*pow((x-1),3))-(19+17*x)/pow((x-1),4)*log(x) ) + xAu*xAu*1/6*x*( (-38-261*x+18*x*x-7*x*x*x)/(6*pow((x-1),4))+x*(31+17*x)/pow((x-1),5)*log(x) );

    return what;
}

