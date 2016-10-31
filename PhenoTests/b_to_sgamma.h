      //
//  b_to_sgamma.h
//  
//
//  Created by Valerie Plaus on 7/15/16.
//
//

#ifndef b_to_sgamma_h
#define b_to_sgamma_h

#include <stdio.h>
#include "Particle_struct.h"
void bsg_nlo(  particle *chhiggs,double *BF,double *pull);

//double bsg_nlo(double *xAu,double *xAd,double *xMH,double xnhiggs,double BF,double pull){
    /*c
     c	Calculates the B -> X_S + gamma BF from the Standard model and 2HDM
     c
     c	M. Ciuchini, G. Degrassi, P. Gambino, and G. F. Giudice, Nucl. Phys. B527, 21 (1998).
     c
     c
     c	xAu - Charged Higgs PL coupling (proportional to mu) (Parity - Left)
     c	xAd - Charged Higgs PR coupling (proportional to md) (Parity - Right)
     c	xMH - Charged Higgs mass
     c	xnhiggs - Higgs index
     c
     c	BF - the result BF(B -> X_S + gamma + X)
     c	pull - the number of std. dev. away from experiment
     c		(taking into account experimental uncertainty and 10% theoretical uncertainty)
     c
     c
     c
     c	G. Shaughnessy 5/3/2012
     c
     c
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

double DtermSq(double Q,   particle *chhiggs) ;
    /* Now we're actually getting to a calculation.
     Q  input   energy scale of the calculation.
     */

double gam0eff(int i,int j) ;
    //Why does he have J as an input here?  He doesn't use it.

double Aterm(double Q, double eta,  particle *chhiggs) ;
    /* Correction for low energy bremsstrahlang.
     Q  input   energy of calculation
     */

double mbart(double Q);
    /* Correction for low energy bremsstrahlang.
     Q  input   energy of calculation
     */

double fij(int i,int j) ;
    /* List of the matrix elements for the possible interactions.  Each element listed below corresponds to a particular particle interaction and Wilson coefficient.
     i  input   column index for matrix element
     j  input   row index for matrix element
     */
    
double mbrun(const double Q) ;
    /* The running value of the bottom mass. Eqn 42
     Q   input   energy scale at which the calling function needs the bottom mass.
     */

double alphasbsg(const double Q);
    /* The running value of the strong coupling constant. Eqn 42
     Q   input   energy scale at which the calling function needs the strong coupling.
     */
    
double kap(const double z) ;
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double f(const double z) ;
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double h(const double z) ;
    /* Phase space calculations for the SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double Li2(const double z) ;			//!either positive or negative arguments allowed
    //Undefined in the paper.

double Li2p(const double z)	;		//! positive arguments only
    //Undefined in the paper.

double C1beff(const int i, double eta,   particle *chhiggs);
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
    
double C0beff(const int i, const double eta,   particle chhiggs) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double C0b(const int i, double eta) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double C1Weff(const int i,   particle *chhiggs) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double E(const double x);
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double G7(const double x) ;
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double G8(const double x) ;
    /* SM contributions at LO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double Delta7(const double x) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */
    
double Delta8(const double x) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   ratio of top mass to W mass scale.
     */

double C0Weff(const int i,   particle chhiggs) ;
    //Why does this function exist?

double C0W(const int i,   particle *chhiggs) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     i   input   which coefficient is being calculated.
     */

double F17(const double x);
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */

double F18(const double x) ;
    /* SM contributions at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */

double C0WHiggs(int i,   particle *chhiggs) ;
    /* Correction to the SM contributions by the charged Higgs at LO to the effective Wilson coefficients.
     i   input   which coefficient correction is being calculated.
     
     Need to pass scaled Yukawas Au and Ad.
     */
    
double F27(const double x) ;
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     */

double F28(const double x);
    /* Correction to the SM contributions by the charged Higgs at LO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     coupling.
     */
    
double C1WHeff(int i,   particle *chhiggs);
    // Summation of each effective Wilson coefficient (i) over all Higgses (j).

double EH(const double x,const double xAu,const double xAd);
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */

double G7H(const double x,const double xAu,const double xAd);
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */

double Delta7H(const double x,const double xAu,const double xAd);
    /* Correction to the SM contributions by the charged Higgs at NLO to the effective Wilson coefficients.
     x   input   running mass of the top quark at the W mass, divided by the W mass.
     xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
     xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
     */
double G8H(const double x,const double xAu,const double xAd);
    /* Correction to the SM contributions to the effective Wilson coefficients.
    x   input   running mass of the top quark at the W mass, divided by the W mass.
    xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
    xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
    */

double Delta8H(const double x,const double xAu,const double xAd);
    /* Correction to the SM contributions to the effective Wilson coefficients.
    x   input   running mass of the top quark at the W mass, divided by the W mass.
    xAu input   Effective up-type (PL) Yukawa coupling, scaled by SM up-type Yukawa coupling.
    xAd input   Effective down-type (PR) Yukawa coupling, scaled by SM down-type Yukawa coupling.
    */


#endif /* b_to_sgamma_h */
