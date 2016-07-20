//
//  Phenotests.c
//  
//
//  Created by Valerie Plaus on 7/12/16.
//
//



#include "Phenotests.h"
#include "Particle_Struct.h"


//----------------------------------------------------------c
int coupling_perturb(double Yb,double Yt,double tanbeta1){
/*c-----------------------------------------------------------c
 c                                                           c
 c This subroutine calculates the ratio for the tbH+ couplingc
 c between the Higgs-stew model and the MSSM.  To get the    c
 c MSSM parameters we use \mu = \lambda_1*s1 and b = 100*mu. c
 c \tan(\beta) is set so that the MSSM charged higgs matches c
 c the lightest charged Higgs mass of Higgs-stew.            c
 c                                                           c
 
    Originally created by Mat McCaskey, UW Madison, 2011
    Adapted for C by Valerie Plaus, Wittenberg University, 2016.
 c-----------------------------------------------------------*/

    int higgscheck;

// makes sure that the yukawas that we have a not out of the realm of perturbativity
    if ((Yb>=sqrt(8.0))||(Yt>=sqrt(8.0))) {
        higgscheck = 0;
    }
// makes sure that the ratios of the vacuum expectation values are not outside expected values
    if (tanbeta1>=55.0||tanbeta1<==0.10){
    higgscheck = 0;
    }
    
/* This makes sure that the vevs themselves are not too small (decoupling), which should in the randomization
 if (v1<=2.0||v2<=2.0||v3<=2.0||v4<=2.0) then {
 higgscheck = 0;} */


//      printf ("%d/t%d/t", Yb-sqrt(8.0), Yt-sqrt(8.0));

    return higgscheck;
}

//-----------------------------------------------------------c
int hzz_check(Particle *A, double *vevs){
/* This subroutine checks that the point does not violate the
 Large Electron-Positron Collider (LEP) constraints.  The
 constraints are expressed as a linear interpolation of a data
 set imported during "initialize" or "openfiles" and apply to
 the cross sections calcuated using eigenvalues and eigenvectors
 of the cp even higgs mass matrix mh(i) and heigvec(i,j).
 
 Personally, I think this would be more fun as a single-Higgs check, repeatedly called in the main program, instead of checking them all in here.
 
 Originally created by Mat McCaskey, UW Madison, 2011
 Adapted for C by Valerie Plaus, Wittenberg University, 2016.

 c-----------------------------------------------------------*/
/*  A           input   The particle's attributes.
 
    vevs        input   An array containing the vacuum expectation values (vevs).
 */

// parameters in this subroutine only
    int hzzcheck; //We could change this to a logical, since it just needs to be true or false.  If changed, the function would also need to change from an int to a logical.
    double precision hzzcutoff; //Place holder for calculating the interpolation of the LEP data at the value of the Higgs mass.

// Lep constraint
        if (A.mass>=2.0*4.7) {
            hzz_ratio = (Dot_Prod(*A.eigenvec,*veves))*(Dot_Prod(*A.eigenvec,*veves))/(246.0*246.0)*(BFcpe[i][0]/BFSM[i][0]);//!(1-BFinv_decay(i));
        hzzcutoff = lineint(hmass[mh[i]],zzhmass,zzhcoup,217);
        /*        else if (hmass(mh(i))>=2.0*1.78) then
         hzz_ratio(i) = (v1*heigvec(1,mh(i)) + v2*heigvec(2,mh(i)) + v3*heigvec(3,mh(i)) +  v4*heigvec(4,mh(i)))**2/246.d0**2*(BFcpe(i,3)/BFSM(i,3))! need to import tau tau lep constraint. this array is a guess
         hzzcutoff = lineint(hmass(mh(i)),zzhmass,zzhcouptau,217)
         printf (BFcpe(i,3),BFSM(i,3)) */
    }
    else{
        hzzcheck = 0;
    }
    
    //        if (hmass(mh(i))<=110.0 ) then
    if (mass<=114.0&&hzz_ratio>=hzzcutoff) {
        hzzcheck = 0;
    }
    
    //      if (hmass(mh(i))<=15){
    //      printf (hmass(mh(i)), hzz_ratio(i), hzzcutoff,hzzcheck);
    //  }
    
    return hzzcheck;
    }



//c-----------------------------------------------------------c
    void cross_sections(*A, *vevs){
/*c-----------------------------------------------------------c
c                                                           c
c This subroutine calculates the ratio fo the tbH+ coupling c
c between the Higgs-nkr model and the MSSM.  To get the    c
c MSSM parameters we use \mu = \lambda_1*s1 and b = 100*mu. c
c \tan(\beta) is set so that the MSSM charged higgs matches c
c the lightest charged Higgs mass of Higgs-stew.            c
c                                                           c
c-----------------------------------------------------------c*/

    double precision geff;  //intermediate memory places to make the code more readable

//c tb production vertex, unscaled
// The loop starts at 1 because the lightest charged Higgs eigenvalue should be a 0 (a Goldstone).
        if (A.mass<=3000){
            tbprod = (A.Y_b)*(A.Y_b)+(A.Y_t)*(A.Y_t);
            csprod = (A.Y_s)*(A.Y_s)+(A.Y_c)*(A.Y_c);
            mch_tb_prod = tbprod*lineint(A.mass,MHc_MSSM,cross_MSSM,(sizeof(cross_MSSM)/sizeof(cross_MSSM[0])));
            mch_cs_prod = csprod*lineint(A.mass,MHc_MSSM,cross_MSSM,(sizeof(cross_MSSM)/sizeof(cross_MSSM[0])));
        }
        else {
            tbprod = 0.0;
            csprod = 0.0;
            mch_tb_prod = 0.0;
            mch_cs_prod = 0.0;
        }

//c production rates for h_0 as ratios relative to the SM
        geff = g2*(Dot_Prod(*vevs,*A.eigenvec));
        ggprod[i] = h2glgl(A.mass,A.Y_t,A.Y_b); //! Gabe's code gives prod ratio
        gagaprod[i] = h2gaga(A.mass,A.Y_t,A.Y_b,geff); //! Gabe's code gives prod ratio
//c          mh_gg_prod(i) = ggprod(i)*?
//c          mh_gaga_prod(i) = gagaprod(i)*?
        geff = 0;


//c cross section for h_0 as ratios relative to the SM
// The second indices for the BF__ arrays correspond to the correct column in an imported file.
// 3rd Column: Higgs decay to 2 tau particles           [2]
// 5th Column: Higgs decay to 2 W particles             [4]
// 6th Column: Higgs decay to 2 photons (gamma-gamma)   [5]
// 7th Column: Higgs decay to 4 leptons                 [6]
// 8th Column: Higgs decay to 2 leptons, 2 neutrinos    [7]
        if (A.mass>=100.0){
            xsgg2ww[i]=ggprod[i]* (BFcpe[i][4]/BFSM[i][4]);
            xsgg2gaga[i]=ggprod[i]* (BFcpe[i][5]/BFSM[i][5]);
            xsgg24l[i]=ggprod[i]* (BFcpe[i][6]/BFSM[i][6]);
            xsgg22l[i]=ggprod[i]* (BFcpe[i][7]/BFSM[i][7]);

//c       if (hmass(mh(i)).gt.80d0.and.hmass(mh(i)).lt.600d0) then
            gghprod[i] = ggprod[i]*lineint(hmass[mh[i]],masshgg,hggprod,(sizeof(hggprod)/sizeof(hggprod[0])));
            xsgg2tautau[i]=gghprod[i]*(BFcpe[i][2]); //! the plot limit here is production x BF
/*c           ATLAS-CONF-2011-132
c          else
c            hggprod(i) = 0d0
c            xsgg2tautau(i) = 0
c          endif
*/
        }
        else {
            xsgg24l[i]=0.0;
            xsgg22l[i]=0.0;
            xsgg2ww[i]=0.0;
            xsgg2gaga[i]=0.0;
        }
/*
c output for testing purposes
c        write(*,*) 'yukawas: Y1 =',Y1,', Y2 =',Y2,' Y2p =',Y2p
c        write(*,*) 'MSSM parameters: mu =',mu_mssm,', tan(beta) =',dtan(beta_mssm)
c        write(*,*) 'MSSM yukawas: Y1 =',Y1_mssm,', Y2 =',Y2_mssm
c        write(*,*) 'Coupling ratio = ',ratio
c      endif 
*/

        return;
}
    
///////////////////////////////////////////////////////////////////////
//c-----------------------------------------------------------c
    int atlas_check_preJuly4(*hmass, *mh){
//c-----------------------------------------------------------c


//c parameters in this subroutine only
        double hwwcutoff,hzz2lcutoff,hzz4lcutoff,hgagacutoff,htautaucutoff;
        double atlaswwcheck,atlas2lcheck,atlas4lcheck,atlasgagacheck,atlastautaucheck;


//c Atlas suggested excess
            hwwcutoff=0.0;
            hzz2lcutoff=0.0;
            hzz4lcutoff=0.0;
            hgagacutoff=0.0;
            htautaucutoff=0.0;
            if (A.mass > 110.305 && A.mass < 600){
                hwwcutoff = lineint(A.mass,masshww,hwwlimits,79);
//c	    print *,'ww',i,hmass(mh(i)),xsgg2ww(i),hwwcutoff
                if (xsgg2ww >= hwwcutoff) {
                    atlaswwcheck = 0;
//c          else
//c	    write(16,*) i,hmass(mh(i)),xsgg2ww(i),hwwcutoff
                }
            }
            if (A.mass > 200.345 && A.mass < 600){
                hzz2lcutoff = lineint(A.mass,massh2l2nu,hZZlnulimits,50);
//c	    print *,'zz',i,hmass(mh(i)),xsgg22l(i),hzz2lcutoff
                if (xsgg22l >= hzz2lcutoff) {
                    atlas2lcheck = 0;
                }
            }
            if (A.mass > 110.152 && A.mass < 598.344) {
                hzz4lcutoff = lineint(A.mass,massh4l,hZZ4llimits,130);
//c	    print *,'4l',i,hmass(mh(i)),xsgg24l(i),hzz4lcutoff
                if (xsgg24l >= hzz4lcutoff){
                    atlas4lcheck = 0;
                }
            }
            if (A.mass > 110.457 && A.mass < 150.102) {
                hgagacutoff = lineint(hmass[mh[i]],masshgaga,hgagalimits,67);
//c	    print *,'gaga',i,hmass(mh(i)),xsgg2gaga(i),hgagacutoff
                if (xsgg2gaga >= hgagacutoff){
                    atlasgagacheck = 0;
                }
            }
            if (A.mass > 100 && A.mass < 600) {
                htautaucutoff = lineint(A.mass,masshtautau,htautaulimits,16);
//c	    print *,'tata',i,hmass(mh(i)),xsgg2tautau(i),htautaucutoff
                if (xsgg2tautau >= htautaucutoff) {
                    atlastautaucheck = 0;
                }
        }

//c	  if(atlaswwcheck*atlas2lcheck*atlas4lcheck*atlasgagacheck*atlastautaucheck=0) {
//c	   printf ("failed");
//c    }
//c	  else {
//c	   printf ("passed");
//c	   if(hmass[mh[1]]>=114){pause
//c     }
//c    }

        return atlaswwcheck*atlas2lcheck*atlas4lcheck*atlasgagacheck*atlastautaucheck;
    }

    
// There's no need for two copies for the pre- and post-July 4 data.
    
/*c-----------------------------------------------------------c
subroutine atlas_check_postJuly4()
c-----------------------------------------------------------c
implicit none
include 'manyhiggs.inc'

c parameters in this subroutine only
double precision hwwcutoff,hzz2lcutoff,hzz4lcutoff,hgagacutoff,htautaucutoff


c Atlas suggested excess
do i=1,4
hwwcutoff=0d0
hzz2lcutoff=0d0
hzz4lcutoff=0d0
hgagacutoff=0d0
htautaucutoff=0d0
if (hmass(mh(i)).gt.125.966d0.and.hmass(mh(i)).lt.599.698d0) then
hwwcutoff = lineint(hmass(mh(i)),masshwwJuly4,hwwlimitsJuly4,64)
c	    print *,'ww',i,hmass(mh(i)),xsgg2ww(i),hwwcutoff
if (xsgg2ww(i) .ge. hwwcutoff) then
atlaswwcheck = 0
endif
endif
if (hmass(mh(i)).gt.293.247d0.and.hmass(mh(i)).lt.599.698d0) then
hzz2lcutoff = lineint(hmass(mh(i)),massh2l2nuJuly4,hZZlnulimitsJuly4,49)
c	    print *,'zz',i,hmass(mh(i)),xsgg22l(i),hzz2lcutoff
if (xsgg22l(i) .ge. hzz2lcutoff) then
atlas2lcheck = 0
endif
endif
if (hmass(mh(i)).gt.125.966d0.and.hmass(mh(i)).lt.600.000d0) then
hzz4lcutoff = lineint(hmass(mh(i)),massh4lJuly4,hZZ4llimitsJuly4,179)
c	    print *,'4l',i,hmass(mh(i)),xsgg24l(i),hzz4lcutoff
if (xsgg24l(i) .ge. hzz4lcutoff) then
atlas4lcheck = 0
endif
endif
if (hmass(mh(i)).gt.126.268d0.and.hmass(mh(i)).lt.212.923d0) then
hgagacutoff = lineint(hmass(mh(i)),masshgagaJuly4,hgagalimitsJuly4,60)
c	    print *,'gaga',i,hmass(mh(i)),xsgg2gaga(i),hgagacutoff
if (xsgg2gaga(i) .ge. hgagacutoff) then
atlasgagacheck = 0
endif
endif
if (hmass(mh(i)).gt.100d0.and.hmass(mh(i)).lt.600d0) then
htautaucutoff = lineint(hmass(mh(i)),masshtautau,htautaulimits,16)
c	    print *,'tata',i,hmass(mh(i)),xsgg2tautau(i),htautaucutoff
if (xsgg2tautau(i) .ge. htautaucutoff) then
atlastautaucheck = 0
endif
endif
enddo

c	  if(atlaswwcheck*atlas2lcheck*atlas4lcheck*atlasgagacheck*atlastautaucheck.eq.0) then
c	   print *,'failed'
c	  else
c	   print *,'passed'
c	   if(hmass(mh(1)).ge.114d0) pause
c	  endif

return
end
*/
//c-----------------------------------------------------------c
    int happypoint(){
/*c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c*/

//logical whscan
//common/scantype/whscan

        if (BFt[1]*BFch[1][3] > 0.04) { //top to charm
            happy1 = 0;
        }
        if (BFt[1]*BFch[1][1] > 0.14)  {
            happy2 = 0;
        }

        if(whscan) {
            if (BFch[1][2]*mch_tb_prod[1] < 3) {
                happy3 = 0;
            }
            if (BFch[2][4]*mch_tb_prod[2] < 3) {
                happy3 = 0;
            }
        }
    
//c	if(hmass(mh(1)).gt.114d0) happy3 = 0
//c	if (chmass(mch(2)).gt.600d0) happy3 = 0

        return happy1*happy2*happy3;
    }


//c-----------------------------------------------------------c
    int bsgam_1only(*A){
/*c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c*/

        double Au,Ad,mchtemp,chBF1;//!,chpull;

//c bsg_nlo subroutine

            mchtemp = A.mass;
            Au = (A.Y_t)/(sqrt(2.0)*mt/246); //This scales the calculated Yukawa coupling by the SM coupling.
            Ad = (A.Y_b)/(sqrt(2.0)*mb/246);
//c        write(*,*) mchtemp(i),Yteff(i),Ybeff(i)
       

        bsg_nlo(Au,Ad,mchtemp,1,chBF1,chpull1)
//c        call bsg_nlo(Au,Ad,mchtemp,3,chBF,chpull)
//c make sure the effective Y's are scaled by the SM coupling, which they now are. May 23, 2012

    if(chpull1 >= (1.96))  {
        bsgcheck1 = 0;
    }

//c        if (higgscheck*chiggscheck*cpoddcheck*charginocheck*neutralinocheck.gt.0) then
//c        write(*,*) chpull, bsgcheck
//c        endif

         //For the charged Higgses, we only check the massive states.
            mchtemp = 0.0;
            Au = 0.0;
            Ad = 0.0;


        return bsgcheck1;
    }


//c-----------------------------------------------------------c
//    int bsgam(){
/*c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c*/
/*
    double Au[3],Ad[3],Yteff[3],Ybeff[3],mchtemp[3],chBF3;//!,chpull

//c bsg_nlo subroutine
    for (i=1;i<=3;i++) {
        mchtemp[i] = chmass[mch[i]];
        Yteff[i] = (Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]]);
        Ybeff[i] = Yb*cheigvec[1][mch[i]];
        Au[i] = (Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]])/(sqrt(2.0)*mt/246.0);
        Ad[i] = Yb*cheigvec[1][mch[i]]/(sqrt(2.0)*mb/246.0);
//c        write(*,*) mchtemp(i),Yteff(i),Ybeff(i)
    }

//c        call bsg_nlo(Au(1),Ad(1),mchtemp(1),1,chBF,chpull)
        bsg_nlo(Au,Ad,mchtemp,3,chBF3,chpull3);
//c make sure the effective Y's are scaled by the SM coupling, which they now are. May 23, 2012

    if(chpull3 >= (1.96)) {
        bsgcheck3 = 0;
    }

//c        if (higgscheck*chiggscheck*cpoddcheck*charginocheck*neutralinocheck.gt.0) then
//c        write(*,*) chpull, bsgcheck
//c        endif

    for (i=1;i<=3;i++) {
        mchtemp[i] = 0.0;
        Yteff[i] = 0.0;
        Ybeff[i] = 0.0;
        Au[i] = 0.0;
        Ad[i] = 0.0;
    }


    return bsgcheck3;
    }
*/
