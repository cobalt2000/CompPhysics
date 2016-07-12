//
//  Phenotests.c
//  
//
//  Created by Valerie Plaus on 7/12/16.
//
//



#include "Phenotests.h"


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
int hzz_check(double **heigenvec, double *hmass, double *mh, double *vevs){
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
/*  heigenvec   input   A matrix that contains the arrays of eigenvectors from the 
                        solved CP even Higgs mass matrix.  Their order corresponds to
                        the order of the eigenvalues.
 
    hmass       input   An array that contains the eigenvalues (or masses) of the CP 
                        even Higgses.  These are in no particular order.
 
    mh          input   An array used to order the eigenvalues from smallest to largest.
 */

// parameters in this subroutine only
    int hzzcheck; //We could change this to a logical, since it just needs to be true or false.  If changed, the function would also need to change from an int to a logical.
    double precision hzzcutoff; //Place holder for calculating the interpolation of the LEP data at the value of the Higgs mass.

// Lep constraint
    for (i=0;i<4;i++){
        if (hmass(mh(i))>=2.0*4.7) {
            hzz_ratio(i) = (veves[1]*heigvec(1,mh(i)) + vevs[2]*heigvec(2,mh(i)) + vevs[3]*heigvec(3,mh(i)) + vevs[4]*heigvec(4,mh(i)))*(veves[1]*heigvec(1,mh(i)) + vevs[2]*heigvec(2,mh(i)) + vevs[3]*heigvec(3,mh(i)) + vevs[4]*heigvec(4,mh(i)))/(246.0*246.0)*(BFcpe(i,1)/BFSM(i,1));//!(1-BFinv_decay(i));
        hzzcutoff = lineint(hmass(mh(i)),zzhmass,zzhcoup,217);
        /*        else if (hmass(mh(i))>=2.0*1.78) then
         hzz_ratio(i) = (v1*heigvec(1,mh(i)) + v2*heigvec(2,mh(i)) + v3*heigvec(3,mh(i)) +  v4*heigvec(4,mh(i)))**2/246.d0**2*(BFcpe(i,3)/BFSM(i,3))! need to import tau tau lep constraint. this array is a guess
         hzzcutoff = lineint(hmass(mh(i)),zzhmass,zzhcouptau,217)
         printf (BFcpe(i,3),BFSM(i,3)) */
    }
    else{
        hzzcheck = 0;}
    
    //        if (hmass(mh(i))<=110.0 ) then
    if (hmass(mh(i))<=114.0&&hzz_ratio(i)>=hzzcutoff) then{
        hzzcheck = 0;}
    
    //      if (hmass(mh(i)).le.15d0) then
    //      printf (hmass(mh(i)), hzz_ratio(i), hzzcutoff,hzzcheck);
    //
    
    return hzzcheck;
    }



c-----------------------------------------------------------c
subroutine cross_sections()
c-----------------------------------------------------------c
c                                                           c
c This subroutine calculates the ratio fo the tbH+ coupling c
c between the Higgs-nkr model and the MSSM.  To get the    c
c MSSM parameters we use \mu = \lambda_1*s1 and b = 100*mu. c
c \tan(\beta) is set so that the MSSM charged higgs matches c
c the lightest charged Higgs mass of Higgs-stew.            c
c                                                           c
c-----------------------------------------------------------c
implicit none
include 'manyhiggs.inc'
double precision Yteff,Ybeff,geff

c tb production vertex, unscaled
do i=1,3
if (chmass(mch(i)).le.3d3) then
tbprod(i) = (Yb*cheigvec(1,mch(i)))**2 + (Yt*cheigvec(2,mch(i)) + Ytp*cheigvec(4,mch(i)))**2
csprod(i) = (Ys*cheigvec(1,mch(i)))**2 + (Yc*cheigvec(2,mch(i)) + Ycp*cheigvec(4,mch(i)))**2
mch_tb_prod(i) = tbprod(i)*lineint(chmass(mch(i)),MHc_MSSM,cross_MSSM,71)
mch_cs_prod(i) = csprod(i)*lineint(chmass(mch(i)),MHc_MSSM,cross_MSSM,71)
else
tbprod(i) = 0d0
csprod(i) = 0d0
mch_tb_prod(i) = 0d0
mch_cs_prod(i) = 0d0
endif
enddo

c production rates for h_0 as ratios relative to the SM
do i=1,4
Yteff = (Yt*heigvec(2,mh(i)) + Ytp*heigvec(4,mh(i)))
Ybeff = Yb*heigvec(1,mh(i))
geff = g2*(heigvec(1,mh(i))*v1+heigvec(2,mh(i))*v2+heigvec(3,mh(i))*v3+heigvec(4,mh(i))*v4)
ggprod(i) = h2glgl(hmass(mh(i)),Yteff,Ybeff) ! Gabe's code gives prod ratio
gagaprod(i) = h2gaga(hmass(mh(i)),Yteff,Ybeff,geff) ! Gabe's code gives prod ratio
c          mh_gg_prod(i) = ggprod(i)*?
c          mh_gaga_prod(i) = gagaprod(i)*?
Yteff = 0
Ybeff = 0
geff = 0
enddo


c cross section for h_0 as ratios relative to the SM
do i=1,4
if (hmass(mh(i)).ge.100d0) then
xsgg24l(i)=ggprod(i)* (BFcpe(i,7)/BFSM(i,7))
xsgg22l(i)=ggprod(i)* (BFcpe(i,8)/BFSM(i,8))
xsgg2ww(i)=ggprod(i)* (BFcpe(i,5)/BFSM(i,5))
xsgg2gaga(i)=ggprod(i)* (BFcpe(i,6)/BFSM(i,6))

c       if (hmass(mh(i)).gt.80d0.and.hmass(mh(i)).lt.600d0) then
gghprod(i) = ggprod(i)*lineint(hmass(mh(i)),masshgg,hggprod,30)
xsgg2tautau(i)=gghprod(i)*(BFcpe(i,3)) ! the plot limit here is production x BF
c           ATLAS-CONF-2011-132
c          else
c            hggprod(i) = 0d0
c            xsgg2tautau(i) = 0
c          endif
else
xsgg24l(i)=0d0
xsgg22l(i)=0d0
xsgg2ww(i)=0d0
xsgg2gaga(i)=0d0
endif
enddo

c output for testing purposes
c        write(*,*) 'yukawas: Y1 =',Y1,', Y2 =',Y2,' Y2p =',Y2p
c        write(*,*) 'MSSM parameters: mu =',mu_mssm,', tan(beta) =',dtan(beta_mssm)
c        write(*,*) 'MSSM yukawas: Y1 =',Y1_mssm,', Y2 =',Y2_mssm
c        write(*,*) 'Coupling ratio = ',ratio
c      endif

return
end

c-----------------------------------------------------------c
subroutine atlas_check_preJuly4()
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
if (hmass(mh(i)).gt.110.305d0.and.hmass(mh(i)).lt.600d0) then
hwwcutoff = lineint(hmass(mh(i)),masshww,hwwlimits,79)
c	    print *,'ww',i,hmass(mh(i)),xsgg2ww(i),hwwcutoff
if (xsgg2ww(i) .ge. hwwcutoff) then
atlaswwcheck = 0
c          else
c	    write(16,*) i,hmass(mh(i)),xsgg2ww(i),hwwcutoff
endif
endif
if (hmass(mh(i)).gt.200.345d0.and.hmass(mh(i)).lt.600d0) then
hzz2lcutoff = lineint(hmass(mh(i)),massh2l2nu,hZZlnulimits,50)
c	    print *,'zz',i,hmass(mh(i)),xsgg22l(i),hzz2lcutoff
if (xsgg22l(i) .ge. hzz2lcutoff) then
atlas2lcheck = 0
endif
endif
if (hmass(mh(i)).gt.110.152d0.and.hmass(mh(i)).lt.598.344d0) then
hzz4lcutoff = lineint(hmass(mh(i)),massh4l,hZZ4llimits,130)
c	    print *,'4l',i,hmass(mh(i)),xsgg24l(i),hzz4lcutoff
if (xsgg24l(i) .ge. hzz4lcutoff) then
atlas4lcheck = 0
endif
endif
if (hmass(mh(i)).gt.110.457d0.and.hmass(mh(i)).lt.150.102d0) then
hgagacutoff = lineint(hmass(mh(i)),masshgaga,hgagalimits,67)
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

c-----------------------------------------------------------c
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

c-----------------------------------------------------------c
subroutine happypoint()
c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c
implicit none
include 'manyhiggs.inc'
logical whscan
common/scantype/whscan

if (BFt(1)*BFch(1,3).gt.0.04d0) then
happy1 = 0
endif
if (BFt(1)*BFch(1,1).gt.0.14d0) then
happy2 = 0
endif

if(whscan) then
if (BFch(1,2)*mch_tb_prod(1).lt.3) happy3 = 0
if (BFch(2,4)*mch_tb_prod(2).lt.3) happy3 = 0
endif
c	if(hmass(mh(1)).gt.114d0) happy3 = 0
c	if (chmass(mch(2)).gt.600d0) happy3 = 0

return
end

c-----------------------------------------------------------c
subroutine bsgam_1only()
c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c
implicit none
include 'manyhiggs.inc'
double precision Au(3),Ad(3),Yteff(3),Ybeff(3),mchtemp(3),chBF1!,chpull

c bsg_nlo subroutine
do i=1,3
mchtemp(i) = chmass(mch(i))
Yteff(i) = (Yt*cheigvec(2,mch(i)) + Ytp*cheigvec(4,mch(i)))
Ybeff(i) = Yb*cheigvec(1,mch(i))
Au(i) = (Yt*cheigvec(2,mch(i)) + Ytp*cheigvec(4,mch(i)))/(dsqrt(2d0)*mt/246d0)
Ad(i) = Yb*cheigvec(1,mch(i))/(dsqrt(2d0)*mb/246d0)
c        write(*,*) mchtemp(i),Yteff(i),Ybeff(i)
enddo

call bsg_nlo(Au(1),Ad(1),mchtemp(1),1,chBF1,chpull1)
c        call bsg_nlo(Au,Ad,mchtemp,3,chBF,chpull)
c make sure the effective Y's are scaled by the SM coupling, which they now are. May 23, 2012

if(chpull1.ge.(1.96d0)) then
bsgcheck1 = 0
endif

c        if (higgscheck*chiggscheck*cpoddcheck*charginocheck*neutralinocheck.gt.0) then
c        write(*,*) chpull, bsgcheck
c        endif

do i=1,3
mchtemp(i) = 0d0
Yteff(i) = 0d0
Ybeff(i) = 0d0
Au(i) = 0d0
Ad(i) = 0d0
enddo


return
end

c-----------------------------------------------------------c
subroutine bsgam()
c-----------------------------------------------------------c
c                                                           c
c  This subroutine checks that the charged Higgs branching  c
c  hasn't already been excluded.                            c
c                                                           c
c-----------------------------------------------------------c
implicit none
include 'manyhiggs.inc'
double precision Au(3),Ad(3),Yteff(3),Ybeff(3),mchtemp(3),chBF3!,chpull

c bsg_nlo subroutine
do i=1,3
mchtemp(i) = chmass(mch(i))
Yteff(i) = (Yt*cheigvec(2,mch(i)) + Ytp*cheigvec(4,mch(i)))
Ybeff(i) = Yb*cheigvec(1,mch(i))
Au(i) = (Yt*cheigvec(2,mch(i)) + Ytp*cheigvec(4,mch(i)))/(dsqrt(2d0)*mt/246d0)
Ad(i) = Yb*cheigvec(1,mch(i))/(dsqrt(2d0)*mb/246d0)
c        write(*,*) mchtemp(i),Yteff(i),Ybeff(i)
enddo

c        call bsg_nlo(Au(1),Ad(1),mchtemp(1),1,chBF,chpull)
call bsg_nlo(Au,Ad,mchtemp,3,chBF3,chpull3)
c make sure the effective Y's are scaled by the SM coupling, which they now are. May 23, 2012

if(chpull3.ge.(1.96d0)) then
bsgcheck3 = 0
endif

c        if (higgscheck*chiggscheck*cpoddcheck*charginocheck*neutralinocheck.gt.0) then
c        write(*,*) chpull, bsgcheck
c        endif

do i=1,3
mchtemp(i) = 0d0
Yteff(i) = 0d0
Ybeff(i) = 0d0
Au(i) = 0d0
Ad(i) = 0d0
enddo


return
end
