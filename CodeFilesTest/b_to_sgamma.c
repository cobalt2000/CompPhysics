//
//  b_to_sgamma.c
//  
//
//  Created by Valerie Plaus on 7/15/16.
//
//

#include "b_to_sgamma.h"

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double bsg_nlo(double xAu,double xAd,double xMH,double xnhiggs,double BF,double pull){
/*c
c	Calculates the B -> X_S + gamma BF from the Standard model and 2HDM
c
c	M. Ciuchini, G. Degrassi, P. Gambino, and G. F. Giudice, Nucl. Phys. B527, 21 (1998).
c
c
c	xAu - Charged Higgs PL coupling (proportional to mu)
c	xAd - Charged Higgs PR coupling (proportional to md)
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

//real*8 xAu(maxhiggs),xAd(maxhiggs),xMh(maxhiggs)
    int i;
    double z;

    double BFthref,BFthrefunc;
    double BFexp,BFexpunc;
    double relunc;
    double pull; //The number of standard deviations from


    nhiggs = xnhiggs;
    do i=1,nhiggs{ // I don't even know why he defined the two separate values... after this, they're the same.
        Au(i) = xAu(i);
        Ad(i) = xAd(i);
        MHbsg(i) = xMH(i);
//c        write(*,*) xMH(i),xAu(i),xAd(i)
//c        write(*,*) MHbsg(i),Au(i),Ad(i)
    }

order = 'NLO';

    pibsg = 4*atan(1);

    mcbsg = 1.41;
    mbbsg = 4.75;
    mbartpole = 175;
    mtbsg = 175;
    nf = 5;
    MZbsg = 91.18;
    MWbsg = 80.39;


//c	For comparison w/ the Literature
//c	mb = 4.8d0
//c	mc = mb - 3.39d0
//c	MW = 80.33d0

    mubarb = mbbsg;
    mub = mbbsg;
//c	mub = 5d0

    beta0 = 11-2/3*nf;
    beta1 = 102-38/3*nf;

    gam0m = 8;
    gam1m = 404/3-40/9*nf;


    if(order.eq.'LO') {
        beta1 = 0;
        gam1m = 0;
    }



    eta = alphasbsg(mwbsg)/alphasbsg(mub);

    alphae = 1d0/130.3;
    CKMproductSq = 0.95;

    z = mcbsg**2/mbbsg**2;

    lam1 = 0;	//! drops out anyway

    lam2 = 0.12; //! GeV**2

    delNPc = -1/9*lam2/C0b(7)*(C0b(2) - C0b(1)/6);

    delNPSL = lam1/2 + 3d0*lam2/2*(1-4*(1-z)**4/f(z));
    delNPgam = lam1/2 - 9d0/2*lam2;
    del = 0.99;

//c	write(*,*)(C0beff(i),i=1,8)
//c	write(*,*)C1beff(7)
//c	stop


    BFcenu = 0.1049;	//! from PDG

    BF = BFcenu * CKMproductSq * 6d0*alphae/(pibsg*f(z)*kap(z))*mbrun(mub)**2/mbbsg**2*(DtermSq(mub) + Aterm(mub))*(1d0- delNPSL / mbbsg**2 + delNPgam/mbbsg**2 + delNPc / mcbsg**2);

//c	Comparison with Experiment
//c	Taken from HFAG world average: 1010.1589
    BFexp = 3.55*0.0001; //d-4 in Fortran
    BFexpunc = 0.256*0.0001; //d-4 in Fortran

//c	Taken as the SM reference
    BFthref = 3.62*0.0001; //d-4 in Fortran
    BFthrefunc = 0.33*0.0001; //d-4 in Fortran

    relunc = BFthrefunc/BFthref;


    pull = dabs((BF-BFexp)/dsqrt(BFexpunc**2 + (BF * relunc)**2));




    return;
    }



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function DtermSq(Q)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 Q,sumreal,sumimag
real*8 r1r,r1c,r2r,r2c,r7,r8r,r8c
real*8 logz
real*8 z
real*8 zeta3
real*8 dtermreal,dtermimag

zeta3 = 1.20206d0

z = mcbsg**2/mbbsg**2


logz = dlog(z)


r2r = 2d0/243d0*(-833d0+144d0*pibsg**2*z**1.5d0
                 .	+ (1728d0-180d0*pibsg**2-1296d0*zeta3+(1296d0-324d0*pibsg**2)*logz + 108d0*Log(z)**2+36d0*logz**3)*z
                 .	+ (648d0+72d0*pibsg**2+(432d0-216d0*pibsg**2)*logz+36d0*logz**3)*z**2
                 .	+ (-54d0-84d0*pibsg**2+1092d0*logz-756d0*logz**2)*z**3)
r2c = 16d0*pibsg/81d0*(-5d0+(45d0-3d0*pibsg**2+9d0*logz+9d0*logz**2)*z
                       .	+ (-3d0*pibsg**2+9d0*logz**2)*z**2 + (28d0-12d0*logz)*z**3)

r7 = -10d0/3d0-8d0/9d0*pibsg**2

r8r = -4d0/27d0*(-33d0+2d0*pibsg**2)
r8c = 24d0/27d0*pibsg

r1r = -1d0/6d0*r2r
r1c = -1d0/6d0*r2c




sumreal = 0d0
sumreal = sumreal + C0beff(1)*(r1r + gam0eff(1,7)*dlog(mbbsg/mub))
sumreal = sumreal + C0beff(2)*(r2r + gam0eff(2,7)*dlog(mbbsg/mub))
sumreal = sumreal + C0beff(7)*(r7 + gam0eff(7,7)*dlog(mbbsg/mub))
sumreal = sumreal + C0beff(8)*(r8r + gam0eff(8,7)*dlog(mbbsg/mub))


sumimag = 0d0
sumimag = sumimag + C0beff(1)*r1c
sumimag = sumimag + C0beff(2)*r2c
sumimag = sumimag + C0beff(8)*r8c


dtermreal = C0beff(7) + alphasbsg(mub)/(4d0*pibsg)*(C1beff(7) + sumreal)
dtermimag = alphasbsg(mub)/(4d0*pibsg)*(sumimag)

DtermSq = (dtermreal**2 + dtermimag**2)

//c	write(*,*)C0beff(7),alphas(mub)/(4d0*pi)*C1beff(7),sumreal,sumimag

return
end



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function gam0eff(i,j)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
integer i,j

if(i.eq.1) gam0eff = -208d0/243d0
if(i.eq.2) gam0eff = 416d0/81d0
if(i.eq.3) gam0eff = -176d0/81d0
if(i.eq.4) gam0eff = -152d0/243d0
if(i.eq.5) gam0eff = -6272d0/81d0
if(i.eq.6) gam0eff = 4624d0/243d0
if(i.eq.7) gam0eff = 32d0/3d0
if(i.eq.8) gam0eff = -32d0/9d0


gam0eff = gam0eff*alphasbsg(mub)/(4d0*pibsg)

return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function Aterm(Q)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 Q,sum
integer i,j

sum = 0d0
do i=1,8
do j=i,8
sum = sum + C0beff(i)*C0beff(j)*fij(i,j)
enddo
//c	 write(*,*)i,C0beff(i)
enddo

//c	write(*,*)sum

Aterm = (dexp(-alphasbsg(mub)*dlog(del)*(7d0+2d0*dlog(del))/(3d0*pibsg))-1d0)*
.	C0beff(7)**2 +
.	alphasbsg(mub)/pibsg * sum


return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function mbart(Q)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 Q


mbart = mbartpole * (alphasbsg(Q)/alphasbsg(mtbsg))**(gam0m/(2d0*beta0))*
.	(1d0 + alphasbsg(mtbsg)/(4d0*pibsg)*gam0m/(2d0*beta0)*(gam1m/gam0m-beta1/beta0)*
     .		(alphasbsg(Q)/alphasbsg(mtbsg)-1d0))

//c	write(*,*)mbart,mbartpole,alphas(Q),gam0m,beta0,gam1m,beta1,nf
return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function fij(i,j)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
integer i,j
include 'bsg_nlo.inc'

if(del.ne.0.99d0) then
write(*,*)'Delta not 0.99'
stop
endif

if(i.eq.1.and.j.eq.1) fij =  0.0009d0
if(i.eq.1.and.j.eq.2) fij = -0.0113d0
if(i.eq.1.and.j.eq.3) fij = -0.0035d0
if(i.eq.1.and.j.eq.4) fij =  0.0006d0
if(i.eq.1.and.j.eq.5) fij = -0.0459d0
if(i.eq.1.and.j.eq.6) fij = -0.0600d0
if(i.eq.1.and.j.eq.7) fij = -0.0030d0
if(i.eq.1.and.j.eq.8) fij =  0.0010d0

if(i.eq.2.and.j.eq.2) fij =  0.0340d0
if(i.eq.2.and.j.eq.3) fij =  0.0210d0
if(i.eq.2.and.j.eq.4) fij = -0.0035d0
if(i.eq.2.and.j.eq.5) fij =  0.2754d0
if(i.eq.2.and.j.eq.6) fij =  0.3599d0
if(i.eq.2.and.j.eq.7) fij =  0.0182d0
if(i.eq.2.and.j.eq.8) fij = -0.0061d0

if(i.eq.3.and.j.eq.3) fij =  0.0140d0
if(i.eq.3.and.j.eq.4) fij = -0.0047d0
if(i.eq.3.and.j.eq.5) fij =  0.3277d0
if(i.eq.3.and.j.eq.6) fij =  0.0666d0
if(i.eq.3.and.j.eq.7) fij =  0.0421d0
if(i.eq.3.and.j.eq.8) fij = -0.0140d0

if(i.eq.4.and.j.eq.4) fij =  0.0088d0
if(i.eq.4.and.j.eq.5) fij = -0.0546d0
if(i.eq.4.and.j.eq.6) fij =  0.1570d0
if(i.eq.4.and.j.eq.7) fij = -0.0070d0
if(i.eq.4.and.j.eq.8) fij =  0.0023d0

if(i.eq.5.and.j.eq.5) fij =  1.9369d0
if(i.eq.5.and.j.eq.6) fij =  0.9506d0
if(i.eq.5.and.j.eq.7) fij =  0.5926d0
if(i.eq.5.and.j.eq.8) fij = -0.1975d0


if(i.eq.6.and.j.eq.6) fij =  0.9305d0
if(i.eq.6.and.j.eq.7) fij =  0.0002d0
if(i.eq.6.and.j.eq.8) fij = -0.0001d0

if(i.eq.7.and.j.eq.7) fij =  3.4211d0
if(i.eq.7.and.j.eq.8) fij =  0.3897d0

if(i.eq.8.and.j.eq.8) fij =  3.2162d0



return
end



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function mbrun(Q)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 Q

mbrun = Q*(1d0-4d0/3d0*alphasbsg(Q)/pibsg)

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function alphasbsg(Q)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 Q,v,alphasMZ


alphasMZ = 0.118d0

v = 1d0+beta0*alphasMZ/(2d0*pibsg)*dlog(Q/MZbsg)
//c	write(*,*)v,beta0,beta1,alphasMZ,pi,q,mz
//c	stop

alphasbsg = alphasMZ / v *(1d0-beta1/beta0*alphasMZ/(4d0*pibsg)*dlog(v)/v)

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function kap(z)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 z

kap = 1d0 - 2d0*alphasbsg(mubarb)/(3d0*pibsg)*h(z)/f(z)

if(order.eq.'LO') kap = 1d0

return
end



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function f(z)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 z


f = 1d0-8d0*z+8d0*z**3-z**4-12d0*z**2*dlog(z)
return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function h(z)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 z
real*8 term1,term2,term3,term4,term5,term6,term7


term1 = -(1d0-z**2)*(25d0/4d0-239d0/3d0*z+25d0/4d0*z**2)
term2 = z*dlog(z)*(20d0+90d0*z-4d0/3d0*z**2+17d0/3d0*z**3)
term3 = z**2*dlog(z)**2*(36d0+z**2)
term4 = (1d0-z**2)*(17d0/3d0-64d0/3d0*z+17d0/3d0*z**2)*dlog(1d0-z)
term5 = -4d0*(1d0+30d0*z**2+z**4)*dlog(z)*dlog(1d0-z)
term6 = -(1d0+16d0*z**2+z**4)*(6d0*Li2(z) - pibsg**2)
term7 = -32d0*z**1.5d0*(1d0+z)*(pibsg**2-4d0*Li2(dsqrt(z))+4d0*Li2(-dsqrt(z))-2d0*dlog(z)*dlog((1d0-dsqrt(z))/(1d0+dsqrt(z))))



h = term1 + term2 + term3 + term4 + term5 + term6 + term7

return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function Li2(z)			!either positive or negative arguments allowerd
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 z

if(z.lt.0d0) then
Li2 = 0.5d0*Li2p(z**2) - Li2p(dabs(z))
return
else
Li2 = Li2p(z)
endif



return
end


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function Li2p(z)			! positive agruments only
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
real*8 z
real*8 sum
integer i

sum = 0d0
z = dabs(z)

if(z.lt.1d0) then
do i=1,50
sum = sum + z**i/(1d0*i)**2
enddo
elseif(z.ge.1d0) then
sum = pibsg**2/3d0-1d0/2d0*dlog(z)**2
do i=1,50
sum = sum - 1d0/(i**2*z**i)
enddo
endif
Li2p = sum

return
end




//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function C1beff(i)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
integer i,j
real*8 sum
real*8 x

x = mbart(mwbsg)**2/mwbsg**2

if(i.eq.7) then
sum = 0d0
do j=1,8
sum = sum + (ei(j)*eta*E(x)+fi(j)+gi(j)*eta)*eta**(ai(j))
enddo

C1beff = eta**(39d0/23d0)*C1Weff(7)+8d0/3d0*(eta**(37d0/23d0)-eta**(39d0/23d0))*C1Weff(8) +
.	(297664d0/14283d0*eta**(16d0/23d0) - 7164416d0/357075d0*eta**(14d0/23d0)
     .	 +256868d0/14283d0*eta**(37d0/23d0) - 6698884d0/357075d0 * eta**(39d0/23d0))*
.	C0W(8) + 37208d0/4761d0*(eta**(39d0/23d0)-eta**(16d0/23d0))*C0W(7) + sum


else
stop 'wtf c1beff'

endif


//c	write(*,*)C1beff
//c	write(*,*)eta**(39d0/23d0)*C1Weff(7)+8d0/3d0*(eta**(37d0/23d0)-eta**(39d0/23d0))*C1Weff(8)
//c	write(*,*)C1beff - (eta**(39d0/23d0)*C1Weff(7)+8d0/3d0*(eta**(37d0/23d0)-eta**(39d0/23d0))*C1Weff(8))
//c	stop

if(order.eq.'LO') C1beff = 0d0

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function C0beff(i)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
integer i,j
real*8 sum




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

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0b(int i) {
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
            .	+ 0.0237*pow(eta,0.4086) - 0.0173*pow(eta,(-0.423))
            .	- 0.1336*pow(eta,(-0.8994)) - 0.0316*pow(eta,(0.1456));
            break;
        case 5: thing = 1/108*pow(eta,(-12/23)) - 1/126*pow(eta,(6/23))
            .	+ 0.0094*pow(eta,0.4086) - 0.0100*pow(eta,(-0.423))
            .	+ 0.0010*pow(eta,(-0.8994)) - 0.0017*pow(eta,(0.1456));
            break;
        case 6: thing = -1/36*pow(eta,(-12/23)) - 1/84*pow(eta,(6/23))
            .	+ 0.0108*pow(eta,0.4086) + 0.0163*pow(eta,(-0.423))
            .	+ 0.0103*pow(eta,(-0.8994)) + 0.0023*pow(eta,(0.1456));
            break;
    }
    
    return thing;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C1Weff(int i) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double thing;
    double q,x;

    q = mwbsg;
    x = mbart(q)**2/mwbsg**2;

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
    switch(i){
        case 1: thing = 15+6*log(q**2/mwbsg**2);
            break;
        case 4: thing = E(x) - 2/3 + 2/3*log(q*q/(mwbsg*mwbsg)) + C1WHeff(i);
            break;
        case 7: thing = G7(x) + Delta7(x)*log(q*q/(mwbsg*mwbsg)) + C1WHeff(i);
            break;
        case 8: thing = G8(x) + Delta8(x)*log(q*q/(mwbsg*mwbsg)) + C1WHeff(i);
            break;
        default: break;
    }

    return thing;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double E(double x){
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(-18+11*x+x*x)/(12*(x-1)*(x-1)*(x-1))
            + x**2*(15-16*x+4*x*x)/(6*(x-1)*x-1)*x-1)*x-1))*log(x)
            - 2/3*log(x);


    return what;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G7(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;
    double term1,term2,term3,term4;


    term1 = (-16*x*x*x*x-122*x*x*x+80*x*x-8*x)/(9*pow((x-1),4))*Li2(1-1/x);
    term2 = (6*x*x*x*x+46*x*x*x-28*x*x)/(3*(x-1)**5)*pow(log(x),2);
    term3 = (-102*pow(x,5)-588*x*x*x*x-2262*x*x*x+3244*x*x-1364*x+208)/(81*pow((x-1),5))*log(x);
    term4 = (1646*x**4+12205*x**3-10740*x**2+2509*x-436)/(486*(x-1)**4);

    what = term1 + term2 + term3 + term4;


return
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G8(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;
    double term1,term2,term3,term4;


    term1 = (-4*x*x*x*x+40*x*x*x+41*x*x+x)/(6*(x-1)*(x-1)*(x-1)*(x-1))*Li2(1-1/x);
    term2 = (-17*x*x*x-31*x*x)/(2*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*pow(log(x),2);
    term3 = (-210*x*x*x*x*x+1086*x*x*x*x+4893*x*x*x+2857*x*x-1994*x+280)/(216*(x-1)**5)*log(x);
    term4 = (737*x*x*x*x-14102*x*x*x-28209*x*x+610*x-508)/(1296*(x-1)*(x-1)*(x-1)*(x-1));

    what = term1 + term2 + term3 + term4;

return
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta7(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;

    whatever = (208-1111*x+1086*x*x+383*x*x*x+82*x*x*x*x)/(81*(x-1)*(x-1)*(x-1)*(x-1))
            + (2*x*x*(14-23*x-3*x*x))/(3*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double Delta8(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;


    whatever = (140-902*x-1509*x*x-398*x*x*x+77*x*x*x*x)/(108*(x-1)*(x-1)*(x-1)*(x-1))
        + (x*x*(31+17*x))/(2*(x-1)*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function C0Weff(i)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
integer i

C0Weff = C0W(i)

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function C0W(i)
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
include 'bsg_nlo.inc'
integer i
real*8 Q,x

q = mwbsg

x = mbart(Q)**2/MWbsg**2

//c	write(*,*)x,mbart(Q),mw,q
//c	stop
if(i.eq.2) then
C0W = 1d0
return
endif

if(i.le.6) then
C0W = 0d0
return
endif

if(i.eq.7) then
C0W = F17(x)
endif

if(i.eq.8) then
C0W = F18(x)
endif

C0W = C0W + C0WHiggs(i)

return
end

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F17(double x){
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(7-5*x-8*x*x)/(24*(x-1)*(x-1)*(x-1))
            + x*x*(3*x-2)/(4*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F18(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double what;

    what = x*(2+5*x-x*x)/(8*(x-1)*(x-1)*(x-1))
            - 3*x*x/(4*(x-1)*(x-1)*(x-1)*(x-1))*log(x);


    return what;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C0WHiggs(int i) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    
    int j;
    double sum;
    double x;


    sum = 0;
    for (j=1,nhiggs){
        x = mbart(mwbsg)**2/mhbsg(j)**2;
        if(i.eq.7) sum = sum + Au(j)**2/3d0*F17(x) - Au(j)*Ad(j)*F27(x);
        if(i.eq.8) sum = sum + Au(j)**2/3d0*F18(x) - Au(j)*Ad(j)*F28(x);
    }

//C0WHiggs = sum
    return sum;
    }

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F27(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;

    whatever = x*(3-5*x)/(12*(x-1)*(x-1)) + x*(3*x-2)/(6*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double F28(double x) {
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
    whatever = x*(3-x)/(4*(x-1)*(x-1)) - x/(2*(x-1)*(x-1)*(x-1))*log(x);


    return whatever;
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double C1WHeff(int i){
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    int j;
    double x;
    double sum =0;

//    whatever = 0;
    return; //Why is this here?  It returns 0, then a double, if there's a double? Why not wait?
//    sum = 0;
    for (j=1;nhiggs;j++){ //I believe that this is where the function starts adding up contributions from multiple Higgs, hence the index "nhiggs".
        x = mbart(mwbsg)**2/mhbsg(j)**2;

        if(i.eq.4) sum = sum + EH(x,Au(j),Ad(j));
        if(i.eq.7) sum = sum + G7H(x,Au(j),Ad(j))
            + Delta7H(x,Au(j),Ad(j))*log(mwbsg**2/MHbsg(j)**2)
            - 4/9*EH(x,Au(j),Ad(j));
        if(i.eq.8) sum = sum + G8H(x,Au(j),Ad(j))
            + Delta8H(x,Au(j),Ad(j))*log(mwbsg**2/MHbsg(j)**2)
            - 1/6*EH(x,Au(j),Ad(j));
    }

//    whatever = sum;

    return sum;
    }


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double EH(const double x,const double xAu,const double xAd){
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
//real*8 xAu,xAd,x

    whatever = xAu*xAu*(x*(16-29*x+7*x*x)/(36*pow((x-1),3))+x*(3*x-2)/(6*pow((x-1),4))*log(x));


    return whatever;
}



//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double G7H(const double x,const double xAu,const double xAd){
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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    double whatever;
    double xAu,xAd,x;

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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

//    double xAu,xAd,x;
    double what;


    what = xAu*xAd*1/3*x*( (81-16*x+7*x*x)/(2*pow((x-1),3))-(19+17*x)/pow((x-1),4)*log(x) ) + xAu*xAu*1/6*x*( (-38-261*x+18*x*x-7*x*x*x)/(6*pow((x-1),4))+x*(31+17*x)/pow((x-1),5)*log(x) );

    return what;
}

