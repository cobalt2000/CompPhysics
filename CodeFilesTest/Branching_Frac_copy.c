//
//  Branching_Frac.c
//  
//
//  Created by Valerie Plaus on 7/14/16.
//
//  Updated by Valerie Plaus on 3/28/17.
//  Updated by Valerie Plaus on 6/16/17.

#include "Branching_Frac.h"
//c-----------------------------------------------------------c
branchingfractions(cpehiggs,cpohiggs,chhiggs){
/*c-----------------------------------------------------------c
c                                                           c
c  This subroutine calculates the branching fractions       c
c  of the charged higgses.                                  c
c                                                           c
c-----------------------------------------------------------c*/

    //instead of initializing all the elements to zero, use "memset"
    // i.e. memset(A, 0, sizeof(double).n)
    //and memcopy(*A, *B, sizeof(A)), IF sizeof(A)=sizeof(double).n
    // where B is copied into A (B is the source, A is the new location)

// New challenge: number of particles is unknown (fixed for now), and number of decay chains is unknown
// Each particle struct has an entry for storing the number of decay chains, and a pointer to an array with their widths
// In the current set up, there are
//    10? decay chains for each of the 3 charged higgs particles
//      (3 for the CP even higgses, 3 for the CP odd higgses, 2 for the lighter charged higges, tb, and cs)
//    8 decay chains for the 4 CP even higgs, (bb, cc, tau tau, mu mu, WW, ga ga, ZZ->2l 2nu, ZZ->4l)
//    and 12 decay chains for the 3 CP odd higgs (bb, cc, tau tau, mu mu, WW, ga ga, ZZ->2l 2nu, ZZ->4l,2 heavier CP odd to 2 lighter ones and any CP even (2x4 options for the heaviest CP odd, 1x4 for the next heaviest))

    
//c parameters used in this subroutine only
//    double precision Yteff,Ybeff,geff;
//    double decaych4(6), decaych5(4), decaych6(2), interestingcut, placeholder, decaych7(2), decaych8(2);
//    int ch1check, ch2check, i1, i2, i3;
//    double precision yuktemp, coupling;
//    double complex yukAtemp;

    double x;
    int cpen; // for finding the length of the decay arrays for the cp even higgs
    int cpon;
    int chn;
    int i,j,k;
    cpen=9;
    cpon=9+cpohiggs[0].evec_size;
    chn=4+cpehiggs[0].evec_size+ cpohiggs[0].evec_size+(cpehiggs[0].evec_size+cpohiggs[0].evec_size)*chhiggs[0].evecsize;
    for (i=0,i<cpehiggs[0].evec_size;i++){
        memset(cpehiggs[i].decay,0,sizeof(double)*cpen);
        memset(cpehiggs[i].branching_frac,0,sizeof(double)*cpen+1);
        memset(cpohiggs[i].decay,0,sizeof(double)*cpon);
        memset(cpohiggs[i].branching_frac,0,sizeof(double)*cpon+1);
    }
    for (i=0,i<chhiggs[0].evec_size;i++){
        memset(chhiggs[i].decay,0,sizeof(double)*chn);
        memset(chhiggs[i].branching_frac,0,sizeof(double)*(chn+1));
    }
    
//c Reinitializing the decay widths.
//c charged higgs widths
// these branching fractions should be initilized inside the structs
/*    do i=1,3 {
        totaldecaych[i] = 0.0; //chhiggs[i].decay[0]
        decaychtoh1[i] = 0.0; //chhiggs[i].decay
        decaychtoa1[i] = 0.0; //chhiggs[i].decay
        decaychtoh2[i] = 0.0; //chhiggs[i].decay
        decaychtoa2[i] = 0.0; //chhiggs[i].decay
        decaychtoh3[i] = 0.0; //chhiggs[i].decay
        do j=1,8 {
            BFch[i][j] = 0.0; //chhiggs[i].branching_frac
            decaych[i][j] = 0.0; //chhiggs[i].decay
        }
        BFchtoh1[i] = 0.0; //chhiggs[i].branching_frac
        BFchtoa1[i] = 0.0; //chhiggs[i].branching_frac
        BFchtoh2[i] = 0.0; //chhiggs[i].branching_frac
        BFchtoa2[i] = 0.0; //chhiggs[i].branching_frac
        BFchtoh3[i] = 0.0; //chhiggs[i].branching_frac
    }


//c cpe decay widths, ours and SM.  we only check the lightest/kinematically interesting ones
    do i=1,4 {
        totaldecaycpe[i] = 0.0; //cpehiggs[i].decay[0]
        totaldecaySM[i] = 0.0; //?
        BFinv_decay[i] = 0.0; //cpehiggs[i].decay[9]
        inv_decay[i] = 0.0; //cpehiggs[i].decay[9]
        do j=1,7 {
            decaycpe[i][j] = 0.0; //cpehiggs[i].decay[j]
            decaySM[i][j] = 0.0; //cpehiggs[i].decay[j]
        }
        do j=1,8 {
            BFcpe[i][j] = 0.0; //cpehiggs[i].branching_frac[j]
            BFSM[i][j] = 0.0; //cpehiggs[i].branching_frac[j]
        }
    }

//c cpo decay widths.
    do i=1,3 {
        totaldecaycpo[i] = 0.0; //cpohiggs[i].decay[0]
        do j=1,12 {
            BFcpo[i][j] = 0.0; //cpohiggs[i].branching_frac[j]
            decaycpo[i][j] = 0.0; //cpohiggs[i].decay[j]
        }
    }


//c      write(*,*) 'Running the subroutine'
//c calculating the branching fractions of the massive charged Higgs
    do i=1,3 {

//c resetting some mumbo jumbo
        do j=1,2 {
            decaych8[j] = 0.0;
            decaych7[j] = 0.0;
            decaych6[j] = 0.0;
        }
        do j=1,4 {
            decaych5[j] = 0.0;
        }
        do j=1,6 {
            decaych4[j] = 0.0;
        }
*/
//c decay to cs
//c      write(*,*) i
/*        decaych[i][1] = 3.0*pf(chmass[mch[i]],1.42,0.104)*
        (
         ( (Ys*cheigvec[1][mch[i]])*(Ys*cheigvec[1][mch[i]]) + (Yc*cheigvec[2][mch[i]] + Ycp*cheigvec[4][mch[i]] )*(Yc*cheigvec[2][mch[i]] + Ycp*cheigvec[4][mch[i]] ))*
             (chmass[mch[i]]*chmass[mch[i]]-1.42*1.42 - 0.104*0.104) +
             4. *1.42 *0.104 *(Ys*cheigvec[1][mch[i]])*(Yc*cheigvec[2][mch[i]] + Ycp*cheigvec[4][mch[i]])
        )
        /(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
 */
        chhiggs[i].decay[1]=3.0*pf(chhiggs[i].mass,1.42,0.104)*
            ((chhiggs[i].Y_s*chhiggs[i].Y_s + chhiggs[i].Y_c*chhiggs[i].Y_c)* (chhiggs[i].mass*chhiggs[i].mass-1.42*1.42 - 0.104*0.104)
            +4. *1.42 *0.104 *chhiggs[i].Y_s*chhiggs[i].Y_c)
            /(8. *pi*chhiggs[i].mass*chhiggs[i].mass);

        
        
//c decay to tb
/*        if (chmass[mch[i]].gt.(174.3 +4.7 )) {
            decaych[i][2] = 3. *pf(chmass[mch[i]],174.3 ,4.7 )*
            (
            (chmass[mch[i]]*chmass[mch[i]] - 174.3 *174.3 - 4.7 *4.7) *
            ((Yb*cheigvec[1][mch[i]])*(Yb*cheigvec[1][mch[i]]) + (Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]])*(Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]]))+
            4. *174.3 *4.7 *(Yb*cheigvec[1][mch[i]])*(Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]])
            )
            /(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }
 */
    if (chhigss[i].mass>(174.3+4.7)) {
        chhiggs[i].decay[2]=3.0*pf(chhiggs[i].mass,174.3,4.7)*
    ((chhiggs[i].Y_b*chhiggs[i].Y_b + chhiggs[i].Y_t*chhiggs[i].Y_t)* (chhiggs[i].mass*chhiggs[i].mass-174.3*174.3 - 4.7*4.7)
    +4. *1.42 *0.104 *chhiggs[i].Y_s*chhiggs[i].Y_c)
    /(8. *pi*chhiggs[i].mass*chhiggs[i].mass);
    }
    

//c decay to tau nu
/*        decaych[i][3] = pf(chmass[mch[i]],1.78 ,0 ) *
        (chmass[mch[i]]*chmass[mch[i]] - 1.78 *1.78) *
        (Ytau * cheigvec[1][mch[i]])*(Ytau * cheigvec[1][mch[i]])/
        (8. *pi*chmass[mch[i]]*chmass[mch[i]]);
*/
        chhiggs[i].decay[3]= pf(chhiggs[i].mass,1.78 ,0 ) * (chhiggs[i].mass*chhiggs[i].mass -1.78*1.78) * chhiggs[i].Y_tau*chhiggs[i].Y_tau
        /(8. *pi* chhiggs[i].mass*chhiggs[i].mass);
 

//c decay to h^0 W: there are 4 possible h^0, which becomes decay array elements 4-7
/*        do j=1,4 {
            if (chmass[mch[i]].gt.(80.4  + hmass[mh[j]])) {
                decaych4[j] = pf(chmass[mch[i]],hmass[mh[j]],80.4 ) *
                (heigvec[1][mh[j]]*cheigvec[1][mch[i]] - heigvec[2][mh[j]]*cheigvec[2][mch[i]] +
                 heigvec[3][mh[j]]*cheigvec[3][mch[i]] - heigvec[4][mh[j]]*cheigvec[4][mch[i]])*
                (heigvec[1][mh[j]]*cheigvec[1][mch[i]] - heigvec[2][mh[j]]*cheigvec[2][mch[i]] +
                heigvec[3][mh[j]]*cheigvec[3][mch[i]] - heigvec[4][mh[j]]*cheigvec[4][mch[i]]) *
                StoSV(chmass[mch[i]],hmass[mh[j]],80.4 ) *
                (g2/2. )*(g2/2.)/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
            }
        }
 */      for (j=0;j<=3;j++) {
            if (chhiggs[i].mass>(80.4  + cpehiggs[j].mass)) {
            x = W_H_rotation(cpehiggs[j].evec, chhiggs[i].evec);
            chhiggs[i].decay[4+j] = pf(chhiggs[i].mass,cpehiggs[j].mass,80.4 ) * x * x *
            StoSV(chhiggs[i].mass,cpehiggs[j].mass,80.4 ) *
            (g2/2. )*(g2/2.)/(8. *pi*chhiggs[i].mass*chhiggs[i].mass);
            x=0;
            }
        }
 
        
//c decay to A^0 W: there are 3 possible A^0, which become elements 8-10
//  The count starts at the second particle (j=1) because the first is massless
/*        do j=1,3 {
            if (chmass[mch[i]].gt.(80.4  + cpomass[ma[j]])) {
                decaych5[j] = pf(chmass[mch[i]],cpomass[ma[j]],80.4 ) *
                (cpoeigvec[1][ma[j]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[j]]*cheigvec[2][mch[i]]
                + cpoeigvec[3][ma[j]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[j]]*cheigvec[4][mch[i]])*
                (cpoeigvec[1][ma[j]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[j]]*cheigvec[2][mch[i]]
                 + cpoeigvec[3][ma[j]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[j]]*cheigvec[4][mch[i]])*
                StoSV(chmass[mch[i]],cpomass[ma[j]],80.4 )*
                (g2/2.)*(g2/2.)/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
            }
        }
*/
    
    for (j=0;j<=3;j++) { //should exclude lightest CP odd, because it's the Goldstone, but we'll leave it in for testing purposes
            if (chhiggs[i].mass>(80.4  + cpohiggs[j].mass)) {
            x = dot_prod(cpohiggs[j].evec, chhiggs[i].evec);
            chhiggs[i].decay[4+cpehiggs[0].evec_size+j] = pf(chhiggs[i].mass,cpohiggs[j].mass,80.4 ) * x * x
            StoSV(chhiggs[i].mass,cpohiggs[j].mass,80.4 ) *
            (g2/2. )*(g2/2.)/(8. *pi*chhiggs[i].mass*chhiggs[i].mass);
            x=0;
            }
         }

        
//c decay to h^+ Z: there are 2 possible h^+_j for i>j
/*        do j=1,2 {
            if (chmass[mch[i]].gt.(zmass[2] + chmass[mch[j]])) {
                decaych6[j] = pf(chmass[mch[i]],chmass[mch[j]],zmass[2]) *
                (cheigvec[1][mch[j]]*cheigvec[1][mch[i]] + cheigvec[2][mch[j]]*cheigvec[2][mch[i]] +
                cheigvec[3][mch[j]]*cheigvec[3][mch[i]] + cheigvec[4][mch[j]]*cheigvec[4][mch[i]])*
                (cheigvec[1][mch[j]]*cheigvec[1][mch[i]] + cheigvec[2][mch[j]]*cheigvec[2][mch[i]] +
                 cheigvec[3][mch[j]]*cheigvec[3][mch[i]] + cheigvec[4][mch[j]]*cheigvec[4][mch[i]])*
                StoSV(chmass[mch[i]],chmass[mch[j]],zmass[2])*
                ((g2*g2-g1*g1)/(G*2 ))*((g2*g2-g1*g1)/(G*2 ))/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
            }
        }
*/
    for (j=0;j<=3;j++) { //j should include only the lightest massive charged higgs, but we'll leave it for testing purposes'
            if (chhiggs[i].mass>(zmass[2] + chhiggs[j].mass)) {
                x = dot_prod(chhiggs[j].evec, chhiggs[i].evec);
                chhiggs[i].decay[4+cpehiggs[0].evec_size+cpohiggs[0].evec_size+j] = pf(chhiggs[i].mass,chhiggs[j].mass,zmass[2]) *x*x*
                StoSV(chhiggs[i].mass,chhiggs[j].mass,zmass[2])*
                ((g2*g2-g1*g1)/(G*2 ))*((g2*g2-g1*g1)/(G*2 ))/(8. *pi*chhiggs[i].mass*chhiggs[i].mass);
                x=0;
            }
        }

/* Temporarily omitted because this channel was suppressed by the large masses involved
//c decay to h^+ h0
        do j=1,2 {
            do k=1,4 {
                if (chmass[mch[i]].gt.(chmass[mch[j]] + hmass[mh[k]])) {
                    do i1=1,4 {
                        do i2=1,4 {
                            do i3=1,4 {
                                yuktemp=yuktemp+yuk[i1][i2][i3]*cheigvec[i1][mch[i]]*cheigvec[i2][mch[j]]*heigvec[i3][mh[k]]; //three way Higgs vertex
                            }
                        }
                    }
                    decaych7[j] = pf(chmass[mch[i]],chmass[mch[j]],hmass[mh[k]]) *
                        (yuktemp)*(yuktemp) /(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
                }
                yuktemp=0;
            }
        }

 */
/* Temporarily omitted: previous scans showed that this was a suppressed channel due to the alignment of A0 and H+
//c decay to h^+ A0
        do j=1,2 {
            do k=1,4 {
                if (chmass[mch[i]].gt.(chmass[mch[j]]+cpomass[ma[k]])) {
                    do i1=1,4 {
                        do i2=1,4 {
                            do i3=1,4 {
                                yukAtemp=yukAtemp+yukA[i1][i2][i3]*cheigvec[i1][mch[i]]*cheigvec[i2][mch[j]]*cpoeigvec[i3][ma[k]]; //three way Higgs vertex
                            }
                        }
                    }
                    decaych8[j] = pf(chmass[mch[i]],chmass[mch[j]],cpomass[ma[k]]) *
                    (yukAtemp)*conjg(yukAtemp)/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
                }
                yukAtemp=0;
            }
        }
*/
/*      Not needed with the current decay array.
        decaych[i][4] = decaych4[1]+decaych4[2]+decaych4[3]+decaych4[4];
        decaych[i][5] = decaych5[1]+decaych5[2]+decaych5[3]+decaych5[4];
        decaych[i][6] = decaych6[1]+decaych6[2];
        decaych[i][7] = decaych7[1]+decaych7[2];
        decaych[i][8] = decaych8[1]+decaych8[2];
        totaldecaych[i] = decaych[i][1] + decaych[i][2] + decaych[i][3] + decaych[i][4] + decaych[i][5] + decaych[i][6]
         + decaych[i][7] + decaych[i][8];
*/
//c      While the charged Higgs may decay through all these modes, it may also be light enough
//c      to be a decay product of the top.  So here is a test for that decay channel.
        if (174.3  .gt. (chmass[mch[i]] + 4.7 )) {
            decayt[i] = pf(174.3 ,chmass[mch[i]],4.7 )*
            (
            (-chmass[mch[i]]*chmass[mch[i]] + 174.3 *174.3 - 4.7 *4.7) *
            ((Yb*cheigvec[1][mch[i]])*(Yb*cheigvec[1][mch[i]]) + (Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]])*(Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]]))+
            4. *chmass[mch[i]]*4.7 *(Yb*cheigvec[1][mch[i]])*(Yt*cheigvec[2][mch[i]] + Ytp*cheigvec[4][mch[i]])
            )
            /(8. *pi*174.3 *174.3)
            BFt[i]=decayt[i] / (decayt[i] + 1.508 );
        }
/* This all appears to be a duplicate of calculations above.
//c decay to lightest h^0 W
        if (chmass[mch[i]].gt.(80.4  + hmass[mh[1]])) {
            decaychtoh1[i]=pf(chmass[mch[i]],hmass[mh[1]],80.4 ) *
           (heigvec[1][mh[1]]*cheigvec[1][mch[i]] - heigvec[2][mh[1]]*cheigvec[2][mch[i]] +
            heigvec[3][mh[1]]*cheigvec[3][mch[i]] - heigvec[4][mh[1]]*cheigvec[4][mch[i]])*
            (heigvec[1][mh[1]]*cheigvec[1][mch[i]] - heigvec[2][mh[1]]*cheigvec[2][mch[i]] +
             heigvec[3][mh[1]]*cheigvec[3][mch[i]] - heigvec[4][mh[1]]*cheigvec[4][mch[i]])*
            StoSV(chmass[mch[i]],hmass[mh[1]],80.4 ) *
            (g2/2. )*(g2/2. )/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }

//c decay to second lightest h^0 W
        if (chmass[mch[i]].gt.(80.4  + hmass[mh[2]])) {
            decaychtoh2[i]=pf(chmass[mch[i]],hmass[mh[2]],80.4 ) *
            (heigvec[1][mh[2]]*cheigvec[1][mch[i]] - heigvec[2][mh[2]]*cheigvec[2][mch[i]] +
            heigvec[3][mh[2]]*cheigvec[3][mch[i]] - heigvec[4][mh[2]]*cheigvec[4][mch[i]])*
            (heigvec[1][mh[2]]*cheigvec[1][mch[i]] - heigvec[2][mh[2]]*cheigvec[2][mch[i]] +
             heigvec[3][mh[2]]*cheigvec[3][mch[i]] - heigvec[4][mh[2]]*cheigvec[4][mch[i]])*
            StoSV(chmass[mch[i]],hmass[mh[2]],80.4 ) *
            (g2/2. )*(g2/2. )/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }

//c decay to third lightest h^0 W
        if (chmass[mch[i]].gt.(80.4  + hmass[mh[3]])) {
            decaychtoh3[i]=pf(chmass[mch[i]],hmass[mh[3]],80.4 ) *
            (heigvec[1][mh[3]]*cheigvec[1][mch[i]] - heigvec[2][mh[3]]*cheigvec[2][mch[i]] +
            heigvec[3][mh[3]]*cheigvec[3][mch[i]] - heigvec[4][mh[3]]*cheigvec[4][mch[i]])*
            (heigvec[1][mh[3]]*cheigvec[1][mch[i]] - heigvec[2][mh[3]]*cheigvec[2][mch[i]] +
             heigvec[3][mh[3]]*cheigvec[3][mch[i]] - heigvec[4][mh[3]]*cheigvec[4][mch[i]])*
            StoSV(chmass[mch[i]],hmass[mh[3]],80.4 ) *
            (g2/2. )*(g2/2. )/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }

//c decay to lightest A^0 W
        if (chmass[mch[i]].gt.(80.4  + cpomass[ma[1]])) {
            decaychtoa1[i]=pf(chmass[mch[i]],cpomass[ma[1]],80.4 )  *
            (cpoeigvec[1][ma[1]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[1]]*cheigvec[2][mch[i]]
            + cpoeigvec[3][ma[1]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[1]]*cheigvec[4][mch[i]])*
            (cpoeigvec[1][ma[1]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[1]]*cheigvec[2][mch[i]]
             + cpoeigvec[3][ma[1]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[1]]*cheigvec[4][mch[i]])*
            StoSV(chmass[mch[i]],cpomass[ma[1]],80.4 )*
            (g2/2. )*(g2/2. )/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }
//c decay to second lightest A^0 W
        if (chmass[mch[i]].gt.(80.4  + cpomass[ma[2]])) {
            decaychtoa2[i]=pf(chmass[mch[i]],cpomass[ma[2]],80.4 )  *
            (cpoeigvec[1][ma[2]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[2]]*cheigvec[2][mch[i]]
            + cpoeigvec[3][ma[2]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[2]]*cheigvec[4][mch[i]])*
            (cpoeigvec[1][ma[2]]*cheigvec[1][mch[i]] + cpoeigvec[2][ma[2]]*cheigvec[2][mch[i]]
             + cpoeigvec[3][ma[2]]*cheigvec[3][mch[i]] + cpoeigvec[4][ma[2]]*cheigvec[4][mch[i]])*
            StoSV(chmass[mch[i]],cpomass[ma[2]],80.4 )*
            (g2/2. )*(g2/2. )/(8. *pi*chmass[mch[i]]*chmass[mch[i]]);
        }


        do j=1,6 {
            BFch[i][j] = decaych[i][j] / totaldecaych[i];
        }
        BFchtoh1[i]=decaychtoh1[i]/totaldecaych[i];
        BFchtoh2[i]=decaychtoh2[i]/totaldecaych[i];
        BFchtoh3[i]=decaychtoh3[i]/totaldecaych[i];
        BFchtoa1[i]=decaychtoa1[i]/totaldecaych[i];
        BFchtoa2[i]=decaychtoa2[i]/totaldecaych[i];
//c       write(*,*) decaychtoh1(i),decaychtoh2(i),decaychtoh3(i)
//c       write(*,*) BFchtoh1(i),BFchtoh2(i),BFchtoh3(i)

        }// ??? must be the end of the i loop
*/
//c CPE Higgs to bb, cc, tau tau, mu mu, WW, ga ga, ZZ->2l 2nu, ZZ->4l
    for (i=0;i<=3,i++) {

            if (cpehiggs[i].mass>=(2 *4.7 )) {
                cpehiggs[i].decay[1] = 3. *cpehiggs[i].mass/(16. *pi) * pow((1. - 4. *4.7 *4.7/(cpehiggs[i].mass*cpehiggs[i].mass)),(1.5) ) *
                (cpehiggs[i].Y_b)*(cpehiggs[i].Y_b);
            }
            cpehiggs[i].decay[2] = 3. *cpehiggs[i].mass/(16. *pi) * pow((1.  - 4. *1.42 *1.42/(cpehiggs[i].mass*cpehiggs[i].mass)),(1.5)) *
            (cpehiggs[i].Y_c)*(cpehiggs[i].Y_c);
            if (cpehiggs[i].mass>=(2 *1.78 )) {
                cpehiggs[i].decay[3] = cpehiggs[i].mass/(16. *pi) * pow((1.  - 4. *1.78 *1.78/(cpehiggs[i].mass*cpehiggs[i].mass)),(1.5 )) *
                (cpehiggs[i].Y_tau)*(cpehiggs[i].Y_tau);
            }
            decaycpe[i][4] = cpehiggs[i].mass/(16. *pi) * pow((1. - 4. *0.106 *0.106/(cpehiggs[i].mass*cpehiggs[i].mass)),(1.5 )) *
            (cpehiggs[i].Y_mu)*(cpehiggs[i].Y_mu);

            if (cpehiggs[i].mass>=100 ) {

//c effective couplings are based on the composition of the mass state with respect to the gauge basis
//                Yteff = (Yt*heigvec[2][mh[i]] + Ytp*heigvec[4][mh[i]]);
//                Ybeff = Yb*heigvec[1][mh[i]];
//                geff = g2*(heigvec[1][mh[i]]*v1+heigvec[2][mh[i]]*v2+heigvec[3][mh[i]]*v3+heigvec[4][mh[i]]*v4);
                geff = g2*dot_prod(cpehiggs[i].evec,vevs);
/*                decaycpe[i][5] = (heigvec[1][mh[i]]*v1+heigvec[2][mh[i]]*v2+heigvec[3][mh[i]]*v3+heigvec[4][mh[i]]*v4)*
                    (heigvec[1][mh[i]]*v1+heigvec[2][mh[i]]*v2+heigvec[3][mh[i]]*v3+heigvec[4][mh[i]]*v4)*
                    lineint(hmass(mh(i)),widthmass,widthww,400)/(246 *246);
 */
                cpehiggs[i].decay[5] = dot_prod(cpehiggs[i].evec,vevs)*dot_prod(cpehiggs[i].evec,vevs)*lineint(cpehiggs[i].mass,widthmass,widthww,400)/(246 *246);
/*                decaycpe[i][6] = h2gaga(hmass[mh[i]],Yteff,Ybeff,geff) *lineint(hmass[mh[i]],widthmass,widthgaga,400);
 */
                cpehiggs[i].decay[6] = h2gaga(cpehiggs[i].mass,cpehiggs[i].Y_t,cpehiggs[i].Y_b,geff) *lineint(cpehiggs[i].mass,widthmass,widthgaga,400);
/*                decaycpe[i][7] = (heigvec[1][mh[i]]*v1+heigvec[2][mh[i]]*v2+heigvec[3][mh[i]]*v3+heigvec[4][mh[i]]*v4)*
                (heigvec[1][mh[i]]*v1+heigvec[2][mh[i]]*v2+heigvec[3][mh[i]]*v3+heigvec[4][mh[i]]*v4)*
                lineint(hmass[mh[i]],widthmass,widthzz,400)/(246 *246);
                */
                cpehiggs[i].decay[7] = dot_prod(cpehiggs[i].evec,vevs)*dot_prod(cpehiggs[i].evec,vevs)*
                    lineint(cpehiggs[i].mass,widthmass,widthzz,400)/(246 *246);

//                Yteff = 0;
//                Ybeff = 0;
                geff = 0;
            }
            else {
                cpehiggs[i].decay[5] = 0;
                cpehiggs[i].decay[6] = 0;
                cpehiggs[i].decay[7] = 0;
            }


//c Here we calculate the partial decay width of each CP even higgs to a pair of the lightest neutralino, which would be invisible to the detectors.
//c Smartin basis: B, W3, H1, H2, H3, H4
//c map old bases to SMartin: 1->3,2->4,3->5,4->6,5->2,6->1
//c old basis: H1, H2, H3, H4, W3, B

//c gauge terms
            coupling = 0. ;

            if (hmass[mh[i]].gt.(2. *chi0mass[mchi0[1]])) {
                coupling = 0. ;
//c        write(*,*) coupling, hmass(mh(i)), chi0mass(mchi0(1)), chi0mass(1)
                coupling = coupling + g1*(-heigvec[1][mh[i]]*chi0vec[3][mchi0[1]] + heigvec[2][mh[i]]*chi0vec[4][mchi0[1]] -
                                          heigvec[3][mh[i]]*chi0vec[5][mchi0[1]] + heigvec[4][mh[i]]*chi0vec[6][mchi0[1]])*chi0vec[1][mchi0[1]];
                coupling = coupling + g2*(heigvec[1][mh[i]]*chi0vec[3][mchi0[1]] - heigvec[2][mh[i]]*chi0vec[4][mchi0[1]] +
                                          heigvec[3][mh[i]]*chi0vec[5][mchi0[1]] - heigvec[4][mh[i]]*chi0vec[6][mchi0[1]])*chi0vec[2][mchi0[1]];

                cpehiggs[i].decay[8] = coupling*coupling*hmass[mh[i]]/(16. *pi)*pow((1. -4. *chi0mass[mchi0[1]]*chi0mass[mchi0[1]]/(hmass[mh[i]]*hmass[mh[i]])),(1.5 ));
            }
            else {
                cpehiggs[i].decay[8] = 0. ;
            }

            cpehiggs[i].decay[0] = cpehiggs[i].decay[1] + cpehiggs[i].decay[2] + cpehiggs[i].decay[3] + cpehiggs[i].decay[4] +
            cpehiggs[i].decay[5] + cpehiggs[i].decay[6] + cpehiggs[i].decay[7] + cpehiggs[i].decay[8];
//c      totaldecaycpe(i)=decaycpe(i,1)+decaycpe(i,2)+decaycpe(i,3)+decaycpe(i,4)+decaycpe(i,5)+decaycpe(i,6)+decaycpe(i,7)
//c     . +inv_decay(i)
        for (j=1;j<=6;j++) {
            cpehiggs[i].branching_frac[j] = cpehiggs[i].decay[j] / cpehiggs[i].decay[0] ;
            }
            cpehiggs[i].branching_frac[7] = cpehiggs[i].decay[7]/cpehiggs[i].decay[0] *(2*(3.3658e-2))*(2*(3.3658e-2)); //!ZZ to 4 leptons (e+mu+tau=3)
            cpehiggs[i].branching_frac[8] = cpehiggs[i].decay[7]/cpehiggs[i].decay[0] *(20e-2)*(2*(3.3658e-2)); //!ZZ to 2 leptons 2 nus
            cpehiggs[i].branching_frac[9] = cpehiggs[i].decay[8]/cpehiggs[i].decay[0];

        }
//c      write(*,*) hmass(mh(1)), decaycpe(1,1),decaycpe(1,2),decaycpe(1,3), decaycpe(1,4)
//c      write(*,*) , decaycpe(1,5) , decaycpe(1,6) , decaycpe(1,7) , inv_decay(1), totaldecaycpe(1)
//c      write(*,*) hmass(mh(1)), BFcpe(1,1),BFcpe(1,2),BFcpe(1,3),BFcpe(1,4),BFcpe(1,1)+BFcpe(1,2)+BFcpe(1,3)+BFcpe(1,4)

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//c decay of the cp odd Higgs to bb,cc,tau tau and mu mu
//c now adding the Ai -> Aj hk mode
    for (i=0;i<=3;i++) { //Goldstone boson is i=0, massive CP odd higgs are i=1,2,3. i=0 should have no contribution!
//            decaycpo[i][1] = 3. *cpomass[ma[i]]/(16. *pi) * pow((1.  - 4. *4.7 *4.7/(cpomass[ma[i]]*cpomass[ma[i]])),(1.5)) *
//            (Yb*cpoeigvec[1][ma[i]])*(Yb*cpoeigvec[1][ma[i]]);
        cpohiggs[i].decay[1] = 3. *cpohiggs[i].mass/(16. *pi) * pow((1. - 4. *4.7 *4.7/(cpohiggs[i].mass*cpohiggs[i].mass)),(1.5) ) *
        (cpohiggs[i].Y_b)*(cpohiggs[i].Y_b);
//        decaycpo[i][2] = 3. *cpomass[ma[i]]/(16. *pi) * pow((1.  - 4. *1.42 *1.42/(cpomass[ma[i]]*cpomass[ma[i]])),(1.5)) *
//            (Yc*cpoeigvec[2][ma[i]] + Ycp*cpoeigvec[4][ma[i]])*(Yc*cpoeigvec[2][ma[i]] + Ycp*cpoeigvec[4][ma[i]]);
       cpohiggs[i].decay[2] = 3. *cpohiggs[i].mass/(16. *pi) * pow((1.  - 4. *1.42 *1.42/(cpohiggs[i].mass*cpohiggs[i].mass)),(1.5)) *
        (cpohiggs[i].Y_c)*(cpohiggs[i].Y_c);
//        decaycpo[i][3] = cpomass[ma[i]]/(16. *pi) * pow((1.  - 4. *1.78 *1.78 /cpomass[ma[i]]*cpomass[ma[i]]),(1.5 )) *
//            (Ytau*cpoeigvec[1][ma[i]])*(Ytau*cpoeigvec[1][ma[i]]);
        cpohiggs[i].decay[3] = cpohiggs[i].mass/(16. *pi) * pow((1.  - 4. *1.78 *1.78/(cpohiggs[i].mass*cpohiggs[i].mass)),(1.5 )) *
        (cpohiggs[i].Y_tau)*(cpohiggs[i].Y_tau);
//        decaycpo[i][4] = cpomass[ma[i]]/(16. *pi) * pow((1.  - 4. *0.106 *0.106/cpomass[ma[i]]*cpomass[ma[i]]),(1.5)) *
//            (Ymu*cpoeigvec[1][ma[i]])*(Ymu*cpoeigvec[1][ma[i]]);
            decaycpo[i][4] = cpohiggs[i].mass/(16. *pi) * pow((1. - 4. *0.106 *0.106/(cpohiggs[i].mass*cpohiggs[i].mass)),(1.5 )) *
        (cpohiggs[i].Y_mu)*(cpohiggs[i].Y_mu);
        
        for (j=1;j<=2;j++) {
            for (k=0;k<=3;k++) {
                if (cpohiggs[i].mass>(cpohiggs[j].mass+cpehiggs[i].mass)){
                    x=W_H_rotation(cpohiggs[i].evec,cpohiggs[j].evec);
                    y=W_H_rotation(cpehiggs[k].evec,vevs);
                    cpohiggs[i].decay[(4*j)+k]=pf(cpohiggs[i].mass,cpohiggs[j].mass,cpehiggs[k])*x*x*y*y*
                    (g2*g2+g1*g1)*(g2*g2+g1*g1)/(4. *8. *pi*cpohiggs[i].mass*cpohiggs[i].mass);
                }
/*                if (cpomass[ma[i]].gt.cpomass[ma[j]]+hmass[mh[k]]) {
                    decaycpo[i][(4*j)+k] = pf(cpomass[ma[i]],cpomass[ma[j]],hmass[mh[k]])*
                    (cpoeigvec[1][ma[i]]*cpoeigvec[1][ma[j]] - cpoeigvec[2][ma[i]]*cpoeigvec[2][ma[j]]
                    + cpoeigvec[3][ma[i]]*cpoeigvec[3][ma[j]] - cpoeigvec[4][ma[i]]*cpoeigvec[4][ma[j]])*
                    (cpoeigvec[1][ma[i]]*cpoeigvec[1][ma[j]] - cpoeigvec[2][ma[i]]*cpoeigvec[2][ma[j]]
                    + cpoeigvec[3][ma[i]]*cpoeigvec[3][ma[j]] - cpoeigvec[4][ma[i]]*cpoeigvec[4][ma[j]])*
                    (heigvec[1][mh[k])*v1 - heigvec[2][mh[k]]*v2 + heigvec[3][mh[k]]*v3 - heigvec[4][mh[k]]*v4)*
                    (heigvec[1][mh[k])*v1 - heigvec[2][mh[k]]*v2 + heigvec[3][mh[k]]*v3 - heigvec[4][mh[k]]*v4)*
                    (g2*g2+g1*g1)**2/(4. *8. *pi*cpomass[ma[i]]*cpomass[ma[i]]);
                    }*/
                }
            }
        for (j=0;j<cpon;j++){
            cpohiggs[i].decay[0]= cpohiggs[i].decay[0]+cpohiggs[i].decay[j];
        }
/*            do j=1,12 {
                totaldecaycpo[i] = totaldecaycpo[i] + decaycpo[i][j];
            }*/
        for (j=0;j<cpon;j++){
            cpohiggs[i].branching_frac[j]= cpohiggs[i].decay[j]/cpohiggs[i].decay[0];
        }
/*        do j=1,12 {
                BFcpo[i][j] = decaycpo[i][j] / totaldecaycpo[i];
            }*/
        }

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//c SM Higgs to bb, cc, tau tau, mu mu, WW, ga ga, ZZ->2l 2nu, ZZ->4l
        do i=1,4 {
            if (hmass[mh[i]].ge.(2 *4.7 )) {
                decaySM[i][1] = 3. *hmass[mh[i]]/(16. *pi) * pow((1.- 4. *4.7*4.7/(hmass[mh[i]]*hmass[mh[i]])),(1.5)) *
                (dsqrt(2 )*4.7 /246 )**2;
            }
            decaySM[i][2] = 3. *hmass[mh[i]]/(16. *pi) * pow((1.- 4. *1.42*1.42/(hmass[mh[i]]*hmass[mh[i]])),(1.5)) *
            (dsqrt(2 )*1.42 /246 )*(dsqrt(2 )*1.42 /246 );
            if (hmass[mh[i]].ge.(2 *1.78 )) {
                decaySM[i][3] = hmass[mh[i]]/(16. *pi) * pow((1.- 4. *1.78*1.78/(hmass[mh[i]]*hmass[mh[i]])),(1.5)) *
                (dsqrt(2 )*1.78 /246 )*(dsqrt(2 )*1.78 /246 );
            }
            decaySM[i][4] = hmass[mh[i]]/(16. *pi) * pow((1.- 4. *0.106*0.106/(hmass[mh[i]]*hmass[mh[i]])),(1.5)) *
            (dsqrt(2 )*0.106 /246 )*(dsqrt(2 )*0.106 /246 );

            if (hmass[mh[i]].ge.100 ) {

                decaySM[i][5] = lineint(hmass[mh[i]],widthmass,widthww,400);
                decaySM[i][6] = lineint(hmass[mh[i]],widthmass,widthgaga,400);
                decaySM[i][7] = lineint(hmass[mh[i]],widthmass,widthzz,400);
            }
            else {
                decaySM[i][5] = 0;
                decaySM[i][6] = 0;
                decaySM[i][7] = 0;
            }

            totaldecaySM[i] = decaySM[i][1] + decaySM[i][2] + decaySM[i][3] + decaySM[i][4] +
            decaySM[i][5] + decaySM[i][6] + decaySM[i][7];

            do j=1,6 {
                BFSM[i][j] = decaySM[i][j] / totaldecaySM[i];
            }
            BFSM[i][7] = decaySM[i][7]/totaldecaySM[i] *(3*(3.3658e-2))*(3*(3.3658e-2)); //!ZZ to 4 leptons
            BFSM[i][8] = decaySM[i][7]/totaldecaySM[i] *(20e-2)*(3*(3.3658e-2)); //!ZZ to 2 leptons 2 nus
//c      write(*,*) hmass(mh(i)), decaySM(i,1),decaySM(i,2),decaySM(i,3), decaySM(i,4)
//c      write(*,*) decaySM(i,5) , decaySM(i,6) , decaySM(i,7)
//c      write(*,*) totaldecaySM(i)
//c      write(*,*) decaySM(i,1)+decaySM(i,2)+decaySM(i,3)+ decaySM(i,4)+decaySM(i,5)+decaySM(i,6)+decaySM(i,7)

        }
//c      write(*,*) hmass(mh(1)), decaySM(1,1),decaySM(1,2),decaySM(1,3), decaySM(1,4)
//c      write(*,*) decaySM(1,5) , decaySM(1,6) , decaySM(1,7)
//c      write(*,*) totaldecaySM(1)
//c      write(*,*) decaySM(1,1)+decaySM(1,2)+decaySM(1,3)+ decaySM(1,4)+decaySM(1,5)+decaySM(1,6)+decaySM(1,7)
//c      write(*,*) totaldecaySM

//cccccccccccccccccccccccccccccccccccccc

    ch1check = 1;
    ch2check = 1;

//c checks to see if H1+ production *BF(H+ -> t b~) is interesting
        if ((chmass[mch[1]].gt.180. ).and.(chmass[mch[1]].lt.1000. )) {
            interestingcut = lineint(chmass[mch[1]],interestingmass,interestingcross,6);
            if ((mch_tb_prod[1]*BFch[1][2]) .lt. interestingcut) {
                ch1check = 0;
            }
        }

//c checks to see if H2+ production *BF(H2+ -> t b~) is interesting
        if ((chmass[mch[2]].gt.180. ).and.(chmass[mch[2]].lt.1000. )) {
            interestingcut = lineint(chmass[mch[2]],interestingmass,interestingcross,6);
            if ((mch_tb_prod[2]*BFch[2][2]*(1. -BFt[1])*(1. -BFt[1])) .lt. interestingcut) {
                ch2check = 0;
            }
        }

//c kills the point if neither H1+ or H2+ production*BF isn't interesting
//c      if ((ch1check .eq. 0).and.(ch2check .eq. 0)) then
//c        happy1 = 0
//c      endif


    return;
    }
