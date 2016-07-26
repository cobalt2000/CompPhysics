//
//  Particle_struct.h
//  
//
//  Created by Valerie Plaus on 7/19/16.
//
//

#ifndef Particle_struct_h
#define Particle_struct_h

#include <stdio.h>

void Calc_Yukawa (double *evec, double *Y, double *Y1, double *Y2);
/*
 evec   input   Pointer to the array for the eigenvector composition of the particle.
 Y      input   Pointer to the array for the yukawas for whichever family, up/down, charm/strange, top/bottom (first, second, and third families).
 Y1     output  Places the effective Yukawa coupling for the down-type quark in the struct for the particle.
 Y2     output  Places the effective Yukawa coupling for the up-type quark in the struct for the particle.
 */

void Fill_Struct (const double m, double *A, double *Y1, double *Y2, double *Y3, Particle *Higgs);
/*
 m      input   The eigenvalue of the mass state.
 A      input   Pointer to the eigenvector of the mass state.
 Y1     input   Pointer to the array holding the first family Yukawa couplings.
 Y2     input   Pointer to the array holding the second family Yukawa couplings.
 Y3     input   Pointer to the array holding the third family Yukawa couplings.
 Higgs  output  Pointer to the struct for the particle in question.
 */


#endif /* Particle_struct_h */
