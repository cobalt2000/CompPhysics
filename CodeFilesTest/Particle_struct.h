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

void calc_yukawa (struct particle *Higgs, double *Y1, double *Y2, double *Y3, double *YL);
    /*
     Higgs   input/output   Pointer to the struct for the particle.
     Y1      input   Pointer to the array for the yukawas for 1st family, down/up.
     Y2      input   Pointer to the array for the yukawas for 2nd family, strange/charm.
     Y3      input   Pointer to the array for the yukawas for 3rd family, bottom/top.
     */

void fill_struct (const double m, double *A, double *Y1, double *Y2, double *Y3, double *YL, struct particle *Higgs);
/*
 m      input   The eigenvalue of the mass state.
 A      input   Pointer to the eigenvector of the mass state.
 Y1     input   Pointer to the array holding the first family Yukawa couplings.
 Y2     input   Pointer to the array holding the second family Yukawa couplings.
 Y3     input   Pointer to the array holding the third family Yukawa couplings.
 Higgs  output  Pointer to the struct for the particle in question.
 */

double dot_prod (double *A, double *B);
/*
The two vectors, A and B, are contracted using the standard inner/dot product. 
If the vectors aren't the same length, returns print-to-screen statement and "0".
*/

double W_H_rotation (double *A, double *B);
/*
 This is just to clarify the basis difference between the CP even and charged Higgs couplings to the W.
 This matches the MSSM Feynman rules in https://arxiv.org/pdf/hep-ph/9511250.pdf
*/

#endif /* Particle_struct_h */
