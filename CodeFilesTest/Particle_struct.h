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
#include <stdlib.h>
#include <math.h>

typedef struct {
    double mass; //The mass eigenvalue of this particle state.
    double Y_d; //The effective Yukawa couplings for down with this particle.
    double Y_u; //The effective Yukawa couplings for up with this particle.
    double Y_s; //The effective Yukawa couplings for strange with this particle.
    double Y_c; //The effective Yukawa couplings for charm with this particle.
    double Y_b; //The effective Yukawa couplings for bottom with this particle.
    double Y_t; //The effective Yukawa couplings for top with this particle.
    double Y_tau; //The effective Yukawa couplings for tau with this particle.
    double Y_mu; //The effective Yukawa couplings for muon with this particle.
    int evec_size;  //This tells us the number of elements in the eigenvector.
    double* evec; //The vector composition of this particle state.
    int bf_size; //This tells us the number of decay modes of the given particle.
    int decay_size; //This tells us the number of decay modes of the given particle.
    double* decay;  //The array containing the various decay modes of the particle.
    double* branching_frac;  //The array containing the various decay modes of the particle.
} particle;

void calc_yukawa (particle *Higgs, double *Y1, double *Y2, double *Y3, double *YL);
    /*
     Higgs   input/output   Pointer to the struct for the particle.
     Y1      input   Pointer to the array for the yukawas for 1st family, down/up.
     Y2      input   Pointer to the array for the yukawas for 2nd family, strange/charm.
     Y3      input   Pointer to the array for the yukawas for 3rd family, bottom/top.
     */

void fill_struct (const double m, double *A, double *Y1, double *Y2, double *Y3, double *YL, particle *Higgs);
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
