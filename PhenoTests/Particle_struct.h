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

typedef struct {
    double mass; //The mass eigenvalue of this particle state.
    double Y_d; //The effective Yukawa couplings for down with this particle.
    double Y_u; //The effective Yukawa couplings for up with this particle.
    double Y_s; //The effective Yukawa couplings for strange with this particle.
    double Y_c; //The effective Yukawa couplings for charm with this particle.
    double Y_b; //The effective Yukawa couplings for bottom with this particle.
    double Y_t; //The effective Yukawa couplings for top with this particle.
    int evec_size;  //This tells us the number of elements in the eigenvector.
    double* eigenvec; //The vector composition of this particle state.
    int bf_size; //This tells us the number of decay modes of the given particle.
    double* branching_frac;  //The array containing the various decay modes of the particle.
} particle;

//Thing to remember: there will need to be a loop creating the array of particles after we cal the eigenvector/eigenvalue function.

void calc_yukawa (particle *Higgs, double *Y1, double *Y2, double *Y3);
    /*
     Higgs   input/output   Pointer to the struct for the particle.
     Y1      input   Pointer to the array for the yukawas for 1st family, down/up.
     Y2      input   Pointer to the array for the yukawas for 2nd family, strange/charm.
     Y3      input   Pointer to the array for the yukawas for 3rd family, bottom/top.
     */

void fill_struct (const double m, double *A, double *Y1, double *Y2, double *Y3, particle *Higgs);
/*
 m      input   The eigenvalue of the mass state.
 A      input   Pointer to the eigenvector of the mass state.
 Y1     input   Pointer to the array holding the first family Yukawa couplings.
 Y2     input   Pointer to the array holding the second family Yukawa couplings.
 Y3     input   Pointer to the array holding the third family Yukawa couplings.
 Higgs  output  Pointer to the struct for the particle in question.
 */

double dot_prod (double *A, double *B);


#endif /* Particle_struct_h */
