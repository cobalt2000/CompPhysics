//
//  Particle_struct.c
//  
//
//  Created by Valerie Plaus on 7/19/16.
//
//

#include "Particle_struct.h"

datatype struct {
    double  mass; //The mass eigenvalue of this particle state.
    double* eigenvec; //The vector composition of this particle state.
    double Y_u; //The effective Yukawa couplings for up with this particle.
    double Y_d; //The effective Yukawa couplings for down with this particle.
    double Y_c; //The effective Yukawa couplings for charm with this particle.
    double Y_s; //The effective Yukawa couplings for strange with this particle.
    double Y_t; //The effective Yukawa couplings for top with this particle.
    double Y_b; //The effective Yukawa couplings for bottom with this particle.
} particle;

datatype struct { // These would either be input by the user or randomized.
    double* Y_ud; //The Yukawa couplings for the first family, up and down, for each Higgs field.
    double* Y_cs; //The Yukawa couplings for the second family, the charm and strange quarks, for each Higgs field.
    double* Y_tb;  //The  Yukawa couplings for the third family, the top and bottom quarks, for each Higgs field.
} Yukawa;
