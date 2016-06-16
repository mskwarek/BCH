//
// Created by marcin on 30.05.16.
//

#ifndef BCH_GFFIELD_H
#define BCH_GFFIELD_H

//
// Created by marcin on 01.05.16.
//

#include "Polynomial.h"

class GfField
{
private:
    int             alpha_to[1048576], index_of[1048576], g[548576];
    int cycle[1024][21];
    int k = 0; /* dafuq is this */
    int d = 0; /* jw */
    int rdncy = 0;
    int nocycles = 0;
    int error_code_capability = 0;
    Polynomial *polynomial;
    void int_generate_next_cycle_set(int *size, int cycle_set_index);
    void int_compute_generator_polynomial(int* zeros);
    void int_search_for_cycle_sets_roots(int *size, int* zeros);
    void int_generate_cycle_sets_mod_n(int *size);
public:
    GfField(int poly_degree);
    int* get_generated_poly();
    int get_k();
    int get_n();
    int get_error_code_capability();
    int get_code_length();
    int* get_index_of();
    int* get_alpha();
    int* get_polynomial_form();
    void generate_gf();
    void gen_poly(int capability = 0);
    void print_primitive_polynomial();
    void print_generator_polynomial();
    void print_error_code_capability();
    void print_bch_code_features();
};


#endif //BCH_GFFIELD_H
