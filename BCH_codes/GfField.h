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
    int syn_error = 0;
    void int_generate_next_cycle_set(int *size, int cycle_set_index);
    void int_compute_generator_polynomial(int* zeros);
    void int_search_for_cycle_sets_roots(int *size, int* zeros);
    void int_generate_cycle_set();
    void int_generate_cycle_sets_mod_n(int *size);
public:
    GfField(int poly_degree);
    int* get_generated_poly();
    int get_k();
    int get_n();
    int is_syn_error();
    int get_error_code_capability();
    int get_code_length();
    int* get_index_of();
    int* get_alpha();
    int* get_polynomial_form();
    void generate_gf();
    void gen_poly(int capability = 0);
    void form_syndromes(int *s, int *cx_coefficients);
    void cos_tam(int *l, int u, int elp[][1024]);
    void form_new_elp(int u, int q, int elp[][1024], int *l, int t2, int *x);
    void form_discrepancy(int *s, int u, int *x, int *l, int error_location_polynomial[][1024]);
    void store_new_elp(int *l, int q, int u);
    void print_generator_polynomial();
    void print_error_code_capability();
    void print_bch_code_features();
    void print_syndromes_features(int *s);
};


#endif //BCH_GFFIELD_H
