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
    int k = 0; /* dafuq is this */
    int d = 0; /* jw */
    int error_code_capability = 0;
    Polynomial *polynomial;
    int syn_error = 0;
public:
    GfField(int poly_degree);
    int* get_generated_poly();
    int get_k();
    int is_syn_error();
    int get_error_code_capability();
    int get_code_length();
    int* get_polynomial_form();
    void generate_gf();
    void gen_poly(int capability = 0);
    void form_syndromes(int *s, int *cx_coefficients);
    void cos_tam(int *l, int u, int elp[][1024]);
    void form_new_elp(int u, int q, int elp[][1024], int *l, int t2, int *x);
    void form_discrepancy(int *s, int u, int *x, int *l, int elp[][1024]);
    void store_new_elp(int *l, int q, int u);
    void correct_errors(int elp[][1024], int u, int *l, int q, int *cx_coefficients);
};


#endif //BCH_GFFIELD_H
