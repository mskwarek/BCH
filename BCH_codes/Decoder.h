//
// Created by marcin on 19.05.16.
//

#ifndef BCH_DECODER_H
#define BCH_DECODER_H

#include "GfField.h"
#include <memory>

class Decoder {
private:
    int step_number = 0;
    int q = 0;
    int             error_location_polynomial[1026][1024], x[1026], elp_degree[1026], u_lu[1026];
    int root_counter = 0;
    int s[1025];
    int root[200], loc[200];
    int t2;
    std::shared_ptr<GfField> gf_field;
    void int_compute_error_location_polynomial();
    void int_init_decoder_variables();
    void int_search_for_greatest_one();
    bool int_can_correct_errors();
    void int_correct_errors(int *cx_coefficients);
    int* int_put_elp_into_index_form();
    void int_find_elp_roots();
    void int_correct_error_bits(int *cs_coefficients);
    bool int_is_syndrom_null();
    void int_cos_tam();
    void int_form_new_elp();
    void int_form_discrepancy();
    void int_store_new_elp();
public:
    Decoder(std::shared_ptr<GfField> gf);
    void print_sigma();
    void print_roots();
    void form_syndromes(int* cx_coefficients);
    void try_to_correct_errors(int* cx_coefficients);
    void print_syndromes_features();
};


#endif //BCH_DECODER_H
