//
// Created by marcin on 01.05.16.
//

#ifndef BCH_CODER_H
#define BCH_CODER_H

#include "Decoder.h"

class Coder {
private:
    GfField *gf_field;
    Decoder *decoder;
    bool int_is_source_different_from_coding_result(int source_data, int coding_result);
    void int_code_if_different(int* coding_result);
    void int_code_if_compatible(int* coding_result);
public:
    Coder();
    ~Coder();
    int* get_generated_polynomial();
    int* get_polynomial_form();
    int get_code_length();
    void encode_bch(int *source_data, int *coding_result);
    void decode_bch(int* const cx_coefficients);
    void code_polynomial(int* source_data, int* bb, int* destination);
    void print_decoding_result(int* data, int* recd);
};


#endif //BCH_CODER_H
