//
// Created by marcin on 01.05.16.
//

#ifndef BCH_CODER_H
#define BCH_CODER_H

#include <memory>
#include "Decoder.h"
#include "Encoder.h"

class Coder {
private:
    std::shared_ptr<GfField> gf_field;
    std::unique_ptr<Decoder> decoder;
    std::unique_ptr<Encoder> encoder;
public:
    Coder(int polynomial_degree, int error_correct_capability);
    ~Coder();
    int* get_generated_polynomial();
    int* get_polynomial_form();
    int get_code_length();
    int get_k();
    void print_primitive_polynomial();
    void print_sigma();
    void print_roots();
    void encode_bch(int *source_data, int *coding_result);
    void decode_bch(int* const cx_coefficients);
    void code_polynomial(int* source_data, int* bb, int* destination);
    int get_decoding_error_number(int *data, int *recd);
    void print_syndromes_features();
    void print_gf_features();
};


#endif //BCH_CODER_H
