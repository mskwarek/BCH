//
// Created by marcin on 01.05.16.
//

#ifndef BCH_CODER_H
#define BCH_CODER_H

#include "GfField.cpp"

class Coder {
public:
    void encode_bch(GfField* gf_field, int *source_data, int *coding_result);
    void decode_bch(GfField* gf, int *cx_coefficients);
    void code_polynomial(int code_length, int k, int* source_data, int* bb, int* destination);
};


#endif //BCH_CODER_H
