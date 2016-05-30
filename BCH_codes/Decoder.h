//
// Created by marcin on 19.05.16.
//

#ifndef BCH_DECODER_H
#define BCH_DECODER_H

#include "GfField.h"

class Decoder {
private:
    int i = 0;
    int j = 0;
    int u = 0;
    int q = 0;
    int             error_location_polynomial[1026][1024], x[1026], l[1026], u_lu[1026];
    int s[1025];
    int t2;
    GfField *gf_field;
    void int_correct_errors(int elp[][1024], int u, int *l, int q, int *cx_coefficients);
public:
    Decoder(GfField *gf);
    void form_syndromes(int* cx_coefficients);
    void try_to_correct_errors(int* cx_coefficients);
};


#endif //BCH_DECODER_H
