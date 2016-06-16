//
// Created by marcin on 19.05.16.
//

#ifndef BCH_ENCODER_H
#define BCH_ENCODER_H

#include "GfField.h"
#include <memory>

class Encoder {
public:
    Encoder(std::shared_ptr<GfField> gf);
    void encode(int *source_data, int *coding_result);
private:
    std::shared_ptr<GfField> gf_field;
    bool int_is_source_different_from_coding_result(int source_data, int coding_result);
    void int_code_if_different(int* coding_result);
    void int_code_if_compatible(int* coding_result);
};


#endif //BCH_ENCODER_H
