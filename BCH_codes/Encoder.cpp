//
// Created by marcin on 19.05.16.
//

#include "Encoder.h"

Encoder::Encoder(GfField *gf)
{
    gf_field = gf;
}

void Encoder::encode(int *source_data, int *coding_result)
{
    register int    i;
    int k = gf_field->get_k();


    for (i = 0; i < gf_field->get_code_length() - k; i++)
        coding_result[i] = 0;
    for (i = k - 1; i >= 0; i--) {
        if (int_is_source_different_from_coding_result(source_data[i], coding_result[gf_field->get_code_length() - k - 1]))
        {
            int_code_if_different(coding_result);
        }
        else
        {
            int_code_if_compatible(coding_result);
        }
    }
}

bool Encoder::int_is_source_different_from_coding_result(int source_data, int coding_result)
{
    register int feedback = source_data ^ coding_result;
    return (bool)(feedback != 0);
}

void Encoder::int_code_if_different(int* coding_result)
{
    int* g = gf_field->get_generated_poly();
    int k = gf_field->get_k();
    int j = 0;
    for (j = gf_field->get_code_length() - k - 1; j > 0; j--)
        if (g[j] != 0)
            coding_result[j] = coding_result[j - 1] ^ 1;
        else
            coding_result[j] = coding_result[j - 1];
    coding_result[0] = g[0] && 1;
}

void Encoder::int_code_if_compatible(int* coding_result)
{
    int k = gf_field->get_k();
    int j = 0;
    for (j = gf_field->get_code_length() - k - 1; j > 0; j--)
        coding_result[j] = coding_result[j - 1];
    coding_result[0] = 0;
}