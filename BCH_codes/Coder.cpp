//
// Created by marcin on 01.05.16.
//

#include "Coder.h"

Coder::Coder(int polynomial_degree, int error_correct_capability)
{
    gf_field = new GfField(polynomial_degree);
    gf_field->generate_gf();
    gf_field->gen_poly(error_correct_capability);

    decoder = new Decoder(gf_field);
    encoder = new Encoder(gf_field);
}

Coder::~Coder()
{
    delete gf_field;
    delete encoder;
    delete decoder;
}

int* Coder::get_generated_polynomial()
{
    gf_field->get_generated_poly();
}

int* Coder::get_polynomial_form()
{
    return gf_field->get_polynomial_form();
}

int Coder::get_code_length()
{
    gf_field->get_code_length();
}

void Coder::encode_bch(int *source_data, int *coding_result)
/*
 * Compute redundacy bb[], the coefficients of b(x). The redundancy
 * polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
 * by the generator polynomial g(x).
 */
{
    encoder->encode(source_data, coding_result);
}

void Coder::code_polynomial(int* source_data, int* bb, int* destination)
{
    /*
    * recd[] are the coefficients of c(x) = x**(length-k)*data(x) + b(x)
    */
    int code_length = gf_field->get_code_length();
    int k = gf_field->get_k();
    int i = 0;
    for (i = 0; i < code_length - k; i++)
        destination[i] = bb[i];
    for (i = 0; i < k; i++)
        destination[i + code_length - k] = source_data[i];
}

void Coder::decode_bch(int* const cx_coefficients)
/*
* Simon Rockliff's implementation of Berlekamp's algorithm.
*
* Assume we have received bits in recd[i], i=0..(n-1).
*
* Compute the 2*t syndromes by substituting alpha^i into rec(X) and
* evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
* Then we use the Berlekamp algorithm to find the error location polynomial
* elp[i].
*
* If the degree of the elp is >t, then we cannot correct all the errors, and
* we have detected an uncorrectable error pattern. We output the information
* bits uncorrected.
*
* If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
* to get the roots, hence the inverse roots, the error location numbers.
* This step is usually called "Chien's search".
*
* If the number of errors located is not equal the degree of the elp, then
* the decoder assumes that there are more than t errors and cannot correct
* them, only detect them. We output the information bits uncorrected.
*/
{
    decoder->form_syndromes(cx_coefficients);
    decoder->try_to_correct_errors(cx_coefficients);
}


int Coder::get_decoding_error_number(int *data, int *recd)
{
    int decerror = 0;
    for (int i = gf_field->get_code_length() - gf_field->get_k(); i < gf_field->get_code_length(); i++)
        if (data[i - gf_field->get_code_length() + gf_field->get_k()] != recd[i])
            decerror++;
    if (decerror)
        printf("There were %d decoding errors in message positions\n", decerror);
    else
        printf("Succesful decoding\n");
    return decerror;
}

int Coder::get_k()
{
  return gf_field->get_k();
}
