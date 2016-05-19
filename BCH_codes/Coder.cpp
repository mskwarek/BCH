//
// Created by marcin on 01.05.16.
//

#include "Coder.h"

Coder::Coder()
{
    gf_field = new GfField(10);
    gf_field->generate_gf();
    gf_field->gen_poly(15);
}

Coder::~Coder()
{
    delete gf_field;
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
    register int    i, j;
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

bool Coder::int_is_source_different_from_coding_result(int source_data, int coding_result)
{
    register int feedback = source_data ^ coding_result;
    return (bool)(feedback != 0);
}

void Coder::int_code_if_different(int* coding_result)
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

void Coder::int_code_if_compatible(int* coding_result)
{
    int k = gf_field->get_k();
    int j = 0;
    for (j = gf_field->get_code_length() - k - 1; j > 0; j--)
        coding_result[j] = coding_result[j - 1];
    coding_result[0] = 0;
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
    register int    i, j, u, q, t2;
    int             elp[1026][1024], x[1026], l[1026], u_lu[1026];
    int s[1025];

    t2 = gf_field->get_error_code_capability() *2;
    gf_field->form_syndromes(s, cx_coefficients);

    if (gf_field->is_syn_error())
    {	/* if there are errors, try to correct them */
        /*
         * Compute the error location polynomial via the Berlekamp
         * iterative algorithm. Following the terminology of Lin and
         * Costello's book :   d[u] is the 'mu'th discrepancy, where
         * u='mu'+1 and 'mu' (the Greek letter!) is the step number
         * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
         * the elp at that step, and u_l[u] is the difference between
         * the step number and the degree of the elp.
         */
        /* initialise table entries */

        x[0] = 0;			/* index form */
        x[1] = s[1];		/* index form */
        elp[0][0] = 0;		/* index form */
        elp[1][0] = 1;		/* polynomial form */
        for (i = 1; i < t2; i++) {
            elp[0][i] = -1;	/* index form */
            elp[1][i] = 0;	/* polynomial form */
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        do {
            u++;
            if (x[u] == -1) {
                gf_field->cos_tam(l, u, elp);
            } else
                /*
                 * search for words with greatest u_lu[q] for
                 * which d[q]!=0
                 */
            {
                q = u - 1;
                while ((x[q] == -1) && (q > 0))
                    q--;
                /* have found first non-zero d[q]  */
                if (q > 0) {
                    j = q;
                    do {
                        j--;
                        if ((x[j] != -1) && (u_lu[q] < u_lu[j]))
                            q = j;
                    } while (j > 0);
                }

                gf_field->store_new_elp(l, q, u);
                gf_field->form_new_elp(u, q, elp, l, t2, x);
            }
            u_lu[u + 1] = u - l[u + 1];

            /* form (u+1)th discrepancy */
            if (u < t2) {
                gf_field->form_discrepancy(s, u, x, l, elp);
            }
        } while ((u < t2) && (l[u + 1] <= gf_field->get_error_code_capability()));

        u++;
        if (l[u] <= gf_field->get_error_code_capability()) {/* Can correct errors */
            gf_field->correct_errors(elp, u, l, q, cx_coefficients);
        }
    }
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

void Coder::print_decoding_result(int* data, int* recd)
{
    int i = 0;
    int decerror = 0;
    for (i = gf_field->get_code_length() - gf_field->get_k(); i < gf_field->get_code_length(); i++)
        if (data[i - gf_field->get_code_length() + gf_field->get_k()] != recd[i])
            decerror++;
    if (decerror)
        printf("There were %d decoding errors in message positions\n", decerror);
    else
        printf("Succesful decoding\n");
}
