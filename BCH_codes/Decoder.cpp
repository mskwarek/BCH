//
// Created by marcin on 19.05.16.
//

#include "Decoder.h"

Decoder::Decoder(GfField *gf)
{
    gf_field = gf;
    t2 = gf_field->get_error_code_capability() *2;
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
}

void Decoder::form_syndromes(int* cx_coefficients)
{
    gf_field->form_syndromes(s, cx_coefficients);
}

void Decoder::try_to_correct_errors(int* cx_coefficients)
{
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