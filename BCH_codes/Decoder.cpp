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
    error_location_polynomial[0][0] = 0;		/* index form */
    error_location_polynomial[1][0] = 1;		/* polynomial form */
    for (i = 1; i < t2; i++) {
        error_location_polynomial[0][i] = -1;	/* index form */
        error_location_polynomial[1][i] = 0;	/* polynomial form */
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


        do {
            u++;
            if (x[u] == -1) {
                gf_field->cos_tam(l, u, error_location_polynomial);
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
                gf_field->form_new_elp(u, q, error_location_polynomial, l, t2, x);
            }
            u_lu[u + 1] = u - l[u + 1];

            /* form (u+1)th discrepancy */
            if (u < t2) {
                gf_field->form_discrepancy(s, u, x, l, error_location_polynomial);
            }
        } while ((u < t2) && (l[u + 1] <= gf_field->get_error_code_capability()));

        u++;
        if (l[u] <= gf_field->get_error_code_capability()) {/* Can correct errors */
            int_correct_errors(error_location_polynomial, u, l, q, cx_coefficients);
        }
    }
}

void Decoder::int_correct_errors(int elp[][1024], int u, int *l, int q, int *cx_coefficients)
{
    int             root[200], loc[200], err[1024], reg[201];
    int count = 0;
    int j = 0;
    int i = 0;
    int *index_of = gf_field->get_index_of();
    int *alpha = gf_field->get_alpha();
    /* put elp into index form */
    for (i = 0; i <= l[u]; i++)
        elp[u][i] = index_of[elp[u][i]];



    /* Chien search: find roots of the error location polynomial */
    for (i = 1; i <= l[u]; i++)
        reg[i] = elp[u][i];
    count = 0;
    for (i = 1; i <= gf_field->get_n(); i++) {
        q = 1;
        for (j = 1; j <= l[u]; j++)
            if (reg[j] != -1) {
                reg[j] = (reg[j] + j) % gf_field->get_n();
                q ^= alpha[reg[j]];
            }
        if (!q) {	/* store root and error
                     * location number indices */
            root[count] = i;
            loc[count] = gf_field->get_n() - i;
            count++;
        }
    }
    printf("\n");
    printf("sigma(x) = ");
    for (i = 0; i <= l[u]; i++)
        printf("%3d ", elp[u][i]);
    printf("\n");
    printf("Roots: ");
    for (i = 0; i < count; ++i)
        printf("%3d ", loc[i]);
    printf("\n");
    if (count == l[u])
        /* no. roots = degree of elp hence <= t errors */
        for (i = 0; i < l[u]; i++)
            cx_coefficients[loc[i]] ^= 1;
    else	/* elp has degree >t hence cannot solve */
        printf("Incomplete decoding: errors detected\n");
}