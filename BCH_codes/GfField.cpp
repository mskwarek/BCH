#include "GfField.h"
#include <iostream>

GfField::GfField(int poly_degree)
{
    polynomial = new Polynomial();
    polynomial->read_p(poly_degree);
}

int* GfField::get_generated_poly()
{
    return g;
}

int GfField::get_k()
{
    return k;
}

int GfField::is_syn_error()
{
    return syn_error;
}

int GfField::get_error_code_capability()
{
    return error_code_capability;
}

int GfField::get_code_length()
{
    return polynomial->get_code_length();
}

int* GfField::get_polynomial_form()
{
    return polynomial->get_form();
}

int GfField::get_n()
{
    return polynomial->get_n();
}

int* GfField::get_index_of()
{
    return index_of;
}

int* GfField::get_alpha()
{
    return alpha_to;
}

void GfField::generate_gf()
/*
 * Generate field GF(2**m) from the irreducible polynomial p(X) with
 * coefficients in p[0]..p[m].
 *
 * Lookup tables:
 *   index->polynomial form: alpha_to[] contains j=alpha^i;
 *   polynomial form -> index form:	index_of[j=alpha^i] = i
 *
 * alpha=2 is the primitive element of GF(2**m)
 */
{
    register int    i, mask;
    int m = polynomial->get_degree();
    mask = 1;
    alpha_to[m] = 0;
    for (i = 0; i < m; i++) {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (polynomial->get_form()[i] != 0)
            alpha_to[m] ^= mask;
        mask <<= 1;
    }
    index_of[alpha_to[m]] = m;
    mask >>= 1;
    for (i = m + 1; i < polynomial->get_n(); i++) {
        if (alpha_to[i - 1] >= mask)
            alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
        else
            alpha_to[i] = alpha_to[i - 1] << 1;
        index_of[alpha_to[i]] = i;
    }
    index_of[0] = -1;
}

void GfField::gen_poly(int capability)
/*
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (x+alpha^i), for every i in
 * the above cycle sets.
 */
{
    register int	ii, jj, ll, kaux;
    register int	test, aux, nocycles, root, noterms, rdncy;
    int             cycle[1024][21], size[1024], min[1024], zeros[1024];

    /* Generate cycle sets modulo n, n = 2**m - 1 */
    cycle[0][0] = 0;
    size[0] = 1;
    cycle[1][0] = 1;
    size[1] = 1;
    jj = 1;			/* cycle set index */
    if (polynomial->get_degree() > 9)  {
        printf("Computing cycle sets modulo %d\n", polynomial->get_n());
        printf("(This may take some time)...\n");
    }
    do {
        /* Generate the jj-th cycle set */
        ii = 0;
        do {
            ii++;
            cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % polynomial->get_n();
            size[jj]++;
            aux = (cycle[jj][ii] * 2) % polynomial->get_n();
        } while (aux != cycle[jj][0]);
        /* Next cycle set representative */
        ll = 0;
        do {
            ll++;
            test = 0;
            for (ii = 1; ((ii <= jj) && (!test)); ii++)
                /* Examine previous cycle sets */
                for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
                    if (ll == cycle[ii][kaux])
                        test = 1;
        } while ((test) && (ll < (polynomial->get_n() - 1)));
        if (!(test)) {
            jj++;	/* next cycle set index */
            cycle[jj][0] = ll;
            size[jj] = 1;
        }
    } while (ll < (polynomial->get_n() - 1));
    nocycles = jj;		/* number of cycle sets modulo n */

    error_code_capability = capability;
    if(capability == 0) {
        printf("Enter the error correcting capability, t: ");
        scanf("%d", &error_code_capability);
    }
    std::cout<<"error correcting capability, t: "<<error_code_capability<<std::endl;

    d = 2 * error_code_capability + 1;

    /* Search for roots 1, 2, ..., d-1 in cycle sets */
    kaux = 0;
    rdncy = 0;
    for (ii = 1; ii <= nocycles; ii++) {
        min[kaux] = 0;
        test = 0;
        for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
            for (root = 1; ((root < d) && (!test)); root++)
                if (root == cycle[ii][jj])  {
                    test = 1;
                    min[kaux] = ii;
                }
        if (min[kaux]) {
            rdncy += size[min[kaux]];
            kaux++;
        }
    }
    noterms = kaux;
    kaux = 1;
    for (ii = 0; ii < noterms; ii++)
        for (jj = 0; jj < size[min[ii]]; jj++) {
            zeros[kaux] = cycle[min[ii]][jj];
            kaux++;
        }

    k = polynomial->get_code_length() - rdncy;

    if (k<0)
    {
        printf("Parameters invalid!\n");
        return;
    }

    printf("This is a (%d, %d, %d) binary BCH code\n", polynomial->get_code_length(), k, d);

    /* Compute the generator polynomial */
    g[0] = alpha_to[zeros[1]];
    g[1] = 1;		/* g(x) = (X + zeros[1]) initially */
    for (ii = 2; ii <= rdncy; ii++) {
        g[ii] = 1;
        for (jj = ii - 1; jj > 0; jj--)
            if (g[jj] != 0)
                g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % polynomial->get_n()];
            else
                g[jj] = g[jj - 1];
        g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % polynomial->get_n()];
    }
    printf("Generator polynomial:\ng(x) = ");
    for (ii = 0; ii <= rdncy; ii++) {
        printf("%d", g[ii]);
        if (ii && ((ii % 50) == 0))
            printf("\n");
    }
    printf("\n");
}

void GfField::form_syndromes(int *s, int *cx_coefficients)
{
    int i = 0;
    int j = 0;
    register int t2 = 0;
    t2 = 2 * error_code_capability;
    printf("t2 = %d", t2);
    /* first form the syndromes */
    printf("S(x) = ");
    for (i = 1; i <= t2; i++) {
        s[i] = 0;
        for (j = 0; j < polynomial->get_code_length(); j++)
            if (cx_coefficients[j] != 0)
                s[i] ^= alpha_to[(i * j) % polynomial->get_n()];
        if (s[i] != 0)
            syn_error = 1; /* set error flag if non-zero syndrome */
/*
* Note:    If the code is used only for ERROR DETECTION, then
*          exit program here indicating the presence of errors.
*/
        /* convert syndrome from polynomial form to index form  */
        s[i] = index_of[s[i]];
        printf("%3d ", s[i]);
    }
    printf("\n");
}

void GfField::cos_tam(int *l, int u, int elp[][1024])
{
    int i = 0;
    l[u + 1] = l[u];
    for (i = 0; i <= l[u]; i++) {
        elp[u + 1][i] = elp[u][i];
        elp[u][i] = index_of[elp[u][i]];
    }
}

void GfField::form_new_elp(int u, int q, int elp[][1024], int *l, int t2, int *x)
{
    int i = 0;
    /* form new elp(x) */
    for (i = 0; i < t2; i++)
        elp[u + 1][i] = 0;
    for (i = 0; i <= l[q]; i++)
        if (elp[q][i] != -1)
            elp[u + 1][i + u - q] =
                    alpha_to[(x[u] + polynomial->get_n() - x[q] + elp[q][i]) % polynomial->get_n()];
    for (i = 0; i <= l[u]; i++) {
        elp[u + 1][i] ^= elp[u][i];
        elp[u][i] = index_of[elp[u][i]];
    }
}

void GfField::form_discrepancy(int *s, int u, int *x, int *l, int error_location_polynomial[][1024])
{
    int i = 0;
    /* no discrepancy computed on last iteration */
    if (s[u + 1] != -1)
        x[u + 1] = alpha_to[s[u + 1]];
    else
        x[u + 1] = 0;
    for (i = 1; i <= l[u + 1]; i++)
        if ((s[u + 1 - i] != -1) && (error_location_polynomial[u + 1][i] != 0))
            x[u + 1] ^= alpha_to[(s[u + 1 - i]
                                  + index_of[error_location_polynomial[u + 1][i]]) % polynomial->get_n()];
    /* put d[u+1] into index form */
    x[u + 1] = index_of[x[u + 1]];
}

void GfField::store_new_elp(int *l, int q, int u)
{
    /*
             * have now found q such that d[u]!=0 and
             * u_lu[q] is maximum
             */
    /* store degree of new elp polynomial */
    if (l[u] > l[q] + u - q)
        l[u + 1] = l[u];
    else
        l[u + 1] = l[q] + u - q;
}

