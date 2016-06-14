#include "GfField.h"
#include <iostream>

GfField::GfField(int poly_degree)
{
    polynomial = new Polynomial();
    polynomial->read_p(poly_degree);
    
    cycle[0][0] = 0;
    cycle[1][0] = 1;
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
    int mask;
    int m = polynomial->get_degree();
    mask = 1;
    alpha_to[m] = 0;
    for (int i = 0; i < m; i++) {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (polynomial->get_form()[i] != 0)
            alpha_to[m] ^= mask;
        mask <<= 1;
    }
    index_of[alpha_to[m]] = m;
    mask >>= 1;
    for (int i = m + 1; i < polynomial->get_n(); i++) {
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
    register int	root;
    int             zeros[1024];
    int size[1024];

    error_code_capability = capability;
    d = 2 * error_code_capability + 1;

    int_generate_cycle_sets_mod_n(size);
    int_search_for_cycle_sets_roots(size, zeros);

    k = polynomial->get_code_length() - rdncy;
    if (k<0)
    {
        throw 1;
    }

    int_compute_generator_polynomial(zeros);
}

/* Generate cycle sets modulo n, n = 2**m - 1 */
void GfField::int_generate_cycle_sets_mod_n(int *size)
{
  int test;
  int ll = 0;
  int jj = 1;			/* cycle set index */
  
  size[0] = 1;
  size[1] = 1;
  if (polynomial->get_degree() > 9)  {
    printf("Computing cycle sets modulo %d\n", polynomial->get_n());
    printf("(This may take some time)...\n");
  }
  do {
    ll = 0;
    int_generate_next_cycle_set(size, jj);
    /* Next cycle set representative */
    do {
      ll++;
      test = 0;
      for (int ii = 1; ((ii <= jj) && (!test)); ii++)
	/* Examine previous cycle sets */
	for (int kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
	  if (ll == cycle[ii][kaux])
	    test = 1;
    } while (test && (ll < (polynomial->get_n() - 1)));
    if (!(test)) {
      jj++;	/* next cycle set index */
      cycle[jj][0] = ll;
      size[jj] = 1;
    }
  } while (ll < (polynomial->get_n() - 1));
  nocycles = jj;		/* number of cycle sets modulo n */  
}

/* Generate the jj-th cycle set */
void GfField::int_generate_next_cycle_set(int *size, int cycle_set_index)
{
  int ii = 0;
  int aux;
  do {
    ii++;
    cycle[cycle_set_index][ii] = (cycle[cycle_set_index][ii - 1] * 2) % polynomial->get_n();
    size[cycle_set_index]++;
    aux = (cycle[cycle_set_index][ii] * 2) % polynomial->get_n();
  } while (aux != cycle[cycle_set_index][0]);
}

/* Search for roots 1, 2, ..., d-1 in cycle sets */
void GfField::int_search_for_cycle_sets_roots(int* size, int* zeros)
{
  int  min[1024];
  int  noterms;
  int kaux = 0;
  int test;

  for (int ii = 1; ii <= nocycles; ii++) {
    min[kaux] = 0;
    test = 0;
    for (int jj = 0; ((jj < size[ii]) && (!test)); jj++)
      for (int root = 1; ((root < d) && (!test)); root++)
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
  for (int ii = 0; ii < noterms; ii++)
    for (int jj = 0; jj < size[min[ii]]; jj++) {
      zeros[kaux] = cycle[min[ii]][jj];
      kaux++;
    }
}

void GfField::int_compute_generator_polynomial(int *zeros)
{
    g[0] = alpha_to[zeros[1]];
    g[1] = 1;		/* g(x) = (X + zeros[1]) initially */
    for (int ii = 2; ii <= rdncy; ii++) {
        g[ii] = 1;
        for (int jj = ii - 1; jj > 0; jj--)
            if (g[jj] != 0)
                g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % polynomial->get_n()];
            else
                g[jj] = g[jj - 1];
        g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % polynomial->get_n()];
    }
}

void GfField::print_error_code_capability()
{
  std::cout<<"error correcting capability, t: "<<error_code_capability<<std::endl;
}

void GfField::print_bch_code_features()
{
  printf("This is a (%d, %d, %d) binary BCH code\n", polynomial->get_code_length(), k, d);
}

void GfField::print_generator_polynomial()
{
  printf("Generator polynomial:\ng(x) = ");
  for (int ii = 0; ii <= rdncy; ii++) {
    printf("%d", g[ii]);
    if (ii && ((ii % 50) == 0))
      printf("\n");
  }
  printf("\n");
}

void GfField::form_syndromes(int *s, int *cx_coefficients)
{
    int t2 = 2 * error_code_capability;
    /* first form the syndromes */
    for (int i = 1; i <= t2; i++) {
        s[i] = 0;
        for (int j = 0; j < polynomial->get_code_length(); j++)
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
    }
}

void GfField::print_syndromes_features(int *s)
{
  std::cout<<"t2 = "<<2 * error_code_capability<<std::endl;  
  std::cout<<"S(x) = ";
  for (int i=0; i < 2 * error_code_capability; i++)
    std::cout<<s[i];
  std::cout<<std::endl;
  
}

void GfField::cos_tam(int *l, int u, int elp[][1024])
{
    l[u + 1] = l[u];
    for (int i = 0; i <= l[u]; i++) {
        elp[u + 1][i] = elp[u][i];
        elp[u][i] = index_of[elp[u][i]];
    }
}

void GfField::form_new_elp(int u, int q, int elp[][1024], int *l, int t2, int *x)
{
    /* form new elp(x) */
    for (int i = 0; i < t2; i++)
        elp[u + 1][i] = 0;
    for (int i = 0; i <= l[q]; i++)
        if (elp[q][i] != -1)
	  elp[u + 1][i + u - q] =
	    alpha_to[(x[u] + polynomial->get_n() - x[q] + elp[q][i]) % polynomial->get_n()];
    for (int i = 0; i <= l[u]; i++) {
        elp[u + 1][i] ^= elp[u][i];
        elp[u][i] = index_of[elp[u][i]];
    }
}

void GfField::form_discrepancy(int *s, int u, int *x, int *l, int error_location_polynomial[][1024])
{
    /* no discrepancy computed on last iteration */
    if (s[u + 1] != -1)
        x[u + 1] = alpha_to[s[u + 1]];
    else
        x[u + 1] = 0;
    for (int i = 1; i <= l[u + 1]; i++)
        if ((s[u + 1 - i] != -1) && (error_location_polynomial[u + 1][i] != 0))
            x[u + 1] ^= alpha_to[(s[u + 1 - i]
                                  + index_of[error_location_polynomial[u + 1][i]]) % polynomial->get_n()];
    /* put d[u+1] into index form */
    x[u + 1] = index_of[x[u + 1]];
}

/*
 * have now found q such that d[u]!=0 and
 * u_lu[q] is maximum
 */
/* store degree of new elp polynomial */
void GfField::store_new_elp(int *l, int q, int u)
{
    if (l[u] > l[q] + u - q)
        l[u + 1] = l[u];
    else
        l[u + 1] = l[q] + u - q;
}

