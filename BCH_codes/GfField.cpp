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
{
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

void GfField::int_generate_cycle_sets_mod_n(int *size)
{
  int test;
  int ll = 0;
  int jj = 1;			/* cycle set index */
  
  size[0] = 1;
  size[1] = 1;
//  if (polynomial->get_degree() > 9)  {
//    printf("Computing cycle sets modulo %d\n", polynomial->get_n());
//    printf("(This may take some time)...\n");
//  }
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

void GfField::print_primitive_polynomial()
{
    polynomial->print_primitive_polynomial();
}
