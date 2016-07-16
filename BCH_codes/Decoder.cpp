//
// Created by marcin on 19.05.16.
//

#include "Decoder.h"
#include <iostream>

Decoder::Decoder(GfField *gf)
{
    gf_field = gf;
    t2 = gf_field->get_error_code_capability() *2;
}

void Decoder::form_syndromes(int* cx_coefficients)
{
    gf_field->form_syndromes(s, cx_coefficients);
}

void Decoder::try_to_correct_errors(int* cx_coefficients)
{
  int_init_decoder_variables();
  if (gf_field->is_syn_error())
  {
    int_compute_error_location_polynomial();
    if (this->int_can_correct_errors()) 
    {
      int_correct_errors(cx_coefficients);
    }
  }
}

/*
 * Compute the error location polynomial via the Berlekamp
 * iterative algorithm. Following the terminology of Lin and
 * Costello's book :   d[u] is the 'mu'th discrepancy, where
 * u='mu'+1 and 'mu' (the Greek letter!) is the step number
 * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
 * the elp at that step, and u_l[u] is the difference between
 * the step number and the degree of the elp.
 */
void Decoder::int_compute_error_location_polynomial()
{
  do {
    step_number++;
    if (x[step_number] == -1) 
      {
	gf_field->cos_tam(elp_degree, step_number, error_location_polynomial);
      } 
    else
      {
	int_search_for_greatest_one();
      }
    u_lu[step_number + 1] = step_number - elp_degree[step_number + 1];

    /* form (u+1)th discrepancy */
    if (step_number < t2) {
      gf_field->form_discrepancy(s, step_number, x, elp_degree, error_location_polynomial);
    }
  } while ((step_number < t2) && (elp_degree[step_number + 1] <= gf_field->get_error_code_capability()));
  
  step_number++;
}

/*
 * search for words with greatest u_lu[q] for
 * which d[q]!=0
 */
void Decoder::int_search_for_greatest_one()
{
  int j = 0;
  q = step_number - 1;
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
  
  gf_field->store_new_elp(elp_degree, q, step_number);
  gf_field->form_new_elp(step_number, q, error_location_polynomial, elp_degree, t2, x);
}

bool Decoder::int_can_correct_errors()
{
  return elp_degree[step_number] <= gf_field->get_error_code_capability();
}

void Decoder::int_init_decoder_variables()
{
    x[0] = 0;			/* index form */
    x[1] = s[1];		/* index form */
    error_location_polynomial[0][0] = 0;		/* index form */
    error_location_polynomial[1][0] = 1;		/* polynomial form */
    for (int i = 1; i < t2; i++) {
        error_location_polynomial[0][i] = -1;	/* index form */
        error_location_polynomial[1][i] = 0;	/* polynomial form */
    }
    elp_degree[0] = 0;
    elp_degree[1] = 0;
    u_lu[0] = -1;
    u_lu[1] = 0;
    step_number = 0;
}

void Decoder::int_correct_errors(int *cx_coefficients)
{
    int *index_of = int_put_elp_into_index_form();
    int_find_elp_roots(root, loc);

    print_sigma();
    print_roots();
    try
      {
	int_correct_error_bits(cx_coefficients);
      }
    catch(...)
      {
	/* elp has degree >t hence cannot solve */
        printf("Incomplete decoding: errors detected\n");
      }
}

void Decoder::print_sigma()
{
  std::cout<<std::endl<<"sigma(x) = ";
  for (int i = 0; i <= elp_degree[step_number]; i++)
    std::cout<<error_location_polynomial[step_number][i]<<" ";
  std::cout<<std::endl;
}

void Decoder::print_roots()
{
  std::cout<<"Roots: ";
  for (int i = 0; i < root_counter; ++i)
    std::cout<<loc[i]<<" ";
  std::cout<<std::endl;
}

void Decoder::int_correct_error_bits(int *cx_coefficients)
{
  if (root_counter == elp_degree[step_number])
    /* no. roots = degree of elp hence <= t errors */
    for (int i = 0; i < elp_degree[step_number]; i++)
      cx_coefficients[loc[i]] ^= 1;
  else
    throw 1;
}

int* Decoder::int_put_elp_into_index_form()
{
    int *index_of = gf_field->get_index_of();
    /* put elp into index form */
    for (int i = 0; i <= elp_degree[step_number]; i++)
        error_location_polynomial[step_number][i] = index_of[error_location_polynomial[step_number][i]];
    return index_of;

}

void Decoder::int_find_elp_roots(int* root, int* loc)
{
    int *alpha = gf_field->get_alpha();
    int  reg[201];
 
   /* Chien search: find roots of the error location polynomial */
    for (int i = 1; i <= elp_degree[step_number]; i++)
        reg[i] = error_location_polynomial[step_number][i];
    for (int i = 1; i <= gf_field->get_n(); i++) {
        q = 1;
        for (int j = 1; j <= elp_degree[step_number]; j++)
            if (reg[j] != -1) {
                reg[j] = (reg[j] + j) % gf_field->get_n();
                q ^= alpha[reg[j]];
            }
        if (!q) {	/* store root and error
                     * location number indices */
            root[root_counter] = i;
            loc[root_counter] = gf_field->get_n() - i;
            root_counter++;
        }
    }
}
