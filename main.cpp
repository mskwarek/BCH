/*
 * File:    bch3.c
 * Title:   Encoder/decoder for binary BCH codes in C (Version 3.1)
 * Author:  Robert Morelos-Zaragoza
 * Date:    August 1994
 * Revised: June 13, 1997
 *
 * ===============  Encoder/Decoder for binary BCH codes in C =================
 *
 * Version 1:   Original program. The user provides the generator polynomial
 *              of the code (cumbersome!).
 * Version 2:   Computes the generator polynomial of the code.
 * Version 3:   No need to input the coefficients of a primitive polynomial of
 *              degree m, used to construct the Galois Field GF(2**m). The
 *              program now works for any binary BCH code of length such that:
 *              2**(m-1) - 1 < length <= 2**m - 1
 *
 * Note:        You may have to change the size of the arrays to make it work.
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983.
 *
 * Thanks to Patrick Boyle (pboyle@era.com) for his observation that 'bch2.c'
 * did not work for lengths other than 2**m-1 which led to this new version.
 * Portions of this program are from 'rs.c', a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89. The
 * previous version of the BCH encoder/decoder in C, 'bch2.c', was written by
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) on 5/19/92.
 *
 * NOTE:
 *          The author is not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 *  For more information, suggestions, or other ideas on implementing error
 *  correcting codes, please contact me at:
 *
 *                           Robert Morelos-Zaragoza
 *                           5120 Woodway, Suite 7036
 *                           Houston, Texas 77056
 *
 *                    email: r.morelos-zaragoza@ieee.org
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * == Copyright (c) 1994-7,  Robert Morelos-Zaragoza. All rights reserved.  ==
 *
 * m = order of the Galois field GF(2**m)
 * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
 * length = length of the BCH code
 * t = error correcting capability (max. no. of errors the code corrects)
 * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
 * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
 * p[] = coefficients of a primitive polynomial used to generate GF(2**m)
 * g[] = coefficients of the generator polynomial, g(x)
 * alpha_to [] = log table of GF(2**m)
 * index_of[] = antilog table of GF(2**m)
 * data[] = information bits = coefficients of data polynomial, i(x)
 * bb[] = coefficients of redundancy polynomial x^(length-k) i(x) modulo g(x)
 * numerr = number of errors
 * errpos[] = error positions
 * recd[] = coefficients of the received polynomial
 * decerror = number of decoding errors (in _message_ positions)
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <map>

#include "BCH_codes/Coder.h"
#include "BCH_codes/GfField.h"
#include "BCH_codes/Config.h"

int             i;
int             data[1048576];
int             numerr, errpos[1024];

void int_rand_errors(Coder *c, int err = 0);
void int_rand_errors_position(Coder *c);
void int_print_polynomial(std::string title, std::string polynomial_symbol, int *data, int start, int length);

void get_random_data(int k, int *data_destination)
{
	int             seed;
	/* Randomly generate DATA */
	seed = 131073;
	srandom(seed);
	for (int i = 0; i < k; i++)
		data_destination[i] = ( random() & 65536 ) >> 16;
}


int main()
{
        Config cnf;
	cnf.load("config.xml");
        std::map<std::string, std::string> codes_features = cnf.getDatabaseInfo();
	int* recd = new int [1048576];
	int bb[548576];
	int poly_degree = std::stoi(codes_features["polynomial_degree"]);
	int correction_cap = std::stoi(codes_features["correction_capabilities"]);
	Coder* c = new Coder(poly_degree, correction_cap);
	get_random_data(c->get_k(), data);
	c->encode_bch(data, bb);
	c->code_polynomial(data, bb, recd);
	
	int_print_polynomial("Code polynomial:", "c(x)", recd, 0, c->get_code_length());

	int_rand_errors(c, 9);
	if(9)
	{
		for (i = 0; i < numerr; i++)
			recd[errpos[i]] ^= 1;

		int_print_polynomial("r(x)", "r(x)", recd, 0, c->get_code_length());
	}
	c->decode_bch(recd);


	printf("Results:\n");
	int_print_polynomial("original data:", "send(x)", data, 0, c->get_k());
	int_print_polynomial("recovered data:", "recovered(x)", recd, c->get_code_length() - c->get_k(), c->get_code_length());

	c->get_decoding_error_number(data, recd);

	delete c;
	return 0;
}



void int_rand_errors(Coder *c, int err)
{
	numerr = err;
	if(err==0) {
		printf("Enter the number of errors:\n");
		scanf("%d", &numerr);    /* CHANNEL errors */

		printf("Enter error locations (integers between");
		printf(" 0 and %d): ", c->get_code_length() - 1);
		/*
		 * recd[] are the coefficients of r(x) = c(x) + e(x)
		 */
		for (i = 0; i < numerr; i++)
			scanf("%d", &errpos[i]);
	}
	else
		int_rand_errors_position(c);
}

void int_rand_errors_position(Coder *c)
{
	srand( time(NULL));
	int tmp = 0;
	std::cout<<"ERRORS at: ";
	for (i = 0; i < numerr; )
	{
		int flag = 0;
		tmp = rand() % c->get_code_length();
		for (int j = 0; j < numerr; ++j)
			if(errpos[j] == tmp)
				flag = 1;
		if(!flag) {
			std::cout << tmp << ", ";
			errpos[i] = tmp;
			++i;
		}
	}
	std::cout<<std::endl;
}

void int_print_polynomial(std::string title, std::string polynomial_symbol, int *data, int start, int length)
{
	std::cout<<title<<std::endl;
	std::cout<<polynomial_symbol<<" = ";
	for (i = start; i < length; i++) {
		printf("%1d", data[i]);
		if (i && ((i % 50) == 0))
			printf("\n");
	}
	printf("\n");
}
