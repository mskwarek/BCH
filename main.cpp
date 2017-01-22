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

int i;
int data[1048576];
int numerr, errpos[1024];

void int_rand_errors(Coder *c, int err = 0);
void int_rand_errors_position(Coder *c);
void int_print_polynomial(std::string title, std::string polynomial_symbol, int *data, int start, int length);

void get_random_data(int k, int *data_destination)
{
	int seed = 131073;
	srandom(seed);
	for (int i = 0; i < k; i++)
	  {
		data_destination[i] = ( random() & 65536 ) >> 16;

	  }
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
	int error_num = std::stoi(codes_features["number_of_errors"]);
	Coder* c = new Coder(poly_degree, correction_cap);
	get_random_data(c->get_k(), data);
	c->encode_bch(data, bb);
	c->code_polynomial(data, bb, recd);
	
	int_print_polynomial("Code polynomial:", "c(x)", recd, 0, c->get_code_length());

	int_rand_errors(c, error_num);
	if(error_num)
	{
		for (i = 0; i < numerr; i++)
			recd[errpos[i]] ^= 1;

		int_print_polynomial("r(x)", "r(x)", recd, 0, c->get_code_length());
	}
	c->decode_bch(recd);


	std::cout<<"Results:\n";
	int_print_polynomial("original data:", "send(x)", data, 0, c->get_k());
	int_print_polynomial("recovered data:", "recovered(x)", recd, c->get_code_length() - c->get_k(), c->get_code_length());

	c->get_decoding_error_number(data, recd);

	delete c;
	return 0;
}

void int_rand_errors(Coder *c, int err)
{
	numerr = err;
	if(err==0)
	{
		printf("Enter the number of errors:\n");
		scanf("%d", &numerr);    /* CHANNEL errors */

		printf("Enter error locations (integers between");
		printf(" 0 and %d): ", c->get_code_length() - 1);
		for (i = 0; i < numerr; i++)
		{
			scanf("%d", &errpos[i]);
		}
	}
	else
	{
		int_rand_errors_position(c);
	}
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
		if(!flag)
		{
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
		  std::cout<<std::endl;
	}
	std::cout<<std::endl;
}
