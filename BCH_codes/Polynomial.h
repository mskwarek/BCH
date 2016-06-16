//
// Created by marcin on 30.05.16.
//

#ifndef BCH_POLYNOMIAL_H
#define BCH_POLYNOMIAL_H

#include <stdio.h>

class Polynomial
{
private:
    int degree = 0;
    int code_length = 0;
    int polynomial_form[21] = {0};
    int n = 0; /* jw */
public:
    int* get_form();
    int get_degree();
    int get_code_length();
    int get_n();
    void read_p(int degree_ = 0);
    void print_primitive_polynomial();
};

#endif //BCH_POLYNOMIAL_H
