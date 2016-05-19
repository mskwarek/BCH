//
// Created by marcin on 01.05.16.
//

#include <stdio.h>

class Polynomial
{
private:
    int degree = 0;
    int code_length = 0;
    int polynomial_form[21] = {0};

    int n = 0; /* jw */
public:
    int* get_form()
    {
        return polynomial_form;
    }

    int get_degree()
    {
        return this->degree;
    }

    int get_code_length()
    {
        return code_length;
    }

    int get_n()
    {
        return n;
    }

    void read_p(int degree_ = 0)
/*
 *	Read m, the degree of a primitive polynomial p(x) used to compute the
 *	Galois field GF(2**m). Get precomputed coefficients p[] of p(x). Read
 *	the code length.
 */
    {
        int			i, ninf;

//        printf("bch3: An encoder/decoder for binary BCH codes\n");
//        printf("Copyright (c) 1994-7. Robert Morelos-Zaragoza.\n");
//        printf("This program is free, please read first the copyright notice.\n");
//        printf("\nFirst, enter a value of m such that the code length is\n");
//        printf("2**(m-1) - 1 < length <= 2**m - 1\n\n");
        degree = degree_;
        do {
            if(degree_ == 0) {
                printf("Enter m (between 2 and 20): ");
                scanf("%d", &degree);
            }
        } while ( !(degree>1) || !(degree<21) );
        for (i=1; i<degree; i++)
            polynomial_form[i] = 0;
        polynomial_form[0] = polynomial_form[degree] = 1;
        if (degree == 2)			polynomial_form[1] = 1;
        else if (degree == 3)	polynomial_form[1] = 1;
        else if (degree == 4)	polynomial_form[1] = 1;
        else if (degree == 5)	polynomial_form[2] = 1;
        else if (degree == 6)	polynomial_form[1] = 1;
        else if (degree == 7)	polynomial_form[1] = 1;
        else if (degree == 8)	polynomial_form[4] = polynomial_form[5] = polynomial_form[6] = 1;
        else if (degree == 9)	polynomial_form[4] = 1;
        else if (degree == 10)	polynomial_form[3] = 1;
        else if (degree == 11)	polynomial_form[2] = 1;
        else if (degree == 12)	polynomial_form[3] = polynomial_form[4] = polynomial_form[7] = 1;
        else if (degree == 13)	polynomial_form[1] = polynomial_form[3] = polynomial_form[4] = 1;
        else if (degree == 14)	polynomial_form[1] = polynomial_form[11] = polynomial_form[12] = 1;
        else if (degree == 15)	polynomial_form[1] = 1;
        else if (degree == 16)	polynomial_form[2] = polynomial_form[3] = polynomial_form[5] = 1;
        else if (degree == 17)	polynomial_form[3] = 1;
        else if (degree == 18)	polynomial_form[7] = 1;
        else if (degree == 19)	polynomial_form[1] = polynomial_form[5] = polynomial_form[6] = 1;
        else if (degree == 20)	polynomial_form[3] = 1;
        printf("p(x) = ");
        n = 1;
        for (i = 0; i <= degree; i++) {
            n *= 2;
            printf("%1d", polynomial_form[i]);
        }
        printf("\n");
        n = n / 2 - 1;
        ninf = (n + 1) / 2 - 1;
        /*
        do  {
            printf("Enter code length (%d < length <= %d): ", ninf, n);
            scanf("%d", &code_length);
        } while ( !((code_length <= n)&&(code_length>ninf)) );
         */
        code_length = (n - ninf)/2 + ninf;
    }



};