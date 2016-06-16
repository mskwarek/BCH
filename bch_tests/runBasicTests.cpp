//
// Created by marcin on 01.05.16.
//
#include <gtest/gtest.h>
#include "../BCH_codes/Coder.h"
#include "test_cases.h"

void int_rand_errors_position(Coder *c, int* errpos, int numerr);

class BCHCodingTest : public testing::Test
{
public:

    std::unique_ptr<Coder> c;

    BCHCodingTest()
    {
      c = std::unique_ptr<Coder>(new Coder(10, 15));
    }

    virtual ~BCHCodingTest()
    {}
};

class BCHCodingRandTest : public testing::Test
{
public:

    std::unique_ptr<Coder> c;
    int* recd;
    int bb[548576];
    int* errors;

    BCHCodingRandTest()
    {
        c = std::unique_ptr<Coder>(new Coder(10, 15));
        recd = new int[1048576];
        c->encode_bch(data, bb);
        c->code_polynomial(data, bb, recd);
        errors = new int[c->get_code_length()];
    }

    virtual ~BCHCodingRandTest()
    {
        delete recd;
        delete errors;
    }
};

TEST_F(BCHCodingTest, CheckPolynomialForm)
{
    for(int i = 0; i<11; i++)
        EXPECT_EQ(expected[i], c->get_polynomial_form()[i]);
}

TEST_F(BCHCodingTest, CheckGeneratedPolynomial)
{
    for(int i = 0; i<g_x.size(); i++)
        EXPECT_EQ(g_x[i], c->get_generated_polynomial()[i]);
}

TEST_F(BCHCodingTest, CheckBb)
{
    int bb[548576];
    c->encode_bch(data, bb);
    for(int i = 0; i< c->get_code_length(); i++)
        EXPECT_EQ(expected_bb[i], bb[i]);
}

TEST_F(BCHCodingTest, CheckAfterCoding)
{
    int recd[1048576];
    int bb[548576];
    c->encode_bch(data, bb);

    c->code_polynomial(data, bb, recd);
    for (int i = 0; i < c_x.size(); i++)
        EXPECT_EQ(c_x[i], recd[i]);
}

TEST_F(BCHCodingTest, Check)
{
    int* recd = new int[1048576];
    int bb[548576];
    c->encode_bch(data, bb);
    c->code_polynomial(data, bb, recd);

    for (int i = 0; i < 9; i++)
        recd[error_position[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_EQ(c->get_decoding_error_number(data, recd), 0);
}

TEST_F(BCHCodingRandTest, Rand15Errors)
{
    int_rand_errors_position(c.get(), errors, 15);
    for (int i = 0; i < 15; i++)
        recd[errors[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_EQ(c->get_decoding_error_number(data, recd), 0);
}

TEST_F(BCHCodingRandTest, Rand10Errors)
{
    int_rand_errors_position(c.get(), errors, 10);
    for (int i = 0; i < 10; i++)
        recd[errors[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_EQ(c->get_decoding_error_number(data, recd), 0);
}

TEST_F(BCHCodingRandTest, Rand5Errors)
{
    int_rand_errors_position(c.get(), errors, 10);
    for (int i = 0; i < 5; i++)
        recd[errors[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_EQ(c->get_decoding_error_number(data, recd), 0);
}

TEST_F(BCHCodingRandTest, RandNoErrors)
{
    int_rand_errors_position(c.get(), errors, 0);
    for (int i = 0; i < 0; i++)
        recd[errors[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_EQ(c->get_decoding_error_number(data, recd), 0);
}

TEST_F(BCHCodingRandTest, RandTooMuchErrors)
{
    int_rand_errors_position(c.get(), errors, 20);
    for (int i = 0; i < 20; i++)
        recd[errors[i]] ^= 1;
    c->decode_bch(recd);
    EXPECT_GT(c->get_decoding_error_number(data, recd), 0);
}

void int_rand_errors_position(Coder* c, int* errpos, int numerr)
{
    srand( time(NULL) );
    int tmp = 0;

    for (int i = 0; i < numerr; )
    {
        int flag = 0;
        tmp = rand() % c->get_code_length();
        for (int j = 0; j < numerr; ++j)
            if(errpos[j] == tmp)
                flag = 1;
        if(!flag) {
            errpos[i] = tmp;
            ++i;
        }
    }
}

int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
