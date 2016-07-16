//
// Created by marcin on 01.05.16.
//
#include <gtest/gtest.h>
#include "../BCH_codes/Coder.h"
#include "test_cases.h"

class BCHCodingTest : public testing::Test
{
public:

    Coder *c;

    BCHCodingTest()
    {

      c = new Coder(10, 15);
    }

    virtual ~BCHCodingTest()
    {
        delete c;
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
    //    for (int i = 0; i < 1024; i++)
    //  EXPECT_EQ(data[i], recd[i]);
}

int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
