//
// Created by marcin on 01.05.16.
//
#include <gtest/gtest.h>
#include "Coder.h"
#include "test_cases.h"

class BCHCodingTest : public testing::Test
{
public:
    GfField *gf;
    Coder *c;

    BCHCodingTest()
    {
        gf = new GfField(10);
        c = new Coder();
        gf->generate_gf();
        gf->gen_poly(15);
    }

    virtual ~BCHCodingTest()
    {
        delete gf;
        delete c;
    }
};

TEST_F(BCHCodingTest, CheckPolynomialForm)
{
    for(int i = 0; i<11; i++)
        EXPECT_EQ(expected[i], gf->get_polynomial_form()[i]);
}

TEST_F(BCHCodingTest, CheckGeneratedPolynomial)
{
    for(int i = 0; i<g_x.size(); i++)
        EXPECT_EQ(g_x[i], gf->get_generated_poly()[i]);
}

TEST_F(BCHCodingTest, CheckBb)
{
    int bb[548576];
    c->encode_bch(gf, data, bb);
    for(int i = 0; i<gf->get_code_length(); i++)
        EXPECT_EQ(expected_bb[i], bb[i]);
}

TEST_F(BCHCodingTest, CheckAfterCoding)
{
    int recd[1048576];
    int bb[548576];
    c->encode_bch(gf, data, bb);

    c->code_polynomial(gf->get_code_length(), gf->get_k(), data, bb, recd);
    for (int i = 0; i < c_x.size(); i++)
        EXPECT_EQ(c_x[i], recd[i]);
}

TEST_F(BCHCodingTest, Check)
{
    int recd[1048576];
    int bb[548576];
    c->encode_bch(gf, data, bb);
    c->code_polynomial(gf->get_code_length(), gf->get_k(), data, bb, recd);

    for (int i = 0; i < 9; i++)
        recd[error_position[i]] ^= 1;

    //c->decode_bch(gf, recd);
}

int main(int ac, char* av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}