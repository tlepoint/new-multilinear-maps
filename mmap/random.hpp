/*
random.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <climits>
#include <gmpxx.h>
#include "prng/fastrandombytes.h"

// http://stackoverflow.com/questions/3030099/pragma-in-define-macro
#define STRINGIFY(a) #a


unsigned generate_random_bit();
void generate_random(mpz_class &, unsigned bits);
void generate_centered_random(mpz_class &, unsigned bits);
void generate_random(mpz_class &, const mpz_class &);
unsigned generate_random(unsigned n);


void generate_random(mpz_class &z, unsigned bits)
{
	if (bits == 0)
	{
		z = 0;
		return;
	}
	size_t nbchar = (bits) / (sizeof(unsigned char) * CHAR_BIT) + 1;
	unsigned char *rnd = new unsigned char[nbchar];
	mpz_class mask(2);
	mask <<= bits;
	mask -= 1;
	_Pragma(STRINGIFY(omp critical))
	{
		fastrandombytes(rnd, nbchar);
	}
	mpz_import (z.get_mpz_t(), nbchar, 1, sizeof(unsigned char), 0, 0, rnd);
	delete[] rnd;
	z &= mask;
}

void generate_centered_random(mpz_class &z, unsigned bits)
{
	generate_random(z, bits);
	if (generate_random_bit())
	{
		z *= -1;
	}
}

void generate_random(mpz_class &z, const mpz_class &n)
{
	if (n == 0)
	{
		z = 0;
		return;
	}
	size_t nbchar = mpz_sizeinbase(n.get_mpz_t(), 2) / (sizeof(unsigned char) * CHAR_BIT) + 1;
	unsigned char *rnd = new unsigned char[nbchar];
	mpz_class mask(2);
	mask <<= mpz_sizeinbase(n.get_mpz_t(), 2);
	mask -= 1;
	do
	{
		_Pragma(STRINGIFY(omp critical))
		{
			fastrandombytes(rnd, nbchar);
		}
		mpz_import (z.get_mpz_t(), nbchar, 1, sizeof(unsigned char), 0, 0, rnd);
		z &= mask;
	}
	while (z >= n);
	delete[] rnd;
}

unsigned generate_random_bit()
{
	unsigned char rnd[1];
	_Pragma(STRINGIFY(omp critical))
	{
		fastrandombytes(rnd, 1);
	}
	return rnd[0] & 1;
}

unsigned generate_random(unsigned n)
{
	if (n == 0)
	{
		return 0;
	}

	size_t nbchar = sizeof(unsigned) / sizeof(unsigned char);
	unsigned char *rnd = new unsigned char[nbchar];
	unsigned mask = 1;
	for (unsigned i = n; i > 0; i >>= 1)
	{
		mask <<= 1;
	}
	mask -= 1;

	unsigned result;
	do
	{
		_Pragma(STRINGIFY(omp critical))
		{
			fastrandombytes(rnd, nbchar);
		}
		result = 0;
		for (unsigned i = 0; i < nbchar; i++)
		{
			result <<= CHAR_BIT;
			result &= rnd[i];
		}
		result &= mask;
	}
	while (result >= n);

	delete[] rnd;
	return result;
}

#endif
