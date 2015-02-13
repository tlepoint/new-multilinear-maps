/*
mpz_shoup.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#ifndef MPZ_SHOUP
#define MPZ_SHOUP

void mpz_centered_mod(mpz_t rop, const mpz_t op, const mpz_t mod)
{
	mpz_mod(rop, op, mod);
	mpz_t modDiv;
	mpz_init(modDiv);
	mpz_tdiv_q_ui(modDiv, mod, 2);

	if (mpz_cmp(rop,modDiv)>0)
	{
		mpz_sub(rop, rop, mod);
	}

	mpz_clear(modDiv);
}

void mpz_compute_shoup(mpz_t rop, const mpz_t op, const mpz_t mod)
{
	mpz_mul_2exp(rop, op, mpz_sizeinbase(mod, 2) + 1);
	mpz_tdiv_q(rop, rop, mod);
}

void mpz_mulmod_shoup(mpz_t rop, const mpz_t op1, const mpz_t op2, const mpz_t op2shoup, const mpz_t mod)
{
	mpz_t tmp;
	mpz_init(tmp);

	mpz_mul(tmp, op1, op2shoup);
	mpz_tdiv_q_2exp(tmp, tmp, mpz_sizeinbase(mod, 2) + 1);
	mpz_mul(tmp, tmp, mod);

	mpz_mul(rop, op1, op2);
	mpz_sub(rop, rop, tmp);
	if (mpz_cmp(rop, mod) >= 0)
	{
		mpz_sub(rop, rop, mod);
	}

	mpz_clear(tmp);
}

#endif
