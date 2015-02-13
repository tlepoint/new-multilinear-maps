/*
mmap.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#ifndef MMAP_HPP
#define MMAP_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <chrono>

#include <sys/stat.h>
#include <sys/types.h>
#include <gmpxx.h>
#include <vector>
#include <math.h>

// Includes
#include "random.hpp"
#include "benchmark.hpp"
#include "crt_tree.hpp"
#include "mpz_shoup.hpp"

// http://stackoverflow.com/questions/3030099/pragma-in-define-macro
#define STRINGIFY(a) #a

/*
* What is implemented
*/
namespace mmap
{
struct parameters;
class public_parameters;
class encoding;
}

/*
* HEADERS
*/
namespace mmap
{

struct parameters
{
	unsigned lambda;
	unsigned kappa;
	unsigned n;
	unsigned eta;
	unsigned etaq;
	unsigned rho;
	unsigned alpha;
	unsigned xi;
	unsigned beta;
	unsigned nu;
	unsigned ne;
	unsigned rhof;
	std::string pathname;

	unsigned etap;
	unsigned nbp;
	unsigned delta;

	parameters() {}
	
	bool init(unsigned lambda, unsigned kappa);
	void create(unsigned n, unsigned rho, unsigned etap);

	bool load_mpz(mpz_class &, const char *filename);
	bool load_mpz_array(mpz_class *, unsigned, const char *filename);
	bool save_mpz(mpz_class &, const char *filename);
	bool save_mpz_array(mpz_class *, unsigned, const char *filename);
};

class public_parameters
{
protected:
	mpz_class *p;
	mpz_class *g;
	mpz_class z;
	mpz_class x0;
	mpz_class *CRT_coeffs_p;
	mpz_class *Qp; // Q mod p[i]
	mpz_class *h;
	mpz_class *Np; // N mod p_i
	crt_tree  crt;
	mpz_class **zki; // z^-k [kappa+1][n]
	bool params_loaded;

public:
	mpz_class *x;
	mpz_class **pi;
	mpz_class y;
	mpz_class x0p;
	mpz_class *yp;
	mpz_class N;
	mpz_class p_zt;
	parameters params;

	public_parameters(unsigned lambda, unsigned kappa);
	~public_parameters();

	bool load();

	void samp(mpz_class &rop, unsigned level);
	encoding enc(const encoding &c, unsigned level);
	void enc(mpz_class &rop, unsigned level);
	void rerand(mpz_class &rop, unsigned level);
	mpz_class ext(encoding &, unsigned nu = 0);

protected:
	void encrypt(mpz_class &, mpz_class *values, unsigned level = 0);

	std::vector<unsigned> decrypt(const mpz_class &, unsigned level);
	std::vector<unsigned> decrypt(const encoding &);
	std::vector<unsigned> noises(const mpz_class &, unsigned level);
	std::vector<unsigned> noises(const encoding &);
};

class encoding
{
public:
	public_parameters *pp;
	mpz_class value;
	unsigned level;

	encoding() {};
	encoding(const public_parameters &, unsigned level = 0);
	encoding(const encoding &c) : pp(c.pp), value(c.value), level(c.level) {};
	encoding(public_parameters *pp, const mpz_class &value, unsigned level): pp(pp), value(value), level(level) {};

	void samp(const public_parameters &, unsigned level);
	encoding enc(unsigned level);
	void rerand();

	encoding &operator*=(const encoding &);
	encoding &operator=(const encoding &);
};

}



/*
* IMPLEMENTATION
*/

namespace mmap
{

inline bool does_file_exist (const std::string &name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

inline unsigned nb_of_bits(unsigned x)
{
	unsigned nb(0), xp(x);
	while (xp > 0)
	{
		xp >>= 1;
		nb++;
	}
	return nb;
}

bool parameters::init(unsigned _lambda, unsigned _kappa)
{
	lambda = _lambda;
	kappa = _kappa;
	
	pathname = std::string("lambda") + std::to_string(lambda) + std::string("-kappa") + std::to_string(kappa);
	mkdir(pathname.c_str(), 0777);
	pathname += std::string("/");

	if (does_file_exist((pathname + std::string("parameters.mmap")).c_str()))
	{
		std::ifstream parametersFile;
		parametersFile.open((pathname + std::string("parameters.mmap")).c_str(), std::ios::binary);
		parametersFile >> n;
		parametersFile >> rhof;
		parametersFile >> eta;
		parametersFile >> etaq;
		parametersFile >> rho;
		parametersFile >> alpha;
		parametersFile >> xi;
		parametersFile >> beta;
		parametersFile >> nu;
		parametersFile >> ne;
		parametersFile.close();
		delta = floor(sqrt(n));
		return true;
	}

	return false;
}

bool parameters::load_mpz(mpz_class &a, const char *filename)
{
	if (does_file_exist((pathname + std::string(filename)).c_str()))
	{
		std::ifstream file;
		file.open((pathname + std::string(filename)).c_str(), std::ios::binary);
		std::string a_string;
		file >> a_string;
		a.set_str(a_string, 32);
		file.close();
		return true;
	}
	return false;
}

bool parameters::load_mpz_array(mpz_class *a, unsigned n, const char *filename)
{
	if (does_file_exist((pathname + std::string(filename)).c_str()))
	{
		std::ifstream file;
		file.open((pathname + std::string(filename)).c_str(), std::ios::binary);
		std::string a_string;
		for (unsigned i = 0; i < n; i++)
		{
			file >> a_string;
			a[i].set_str(a_string, 32);
		}
		file.close();
		return true;
	}
	return false;
}

public_parameters::public_parameters(unsigned lambda, unsigned kappa)
{
	params_loaded = params.init(lambda, kappa);
}

bool public_parameters::load()
{
	if (!params_loaded)
	{
		return false;
	}

	// Allocation
	pi = (mpz_class **) malloc(2 * sizeof(mpz_class *));
	pi[0] = new mpz_class[params.delta];
	pi[1] = new mpz_class[params.delta];
	x = new mpz_class[2 * params.lambda];
	yp = new mpz_class[params.ne];

	// Loads public parameters
	bool success = true;
	success &= params.load_mpz(x0p, "x0p.mmap");
	success &= params.load_mpz_array(x, 2 * params.lambda, "x.mmap");
	success &= params.load_mpz_array(pi[0], params.delta, "pi0.mmap");
	success &= params.load_mpz_array(pi[1], params.delta, "pi1.mmap");
	success &= params.load_mpz(y, "y.mmap");
	success &= params.load_mpz(N, "N.mmap");
	success &= params.load_mpz_array(yp, params.ne, "yp.mmap");
	success &= params.load_mpz(p_zt, "p_zt.mmap");

	// Load private parameters
	// p = new mpz_class[params.n];
	// g = new mpz_class[params.n];
	// success &= params.load_mpz_array(p, params.n, ".p.mmap");
	// success &= params.load_mpz_array(g, params.n, ".g.mmap");
	// success &= params.load_mpz(z, ".z.mmap");
	// success &= params.load_mpz(x0, ".x0.mmap");

	return success;
}

public_parameters::~public_parameters()
{
	if (p)
	{
		delete[] p;
	}
	if (g)
	{
		delete[] g;
	}
	if (CRT_coeffs_p)
	{
		delete[] CRT_coeffs_p;
	}
	if (Qp)
	{
		delete[] Qp;
	}
	if (x)
	{
		delete[] x;
	}
	if (pi)
	{
		delete[] pi[0];
		delete[] pi[1];
		free(pi);
	}
	if (Np)
	{
		delete[] Np;
	}
	if (h)
	{
		delete[] h;
	}
	if (yp)
	{
		delete[] yp;
	}
	if (zki)
	{
		for (unsigned k = 0; k < params.kappa + 1; k++)
			delete[] zki[k];
		free(zki);
	}
}

void public_parameters::samp(mpz_class &rop, unsigned level)
{
	assert(level == 0);
	rop = 0;
	for (unsigned i = 0; i < 2 * params.lambda; i++)
	{
		if (generate_random_bit() == 1)
		{
			rop += x[i];
		}
	}
}

void public_parameters::enc(mpz_class &rop, unsigned level)
{
	if (level == 0)
	{
		return;
	}

	for (unsigned i = 0; i < level; i++)
	{
		rop *= y;
		mpz_mod(rop.get_mpz_t(), rop.get_mpz_t(), x0p.get_mpz_t());
	}
}

void public_parameters::rerand(mpz_class &rop, unsigned level)
{
	assert(level == 1);
	mpz_class maskLeft = 0;
	mpz_class maskRight = 0;

	for (unsigned i = 0; i < params.delta; i++)
	{
		if (generate_random_bit() == 1)
		{
			maskLeft += pi[0][i];
		}
		if (generate_random_bit() == 1)
		{
			maskRight += pi[1][i];
		}
	}
	rop += maskLeft * maskRight;
	mpz_mod(rop.get_mpz_t(), rop.get_mpz_t(), x0p.get_mpz_t());
}

void public_parameters::encrypt(mpz_class &rop, mpz_class *values, unsigned level)
{
	if (level > 0)
	{
		for (unsigned j = 0; j < params.n; j++)
		{
			mpz_mul(values[j].get_mpz_t(), values[j].get_mpz_t(), zki[level][j].get_mpz_t());
		}
	}
	rop = crt.do_crt(values);
}

std::vector<unsigned> public_parameters::decrypt(const mpz_class &op, unsigned level)
{
	std::vector<unsigned> message(params.n);

	mpz_class value(op);
	for (unsigned j = 0; j < level; j++)
	{
		value *= z;
	}

	mpz_t coeff;
	mpz_init(coeff);
	for (unsigned i = 0; i < params.n; i++)
	{
		mpz_centered_mod(coeff, value.get_mpz_t(), p[i].get_mpz_t());
		mpz_mod(coeff, coeff, g[i].get_mpz_t());
		message[i] = mpz_get_ui(coeff);
	}
	mpz_clear(coeff);
	return message;
}

std::vector<unsigned> public_parameters::noises(const mpz_class &op, unsigned level)
{
	std::vector<unsigned> noises(params.n);

	mpz_class value(op);
	for (unsigned j = 0; j < level; j++)
	{
		value *= z;
	}

	mpz_t coeff;
	mpz_init(coeff);
	for (unsigned i = 0; i < params.n; i++)
	{
		mpz_centered_mod(coeff, value.get_mpz_t(), p[i].get_mpz_t());
		noises[i] = mpz_sizeinbase(coeff, 2);
	}
	mpz_clear(coeff);
	return noises;
}

std::vector<unsigned> public_parameters::decrypt(const encoding &c)
{
	return decrypt(c.value, c.level);
}

std::vector<unsigned> public_parameters::noises(const encoding &c)
{
	return noises(c.value, c.level);
}

mpz_class public_parameters::ext(encoding &c, unsigned nu)
{
	if (nu == 0)
	{
		nu = c.pp->params.nu;
	}
	mpz_class rop = c.value;
	for (unsigned j = 0; j < c.pp->params.ne; j++)
	{
		mpz_mod(rop.get_mpz_t(), rop.get_mpz_t(), c.pp->yp[c.pp->params.ne - 1 - j].get_mpz_t());
	}
	mpz_mul(rop.get_mpz_t(), rop.get_mpz_t(), p_zt.get_mpz_t());
	mpz_centered_mod(rop.get_mpz_t(), rop.get_mpz_t(), N.get_mpz_t());
	rop >>= (params.xi - nu);
	return rop;
}

encoding::encoding(const public_parameters &_pp, unsigned _level)
{
	samp(_pp, _level);
}

void encoding::samp(const public_parameters &_pp, unsigned _level)
{
	assert(_level == 0);
	level = _level;
	pp = const_cast<public_parameters *>(&_pp);
	pp->samp(value, level);
}

encoding public_parameters::enc(const encoding &c, unsigned _level)
{
	assert(c.level <= _level);
	encoding result(c);
	if (c.level != _level)
	{
		result.pp->enc(result.value, _level - result.level);
		result.level = _level;
	}
	return result;
}

void encoding::rerand()
{
	pp->rerand(value, level);
}

encoding &encoding::operator*=(const encoding &c)
{
	assert(pp == c.pp);
	mpz_mul(value.get_mpz_t(), value.get_mpz_t(), c.value.get_mpz_t());
	mpz_mod(value.get_mpz_t(), value.get_mpz_t(), pp->x0p.get_mpz_t());
	level += c.level;
	return *this;
}

encoding &encoding::operator=(const encoding &c)
{
	pp = c.pp;
	value = mpz_class(c.value);
	level = c.level;
	return *this;
}

}


#endif
