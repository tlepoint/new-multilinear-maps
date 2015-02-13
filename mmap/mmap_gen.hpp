/*
mmap_gen.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#ifndef MMAP_GEN_HPP
#define MMAP_GEN_HPP

#include "mmap.hpp"

// Fplll
#include <fplll.h>

// http://stackoverflow.com/questions/3030099/pragma-in-define-macro
#define STRINGIFY(a) #a

namespace mmap {

class public_parameters_generate : public public_parameters
{
private:
	mpz_class z_invert;
	mpz_class z_mkappa;
	mpz_class *zi;
	mpz_class *ziinv; // z^(-1) mod p_i
	mpz_class x0inv;
	mpz_class x0inv_shoup;
	std::vector<fplll::ZZ_mat<mpz_t>> Mat;
	crt_tree  crtN;

public:
	public_parameters_generate(unsigned lambda, unsigned kappa, unsigned n, unsigned rho, unsigned etap);
	~public_parameters_generate();
	void generate();
	void save();

private:
	void generate_p(unsigned i);
	void generate_g(unsigned i);
	void generate_z();
	void generate_CRT_coeffs(unsigned i);
	void generate_x();
	void generate_pi(unsigned level);
	void generate_y();
	void generate_yp();
	void generate_N();
	void generate_h();
	void generate_p_zt();
	void generate_p_zt_coeff(unsigned i);

	void multisamp(mpz_class **results, mpz_class ***values, unsigned m, unsigned level);
};

}

/*
* IMPLEMENTATION
*/

namespace mmap
{

void parameters::create(unsigned _n, unsigned _rho, unsigned _etap)
{
	n = _n;
	rho = _rho;
	etap = _etap;

	alpha 	= lambda;
	beta 	= lambda;
	nu 		= lambda;

	// Computations in the paper
	rhof 	= kappa * (2 * rho + 2 * alpha + nb_of_bits(2 * lambda)) + rho + nb_of_bits(2 * lambda) + 10;
	eta 	= rhof + 2 * alpha + 2 * beta + lambda + 8 + nu;
	xi 		= n * eta + 2 * eta + 1;
	etaq 	= 2 * eta + alpha;
	ne 		= (etaq - rhof) / (rhof - rho) + 1;

	std::ofstream parametersFile;
	parametersFile.open((pathname + std::string("parameters.mmap")).c_str(), std::ios::binary);
	parametersFile << n << std::endl;
	parametersFile << rhof << std::endl;
	parametersFile << eta << std::endl;
	parametersFile << etaq << std::endl;
	parametersFile << rho << std::endl;
	parametersFile << alpha << std::endl;
	parametersFile << xi << std::endl;
	parametersFile << beta << std::endl;
	parametersFile << nu << std::endl;
	parametersFile << ne << std::endl;
	parametersFile.close();

	std::cout << "(eta = " << eta << ")" << std::endl;
	if (eta%etap < 350)
	{
		std::cout << "WARNING: Your value of 'etap' might yield an insecure system as the p_i's will have factors of size eta\%etap=" << eta%etap << std::endl; 
		std::cout << "\t\t(Also the assertion is crt_tree might fail because some factors can be equal.)" << std::endl;
	}

	nbp = eta / etap + ((eta % etap == 0) ? 0 : 1);
	delta = floor(sqrt(n));
}

public_parameters_generate::public_parameters_generate(unsigned lambda, unsigned kappa, unsigned n, unsigned rho, unsigned etap) : public_parameters(lambda, kappa)
{
	params.create(n, rho, etap);
}

public_parameters_generate::~public_parameters_generate()
{
	if (zi)
	{
		delete[] zi;
	}
	if (ziinv)
	{
		delete[] ziinv;
	}
}

bool parameters::save_mpz(mpz_class &a, const char *filename)
{
	std::ofstream file;
	file.open((pathname + std::string(filename)).c_str(), std::ios::binary | std::ios::trunc);
	file << a.get_str(32) << std::endl;
	file.close();
	return true;
}

bool parameters::save_mpz_array(mpz_class *a, unsigned n, const char *filename)
{
	std::ofstream file;
	file.open((pathname + std::string(filename)).c_str(), std::ios::binary | std::ios::trunc);
	for (unsigned i = 0; i < n; i++)
	{
		file << a[i].get_str(32) << std::endl;
	}
	file.close();
	return true;
}

void public_parameters_generate::generate()
{
	TIME("\tPrimes p",
	// Primes p
	p = new mpz_class[params.n];
	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned i = 0; i < params.n; i++)
	{
		generate_p(i);
	}
	)

	TIME("\tx0, x0p, crt (new)",
	// x0 + crt
	crt.init(p, params.n);
	x0 = crt.get_mod();

	//x0p
	x0p = x0;
	mpz_class q;
	generate_random(q, params.etaq);
	x0p *= q;
	)

	TIME("\tg",
	// Primes g
	g = new mpz_class[params.n];
	for (unsigned i = 0; i < params.n; i++)
	{
		generate_g(i);
	}
	)

	TIME("\tz, z_invert (new)",
	// Mask z & z_invert
	zi = new mpz_class[params.n];
	for (unsigned j = 0; j < params.n; j++)
	{
		generate_random(zi[j], p[j]);
	}
	z = crt.do_crt(zi);

	ziinv = new mpz_class[params.n];
	for (unsigned j = 0; j < params.n; j++)
	{
		mpz_invert(ziinv[j].get_mpz_t(), zi[j].get_mpz_t(), p[j].get_mpz_t());
	}
	z_invert = crt.do_crt(ziinv);

	zki = (mpz_class **) malloc((params.kappa + 1) * sizeof(mpz_class *));
	for (unsigned k = 0; k < params.kappa + 1; k++)
	{
		zki[k] = new mpz_class[params.n];
		for (unsigned j = 0; j < params.n; j++)
		{
			mpz_powm_ui(zki[k][j].get_mpz_t(), ziinv[j].get_mpz_t(),
			            k, p[j].get_mpz_t());
		}
	}

	z_mkappa = crt.do_crt(zki[params.kappa]);
	)

	TIME("\tCRT_coeffs_p and Qp",
	// CRT Coefficients & Shoup coefficients
	CRT_coeffs_p = new mpz_class[params.n];
	Qp = new mpz_class[2 * params.n];
	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned i = 0; i < params.n; i++)
	{
		generate_CRT_coeffs(i);
	}
	)

	TIME("\tN",
	// N
	Np = new mpz_class[2 * params.n];
	generate_N();
	)

	TIME("\th",
	// h
	h = new mpz_class[params.n];
	generate_h();
	)

	TIME("\tGeneration yp, p_zt",
	// yp = enc(0, kappa)
	yp = new mpz_class[params.ne];
	generate_yp();

	// Matrices for LLL
	Mat.assign(params.n, ZZ_mat<mpz_t>(2, 2));
	mpz_class one = 1;
	mpz_class m(N);
	m >>= (params.xi - (2 * (params.eta - 3)));
	mpz_realloc2(m.get_mpz_t(), 2 * (params.eta - 3));
	for (unsigned i = 0; i < params.n; i++)
	{
		mpz_set(Mat[i](0, 0).getData(), one.get_mpz_t());
		mpz_set(Mat[i](1, 1).getData(), m.get_mpz_t());
	}

	// p_zt
	generate_p_zt();
	)

	TIME("\tSampling x's",
	// x's for sampling
	x = new mpz_class[2 * params.lambda];
	generate_x();
	// )
	// TIME("\tRerand pi",
	// pi's for rerandomization
	pi = (mpz_class **) malloc(2 * sizeof(mpz_class *));
	pi[0] = new mpz_class[params.delta];
	pi[1] = new mpz_class[params.delta];
	generate_pi(0);
	generate_pi(1);
	)
	
	TIME("\ty",
	// y = enc(1, 1)
	generate_y();
	)

}

void public_parameters_generate::save()
{
	// Private
	// params.save_mpz_array(p, params.n, ".p.mmap");
	// params.save_mpz(x0, ".x0.mmap");
	// params.save_mpz_array(g, params.n, ".g.mmap");
	// params.save_mpz(z, ".z.mmap");
	// params.save_mpz_array(CRT_coeffs_p, params.n, ".CRT_coeffs_p.mmap");
	// params.save_mpz_array(Qp, 2 * params.n, ".Qp.mmap");
	// params.save_mpz_array(Np, 2 * params.n, ".Np.mmap");
	// params.save_mpz_array(h, params.n, ".h.mmap");

	// Public
	params.save_mpz(x0p, "x0p.mmap");
	params.save_mpz_array(x, 2 * params.lambda, "x.mmap");
	params.save_mpz_array(pi[0], params.delta, "pi0.mmap");
	params.save_mpz_array(pi[1], params.delta, "pi1.mmap");
	params.save_mpz(y, "y.mmap");
	params.save_mpz(N, "N.mmap");
	params.save_mpz_array(yp, params.ne, "yp.mmap");
	params.save_mpz(p_zt, "p_zt.mmap");
}

void generate_prod_prime(mpz_class * array, unsigned i, unsigned size, unsigned chunk_size)
{
	unsigned num = size / chunk_size + ((size%chunk_size == 0) ? 0 : 1);
	mpz_class val;
	array[i] = 1;
	for (unsigned j = 0; j < num; j++)
	{
		generate_random(val, (j < num - 1) ? chunk_size : size - (num - 1)*chunk_size);
		mpz_nextprime(val.get_mpz_t(), val.get_mpz_t());
		array[i] *= val;
	}
}

void public_parameters_generate::generate_p(unsigned i)
{
	// mpz_class val;
	// p[i] = 1;
	// for (unsigned j = 0; j < params.nbp; j++)
	// {
	// 	generate_random(val, (j < params.nbp - 1) ? params.etap : params.eta - (params.nbp - 1)*params.etap);
	// 	mpz_nextprime(val.get_mpz_t(), val.get_mpz_t());
	// 	p[i] *= val;
	// }
	generate_prod_prime(p, i, params.eta, params.etap);
}

void public_parameters_generate::generate_g(unsigned i)
{
	generate_random(g[i], params.alpha);
	mpz_nextprime(g[i].get_mpz_t(), g[i].get_mpz_t());
}

void public_parameters_generate::generate_z()
{
	generate_random(z, x0);
}

void public_parameters_generate::generate_CRT_coeffs(unsigned i)
{
	mpz_divexact(Qp[i].get_mpz_t(), x0.get_mpz_t(), p[i].get_mpz_t()); // Q = x0/p[i]
	mpz_mod(Qp[i].get_mpz_t(), Qp[i].get_mpz_t(), p[i].get_mpz_t());
	mpz_compute_shoup(Qp[i + params.n].get_mpz_t(), Qp[i].get_mpz_t(), p[i].get_mpz_t());

	mpz_realloc2(Qp[i].get_mpz_t(), params.eta + 1);
	mpz_realloc2(Qp[i + params.n].get_mpz_t(), params.eta + 2);

	mpz_invert(CRT_coeffs_p[i].get_mpz_t(), Qp[i].get_mpz_t(), p[i].get_mpz_t());
	mpz_mul(CRT_coeffs_p[i].get_mpz_t(), CRT_coeffs_p[i].get_mpz_t(), z_mkappa.get_mpz_t());
	mpz_mul(CRT_coeffs_p[i].get_mpz_t(), CRT_coeffs_p[i].get_mpz_t(), g[i].get_mpz_t());
	mpz_mod(CRT_coeffs_p[i].get_mpz_t(), CRT_coeffs_p[i].get_mpz_t(), p[i].get_mpz_t());

	mpz_realloc2(CRT_coeffs_p[i].get_mpz_t(), params.eta + 1);
}

void public_parameters_generate::generate_x()
{
	mpz_class **values = (mpz_class **) malloc(2 * params.lambda * sizeof(mpz_class *));
	for (unsigned j = 0; j < 2 * params.lambda; j++)
	{
		values[j] = new mpz_class[params.n];
		for (unsigned i = 0; i < params.n; i++)
		{
			generate_centered_random(values[j][i], params.alpha + params.rho);
			values[j][i] += (values[j][i] < 0) ? p[i] : 0;
		}
	}
	multisamp(&x, &values, 2 * params.lambda, 0);

	for (unsigned j = 0; j < 2 * params.lambda; j++)
		delete[] values[j];
	free(values);
}

void public_parameters_generate::generate_pi(unsigned level)
{
	mpz_class **values = (mpz_class **) malloc(params.delta * sizeof(mpz_class *));
	for (unsigned j = 0; j < params.delta; j++)
	{
		values[j] = new mpz_class[params.n];
		for (unsigned i = 0; i < params.n; i++)
		{
			generate_centered_random(values[j][i], params.rho);
			values[j][i] *= g[i];
			values[j][i] += (values[j][i] < 0) ? p[i] : 0;
		}
	}
	multisamp(&pi[level], &values, params.delta, level);

	for (unsigned j = 0; j < params.delta; j++)
		delete[] values[j];
	free(values);
}

void public_parameters_generate::generate_y()
{
	mpz_class *values = new mpz_class[params.n];
	for (unsigned int i = 0; i < params.n; i++)
	{
		generate_centered_random(values[i], params.rho);
		values[i] *= g[i];
		values[i] += 1;
		values[i] += (values[i] < 0) ? p[i] : 0;
	}
	encrypt(y, values, 1);
	delete[] values;
}

void public_parameters_generate::generate_p_zt_coeff(unsigned i)
{
	// Q = x0/p_i
	mpz_class Q;
	mpz_divexact(Q.get_mpz_t(), x0.get_mpz_t(), p[i].get_mpz_t()); // Q = x0/p[i]

	// u = x_0/p_i * [z^(-kappa)*g_i*(x_0/p_i)]_{p_i}
	mpz_class u = CRT_coeffs_p[i] * Q;

	// The following is apparently less efficient than recomputing Q,
	// etc.
	//
	// mpz_class u = crt.do_crt_oneindex(g[i] * zki[params.kappa][i], i);

	//
	// M=Matrix(ZZ,[[1, MSB(pinv[i]*u[i]]),[0,MSB(N)]])
	//
	mpz_class tmp;
	// tmp = CRT_coeffs_p[i]*Q mod p[i]
	mpz_mulmod_shoup(tmp.get_mpz_t(), CRT_coeffs_p[i].get_mpz_t(), Qp[i].get_mpz_t(), Qp[i + params.n].get_mpz_t(), p[i].get_mpz_t());
	// tmp = tmp * N^(-1) mod p[i]
	mpz_mulmod_shoup(tmp.get_mpz_t(), tmp.get_mpz_t(), Np[i].get_mpz_t(), Np[i + params.n].get_mpz_t(), p[i].get_mpz_t());
	// tmp *= N
	mpz_mul(tmp.get_mpz_t(), tmp.get_mpz_t(), N.get_mpz_t());
	// (u - tmp)/p
	mpz_sub(u.get_mpz_t(), u.get_mpz_t(), tmp.get_mpz_t());
	mpz_divexact(u.get_mpz_t(), u.get_mpz_t(), p[i].get_mpz_t());
	mpz_centered_mod(u.get_mpz_t(), u.get_mpz_t(), N.get_mpz_t());
	u >>= (params.xi - (2 * (params.eta - 3)));
	mpz_realloc2(u.get_mpz_t(), 2 * (params.eta - 3));

	mpz_set(Mat[i](0, 1).getData(), u.get_mpz_t());
	fplll::lllReduction(Mat[i]);

	_Pragma(STRINGIFY(omp critical))
	{
		p_zt += h[i] * Q * mpz_class(Mat[i](0, 0).getData());
	}

	Mat[i].clear();
}

void public_parameters_generate::generate_p_zt()
{
	p_zt = 0;

	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned i = 0; i < params.n; i++)
	{
		generate_p_zt_coeff(i);
	}

	mpz_mod(p_zt.get_mpz_t(), p_zt.get_mpz_t(), N.get_mpz_t());
	mpz_mulmod_shoup(p_zt.get_mpz_t(), p_zt.get_mpz_t(), x0inv.get_mpz_t(), x0inv_shoup.get_mpz_t(), N.get_mpz_t());
}

void public_parameters_generate::generate_yp()
{
	mpz_class *values = new mpz_class[params.n];
	mpz_class q;

	for (unsigned j = 0 ; j < params.ne; j++)
	{
		for (unsigned int i = 0; i < params.n; i++)
		{
			generate_centered_random(values[i], params.rho);
			values[i] *= g[i];
			values[i] += (values[i] < 0) ? p[i] : 0;
		}
		encrypt(yp[j], values, params.kappa);
		generate_random(q, params.rhof + j * (params.rhof - params.rho));
		yp[j] += q * x0;
	}
	delete[] values;
}

void public_parameters_generate::multisamp(mpz_class **results, mpz_class ***values, unsigned m, unsigned level)
{
	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned j = 0; j < m; j++)
	{
		encrypt((*results)[j], (*values)[j], level);
	}
}

void public_parameters_generate::generate_N()
{
	mpz_class *val = new mpz_class[params.n];
	unsigned size_factor = params.xi / params.n + ((params.xi % params.n == 0) ? 0 : 1);
	
	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned i = 0; i<params.n-1; i++)
	{
		generate_prod_prime(val, i, size_factor, params.etap);
	}
	generate_prod_prime(val, params.n-1, params.xi-(params.n-1)*size_factor, params.etap);
	
	N = 1;
	for (unsigned i = 0; i<params.n; i++)
	{
		N *= val[i];
	}

	// Correction of missing bits in N
	if (mpz_sizeinbase(N.get_mpz_t(), 2) < params.xi)
	{
		generate_prod_prime(val, 0, params.xi-mpz_sizeinbase(N.get_mpz_t(), 2), params.etap);
		N *= val[0];
	}

	delete[] val;

	mpz_invert(x0inv.get_mpz_t(), x0.get_mpz_t(), N.get_mpz_t());
	mpz_compute_shoup(x0inv_shoup.get_mpz_t(), x0inv.get_mpz_t(), N.get_mpz_t());

	_Pragma(STRINGIFY(omp parallel for))
	for (unsigned i = 0; i < params.n; i++)
	{
		mpz_invert(Np[i].get_mpz_t(), N.get_mpz_t(), p[i].get_mpz_t());
		mpz_compute_shoup(Np[i + params.n].get_mpz_t(), Np[i].get_mpz_t(), p[i].get_mpz_t());
	}
}

void public_parameters_generate::generate_h()
{
	for (unsigned i = 0; i < params.n; i++)
	{
		generate_random(h[i], params.beta);
	}
}

}


#endif
