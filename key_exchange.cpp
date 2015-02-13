/*
key_echange.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#include "mmap/mmap.hpp"
#include "mmap/benchmark.hpp"

#include <iostream>
#include <getopt.h>

int main (int argc, char **argv)
{
	unsigned lambda = 52, kappa = 6;

	int c;
	opterr = 0;
	while ((c = getopt (argc, argv, "l:k:")) != -1)
	{
		switch (c)
		{
		case 'l':
			lambda = atoi(optarg);
			break;
		case 'k':
			kappa = atoi(optarg);
			break;
		case '?':
			std::cout << "Usage: " << argv[0] << " -l LAMBDA -k KAPPA" << std::endl;
			return 1;
		default:
			abort ();
		}
	}

	std::cout << std::endl;
	std::cout << "---------------------------------------------------------" << std::endl;
	std::cout << "Diffie-Hellman key exchange for N=" << kappa + 1 << " users and lambda=" << lambda << std::endl;

	/*
	* Load public parameters
	*/
	TIME( "Load pp",
	mmap::public_parameters pp(lambda, kappa);
	if (!pp.load())
	{
		std::cout << "ERROR: cannot load public parameters for ";
		std::cout << "lambda=" << lambda << " and kappa=" << kappa << "!" << std::endl;
		std::cout << "Did you generate them with ./generate_pp?" << std::endl;
		exit(0);
	}
	)

	/*
	* Publish
	*/
	mmap::encoding *c0s = new mmap::encoding[kappa + 1];
	mmap::encoding *c1s = new mmap::encoding[kappa + 1];

	// User private encodings
	auto publish_start = std::chrono::system_clock::now();
	TIME( "\tsamp",
	for (unsigned i = 0; i < kappa + 1; i++)
	{
		c0s[i].samp(pp, 0);
	}
	)

	// User public encodings
	TIME( "\tenc",
	for (unsigned i = 0; i < kappa + 1; i++)
	{
		c1s[i] = pp.enc(c0s[i], 1);
	}
	)

	TIME( "\trerand",
	for (unsigned i = 0; i < kappa + 1; i++)
	{
		c1s[i].rerand();
	}
	)

	auto publish_end = std::chrono::system_clock::now();
	std::cout << "Publish (per party): " << get_time_ms(publish_start, publish_end) / (kappa + 1) << "ms" << std::endl;

	/*
	* Keygen
	*/
	auto keygen_start = std::chrono::system_clock::now();
	// User multiplication
	mmap::encoding *ckappas = new mmap::encoding[kappa + 1];
	TIME("\tMult",
	for (unsigned i = 0; i < kappa + 1; i++)
	{
		ckappas[i] = c0s[i];
		for (unsigned j = 0; j < kappa + 1; j++)
		{
			if (j != i)
			{
				ckappas[i] *= c1s[j];
			}
		}
	}
	)

	// Extract
	mpz_class *keys = new mpz_class[kappa + 1];
	TIME("\tExtract",
	for (unsigned i = 0; i < kappa + 1; i++)
	{
		keys[i] = pp.ext(ckappas[i]);
	}
	)
	auto keygen_end = std::chrono::system_clock::now();
	std::cout << "Keygen (per party): " << get_time_ms(keygen_start, keygen_end) / (kappa + 1) << "ms" << std::endl;

	for (unsigned i = 0; i < kappa + 1; i++)
	{
		std::cout << "\t keys[" << i << "] = " << keys[i].get_str() << std::endl;
	}

	delete[] c0s;
	delete[] c1s;
	delete[] ckappas;
	delete[] keys;

	return 0;
}
