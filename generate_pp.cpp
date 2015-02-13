/*
generate_pp.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#include "mmap/mmap.hpp"
#include "mmap/mmap_gen.hpp"
#include "mmap/benchmark.hpp"

#include <iostream>

int main (int argc, char **argv)
{
	unsigned lambda = 52, kappa = 6, n = 540, rho = lambda, etap = 420;

	int c;
	opterr = 0;
	while ((c = getopt (argc, argv, "l:k:n:r:e:")) != -1)
	{
		switch (c)
		{
		case 'l':
			lambda = atoi(optarg);
			break;
		case 'k':
			kappa = atoi(optarg);
			break;
		case 'n':
			n = atoi(optarg);
			break;
		case 'r':
			rho = atoi(optarg);
			break;
		case 'e':
			etap = atoi(optarg);
			break;
		case '?':
			std::cout << "Usage: " << argv[0] << " -l LAMBDA -k KAPPA -n N -r RHO -e ETAP" << std::endl;
			std::cout << "\tETAP (default 420) is the size of the factors in the p_i's and N" << std::endl;
			return 1;
		default:
			abort ();
		}
	}

	// Public parameters
	std::cout << std::endl;
	std::cout << "---------------------------------------------------------" << std::endl;
	std::cout << "Generation of parameters for lambda=" << lambda
	<< ", kappa=" << kappa
	<< ", n=" << n
	<< ", rho=" << rho
	<< " and etap=" << etap << std::endl;

	// Initialization
	mmap::public_parameters_generate pp(lambda, kappa, n, rho, etap);

	// Generation of the parameters
	auto generationtime_start = std::chrono::system_clock::now();
	pp.generate();
	auto generationtime_end = std::chrono::system_clock::now();

	// Save the parameters on disk
	TIME( "Save parameters",
		pp.save();
	)


	std::cout << "Generation time: " << get_time_ms(generationtime_start, generationtime_end) << "ms" << std::endl;
	std::cout << "------------" << std::endl;
	std::cout << "eta=" << pp.params.eta << std::endl;
	std::cout << "rhof=" << pp.params.rhof << std::endl;
	std::cout << "ne=" << pp.params.ne << std::endl;
	std::cout << "------------" << std::endl;

	return 0;
}
