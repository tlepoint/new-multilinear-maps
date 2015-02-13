# Proof-of-Concept Implementation of the "New" Multilinear Maps over the Integers

This is an implementation of the "new" variant of the cryptographic multilinear maps over the Integers.

This variant is described in the following article:

[1] J.-S. Coron, T. Lepoint, M. Tibouchi, "New Multilinear Maps over the Integers". Available on http://eprint.iacr.org

and build upon the following article:

[2] J.-S. Coron, T. Lepoint, M. Tibouchi, "Practical Multilinear Maps over the Integers". Available at http://eprint.iacr.org/2013/183

This C++ proof-of-concept requires:
- *GMP* http://gmplib.org/
- *fplll* (to generate parameters) https://github.com/dstehle/fplll

HOW TO RUN?
-----------

First, modify the Makefile accordingly.

To generate the public parameters (requires fplll)
--------------------------------------------------

```
$ make generate_pp
$ ./generate_pp -l 52 -k 6 -n 540 -r 52 -e 420
``` 

Explaination of the parameters: *-l* is the expected security level (lambda), *-k* is the multilinearity level (kappa), *-n* is the number of primes used, *-r* is the noise size (rho), *-e* is the size of the factors in the primes (etap). By default, running without argument yields the values of above.

This creates a directory *lambda52-kappa6* with all the public parameters of the multilinear map.

Note that all the arguments must be chosen to ensure sufficient security according to [1] and [2]. This is not performed by our implementation.

To run a kappa+1 Diffie-Hellman key exchange
--------------------------------------------

```
$ make key_generation
$ ./key_generation -l 52 -k 6
``` 

This requires the multilinear maps public parameters (created by __generate_pp__) to be stored in a *lambda52-kappa6/* directory. By default, running without argument yields the values of above.

Explaination of the parameters: *-l* is the expected security level (lambda), *-k* is the multilinearity level (kappa).

Note that there is a branch *with-small-parameters* containing small parameters to allow testing without running __generate_pp__.