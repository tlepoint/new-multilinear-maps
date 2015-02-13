/*
crt_tree.hpp version 20150213
Mehdi Tibouchi
Public domain.
*/

#ifndef CRT_TREE_HPP
#define CRT_TREE_HPP

class crt_tree
{
private:
	unsigned  n, n2;
	mpz_class mod;
	mpz_class crt_left;
	mpz_class crt_right;
	crt_tree* left;
	crt_tree* right;

public:
	crt_tree() {
		left  = NULL;
		right = NULL;
		n  = 0;
		n2 = 0;
	}

	crt_tree(const mpz_class* ps, unsigned nps) {
		init(ps, nps);
	}
	~crt_tree() {
	    if(left)
		delete left;
	    if(right)
		delete right;
	}


	void init(const mpz_class* ps, unsigned nps) {
		n  = nps;
		n2 = n/2;
		assert(n > 0);
		if(n == 1) {
			left  = NULL;
			right = NULL;
			mod   = *ps;
		}
		else {
			left  = new crt_tree(ps, n2);
			right = new crt_tree(ps + n2, n - n2);

			mpz_class g(0);
			mpz_gcdext(g.get_mpz_t(),
					crt_right.get_mpz_t(), crt_left.get_mpz_t(),
					left->mod.get_mpz_t(), right->mod.get_mpz_t());
			assert(g == 1);

			crt_left *= right->mod;
			crt_right *= left->mod;

			mod = left->mod * right->mod;
		}
	}

	inline mpz_class& get_mod(void) { return mod; }
	
	mpz_class do_crt(const mpz_class* cs) {
		if(n == 1)
			return *cs;

		mpz_class val_left, val_right, val;
		
		val_left  = left->do_crt(cs);
		val_right = right->do_crt(cs + n2);

		val  = val_left * crt_left + val_right * crt_right;
		mpz_mod(val.get_mpz_t(), val.get_mpz_t(), mod.get_mpz_t());
		return val;
	}
	
	mpz_class do_crt_oneindex(const mpz_class& c, unsigned i) {
		if(n == 1) {
			assert(i==0);
			return c;
		}

		mpz_class val;
		
		if(i < n2) {
			val = crt_left * left->do_crt_oneindex(c, i);
		}
		else {
			val = crt_right * right->do_crt_oneindex(c, i-n2);
		}

		mpz_mod(val.get_mpz_t(), val.get_mpz_t(), mod.get_mpz_t());
		return val;
	}
};

#endif
