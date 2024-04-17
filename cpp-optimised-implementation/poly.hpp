#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>

#include "fp2.hpp"

typedef fp2_elem *fp2X_elem; //array of fp2_elems

// Currently not working....
// Also, NTL resultant computation is hella slow??
// Good enough for our sizes probably
void Poly_mult_naive(fp2X_elem h, fp2X_elem const &f, int f_size, fp2X_elem const &g, int g_size) {
    for (int i = 0 ; i < f_size + g_size - 1 ; i++) {
        fp2_elem coeff = Fp2_zero();
        for (int j = 0 ; j < std::min(i, g_size) ; j++) {
            if (i-j < f_size) {
                coeff = Fp2_add(coeff, Fp2_mul(f[j], g[i-j]));
            }
        }
        h[i] = coeff;
    }
}

void Poly_add(fp2X_elem h, fp2X_elem const &f, int f_size, fp2X_elem const &g, int g_size) {
    if (f_size < g_size) {
        return Poly_add(h, g, g_size, f, f_size);
    }

    for (size_t i = 0 ; i < g_size ; i++) {
        h[i] = Fp2_add(f[i], g[i]);
    }
    for (size_t i = g_size; i < f_size ; i++) {
        h[i] = f[i];
    }
}

void Poly_neg(fp2X_elem h, fp2X_elem const &f, int f_size) {
    for (int i = 0 ; i < f_size ; i++) {
        h[i] = Fp2_negative(f[i]);
    }
}

void Poly_sub(fp2X_elem h, fp2X_elem const &f, int f_size, fp2X_elem const &g, int g_size) {
    fp2_elem g_min[g_size];
    Poly_neg(g_min, g, g_size);
    Poly_add(h, f, f_size, g_min, g_size);
}

void Poly_scale(fp2X_elem h, int f_size, fp2_elem const &a) {
    for (int i = 0 ; i < f_size ; i++) {
        h[i] = Fp2_mul(a, h[i]);
    }
}

void Poly_shift(fp2X_elem h, fp2X_elem const &f, int f_size, int shift) {
    for (size_t i = 0 ; i < shift ; i++) {
        h[i] = Fp2_zero();
    }
    for (int i = shift ; i < f_size + shift ; i++) {
        h[i] = f[i];
    }
}

bool Poly_equal(fp2X_elem const &f, int f_size, fp2X_elem const &g, int g_size) {
    if (f_size != g_size) {
        return false;
    }
    for (size_t i = 0 ; i < f_size ; i++) {
        if (!(Fp2_equal(f[i], g[i])))
            return false;
    }
    return true;
}

// KARATSUBA MULTIPLICATION //
void Poly_mult(fp2X_elem h, fp2X_elem const &f, int n, fp2X_elem const &g, int m) {
    if (n < m) {
        return Poly_mult(h, g, m, f, n);
    }

    if (m == 1) {
        for (int i = 0 ; i < n ; i++) {
            h[i] = Fp2_mul(f[i], g[0]);
        }
        return;
    }

    if (n == 2) {
        fp2_elem t0, t1, t2;
        t0 = Fp2_mul(f[0], g[0]);
        h[0] = t0;

        t1 = Fp2_add(f[0], f[1]);
        t2 = Fp2_add(g[0], g[1]);
        h[2] = Fp2_mul(f[1], g[1]);
        t1 = Fp2_mul(t1, t2);
        t1 = Fp2_sub(t1, t0);
        h[1] = Fp2_sub(t1, h[2]);
        return;
    }
  
	// Cases for f quadratic and ...
	if ((n == 3)) {
        // ... g linear
		if (m == 2) {
		    fp2_elem t0, t1, t2, t3;
			t0 = Fp2_mul(f[0], g[0]);
			t2 = Fp2_mul(f[1], g[1]);
			t3 = Fp2_add(f[0], f[1]);
			t1 = Fp2_add(g[0], g[1]);
			t1 = Fp2_mul(t1, t3);
			t1 = Fp2_sub(t1, t2);
			t1 = Fp2_sub(t1, t0);
			t3 = Fp2_mul(f[2], g[1]);
            h[0] = t0;
            h[1] = t1;
            h[2] = Fp2_add(Fp2_mul(f[2], g[0]), t2);
            h[3] = t3;
			return;
		}
        // ... g quadratic
        if (m == 3) {
		    fp2_elem t0, t1, t2, fg_high[3];
			
            Poly_mult(fg_high, &(f[1]), n-1, &(f[1]), m-1);

			t0 = Fp2_add(f[0], f[1]);
			t1 = Fp2_add(g[0], g[1]);
			t1 = Fp2_mul(t0, t1);
			t0 = Fp2_add(f[0], f[2]);
			t2 = Fp2_add(g[0], g[2]);
			t2 = Fp2_mul(t0, t2);
			
            h[0] = Fp2_mul(f[0], g[0]);

			t1 = Fp2_sub(t1, h[0]);
            h[1] = Fp2_sub(t1, fg_high[0]);
            t2 = Fp2_sub(t2, h[0]);
            t2 = Fp2_sub(t2, fg_high[2]);
            h[2] = Fp2_add(t2, fg_high[0]);
			h[3] = fg_high[1];
            h[4] = fg_high[2];
      
			return;
        }
    }

    int nf = n >> 1;
    int mf = n - nf;
    int mg = m - nf;
    int lenh = n + m - 1;
    // f = f0 + x^nf * f1
    if (m <= nf) {
        fp2_elem fg_low[nf + m - 1], fg_high[mf + m - 1];
        //fp2X_elem f0 = Poly_shorten(f, nf);
        Poly_mult(fg_low, f, nf, g, m);
		Poly_mult(fg_high, &(f[nf]), mf, g, m);
		for(int i = 0; i < nf; i++)
			h[i] = fg_low[i];
    
		for(int i = nf; i < nf + m - 1; i++)
			h[i] = Fp2_add(fg_high[i-nf], fg_low[i]);

		for(int i = nf + m - 1; i < lenh; i++)
			h[i] = fg_high[i - nf];

		return;
    }

    // f = f0 + x^nf * f1, g = g0 + x^nf g1
    fp2_elem f_mid[mf], g_mid[mf], fg_low[2*nf-1], fg_mid[2*mf-1], fg_high[mf+mg-1];
  
	for(int i = 0; i < nf; i++)
		f_mid[i] = Fp2_add(f[i], f[i+nf]);

	if (n & 1)
		f_mid[nf] = f[n-1];
  
	int i = 0;
	while(i < nf && i < mg)	{
		g_mid[i] = Fp2_add(g[i], g[nf+i]);
		i++;
	}
	while(i < nf) {
		g_mid[i] = g[i];
		i++;
	}
	
	Poly_mult(fg_low, f, nf, g, nf);
	Poly_mult(fg_high, &(f[nf]), mf, &(g[nf]), mg);
	
	if ((n & 1) && (mg == mf))	{
		g_mid[nf] = g[m - 1];
		Poly_mult(fg_mid, f_mid, mf, g_mid, mf);
	} else {
	    Poly_mult(fg_mid, f_mid, mf, g_mid, nf);
	}

	for(int i = 0; i < mf + mg - 1; i++)
		fg_mid[i] = Fp2_sub(fg_mid[i], fg_high[i]);

	for(int i = 0; i < 2*nf-1; i++)
		fg_mid[i] = Fp2_sub(fg_mid[i], fg_low[i]);

	for(int i = 0; i < nf; i++)
		h[i] = fg_low[i];

	for(int i = nf; i < 2*nf-1; i++)
		h[i] = Fp2_add(fg_low[i], fg_mid[i-nf]);

	h[2*nf-1] = fg_mid[nf-1];
  
	for(int i = 2*nf; i < 2*nf+mf-1; i++)
		h[i] = Fp2_add(fg_mid[i-nf], fg_high[i-2*nf]);

	for(int i = 2*nf+mf-1; i < lenh; i++)
		h[i] = fg_high[i-2*nf];
	return;
}

/*
std::vector<fp2_elem> Poly_mult_low(size_t deg, std::vector<fp2_elem> const &f, std::vector<fp2_elem> const &g) {
    // Computes a*b mod X^deg
    int n = f.size();
    int m = g.size();

    if (n < m) {
        return Poly_mult_low(deg, g, f);
    }

    if (deg == 1) {
        return {Fp2_mul(f[0], g[0])};
    }

    if (deg >= n + m - 1) {
        return Poly_mult(f, g);
    }

    if (n > deg) {
        return Poly_mult_low(deg, Poly_shorten(f, deg), g);
    }

    if (m > deg) {
        return Poly_mult_low(deg, f, Poly_shorten(g, deg));
    }

    if (m == 1) {
        return Poly_scale(g[0], f);
    }

    if (deg == 2) {
        fp2X_elem h;
        h.push_back(Fp2_mul(f[0], g[0]));
        h.push_back(Fp2_add(Fp2_mul(f[0], g[1]), Fp2_mul(f[1], g[0])));
        return h;
    }

    int lenf1 = n >> 1;
    int leng1 = m >> 1;
    // f(x) = f0(x^2) + x * f1(x^2)
    // g(x) = g0(x^2) + x * g1(x^2)
	int lenf0 = n - lenf1;
    int leng0 = m - leng1;

    // f = f0 + x^nf * f1, g = g0 + x^nf g1
    fp2X_elem f0, f1;
    for (size_t i = 0 ; i < lenf1 ; i ++) {
        f0.push_back(f[2*i]);
        f1.push_back(f[2*i + 1]);
    }
    if (lenf0 > lenf1) {
        f0.push_back(f[n-1]);
    }

    fp2X_elem g0, g1;
    for (size_t i = 0 ; i < leng1 ; i ++) {
        g0.push_back(g[2*i]);
        g1.push_back(g[2*i + 1]);
    }
    if (leng0 > leng1) {
        g0.push_back(g[m-1]);
    }

    assert (g0.size() == leng0);
    assert (f0.size() == lenf0);
    assert (g1.size() == leng1);
    assert (f1.size() == lenf1);

    int l1 = deg >> 1;
    int l0 = deg - l1;

    int n0 = lenf0 + leng0 - 1;
    if (n0 > l0) 
        n0 = l0;

    int n1 = lenf1 + leng1 - 1;
    if (n1 > l1) 
        n1 = l1;

    int n01 = lenf0 + leng0 -1;
	if(n01 > l1)
		n01 = l1;
    
    // Recursive calls
    std::vector<fp2_elem> fg_low = Poly_mult_low(n0, f0, g0);
    std::vector<fp2_elem> fg_high = Poly_mult_low(n1, f1, g1);
    std::vector<fp2_elem> fg_mid = Poly_mult_low(n01, Poly_add(f0, f1), Poly_add(g0, g1));
    //fg_mid = Poly_sub(fg_mid, Poly_add(fg_low, fg_high));

    fp2X_elem h;

    // Combine results
    h.push_back(fg_low[0]);

    for (size_t i = 1 ; i < deg ; i ++) {
        size_t i2 = i/2;
        if (i % 2 == 0) {
            if (i2 < n1) {
                h.push_back(Fp2_add(fg_low[i2], fg_high[i2-1]));
            } else if (2*n1 < n) {
                assert (i2 == n1);
                h.push_back(Fp2_add(fg_low[i2], fg_high[i2-1]));
            }
        } else {
            if (i2 < n1) {
                fp2_elem t = Fp2_sub(fg_mid[i2], fg_low[i2]);
                h.push_back(Fp2_sub(t, fg_high[i2]));
            } else if (i2 < n01) {
                h.push_back(Fp2_sub(fg_mid[i2], fg_low[i2]));
            } else {
                h.push_back(Fp2_zero());
            }
        }
    }

    return h;
}

std::vector<fp2_elem> Xmin(fp2_elem const &a) {
    return {Fp2_negative(a), Fp2_one()};
}

// Product tree (from SCALLOP code)
std::vector<std::vector<fp2X_elem>> product_tree(std::vector<fp2X_elem> const &leaves) {
    if (leaves.empty())
        throw std::logic_error("no leaves");
    std::vector<std::vector<fp2X_elem>> tree{leaves};
    auto prev = leaves; //<- copies leaves to not modify the original list...
    while (prev.size() > 1) {
        std::vector<fp2X_elem> next;
        {
            for (size_t i = 0; i < prev.size()-1; i += 2) {
                next.push_back(Poly_mult(prev[i], prev[i+1]));
            }
            if (prev.size() % 2)
                next.push_back(prev.back());
        }
        tree.insert(tree.begin(), next);
        prev = next;
    }
    return tree;
};


fp2X_elem Reciprocal(fp2X_elem const &f, int n) {
    // computes h, such that h*f = 1 mod x^n
    assert (!(n == 0));
    fp2X_elem h;
    if (n == 1) {
        h.push_back(Fp2_inv(f[0]));
        return h;
    }
    // More easy cases later...
    int m = n - (n>>1);
    fp2X_elem g = Reciprocal(f, m);
    //fp2X_elem t = Poly_mult_middle(m, n, g, f);
    //t = Poly_mult_middle(m, n, g, t);
    fp2X_elem t = Poly_mult_low(n, g, f);
    t = Poly_mult_low(n, g, t);

    for(size_t i = 0; i < m; i++)
        h.push_back(g[i]);
    for(size_t i = m; i < n; i++)
        h.push_back(Fp2_negative(t[i]));
    return h;
}

fp2X_elem Reverse(fp2X_elem const &f) {
    size_t n = f.size();
    fp2X_elem f_rev;
    for (size_t i = 0 ; i < n ; i++) {
        f_rev.push_back(f[n-i-1]);
    }
    return f_rev;
}

fp2X_elem Poly_mod(fp2X_elem const &g, fp2X_elem const &f, fp2X_elem const &f_rev_inv) {
    if (g.size() < f.size())
        return g;

    //fp2X_elem f_rev = Reverse(f);
    //fp2X_elem f_rev_inv = Reciprocal(f_rev, g.size() - f.size() + 1);

    fp2X_elem g_rev = Reverse(g);
    fp2X_elem Q = Poly_mult_low(g.size() - f.size() + 1, f_rev_inv, g_rev);
    fp2X_elem Q_rev = Reverse(Q);
    fp2X_elem t = Poly_mult_low(f.size() - 1, Q_rev, f);

    fp2X_elem h;
    for (size_t i = 0 ; i < f.size() - 1 ; i++) {
        h.push_back(Fp2_sub(g[i], t[i]));
    }

    return h;
}

std::vector<std::vector<fp2X_elem>> rev_inv_tree(std::vector<std::vector<fp2X_elem>> const &tree, fp2X_elem const &g) {
    std::vector<std::vector<fp2X_elem>> rev_inv_tree;
    std::vector<fp2X_elem> prev = {g};
    bool first = true;
    for (auto const &nodes : tree) {
        std::vector<fp2X_elem> next;
        for (size_t i = 0 ; i < nodes.size(); i++) {
            int gsize = prev[i/2].size()-1;
            if (first) {
                first = false;
                gsize += 1;
            }
            fp2X_elem f_rev = Reverse(nodes[i]);
            if (gsize < f_rev.size()) {
                next.push_back(f_rev); // just add this, won't get used anyway

            } else {
                fp2X_elem f_rev_inv = Reciprocal(f_rev, gsize - f_rev.size() + 1);
                next.push_back(f_rev_inv);
            }
        }
        rev_inv_tree.push_back(next);
        prev = nodes;
    }
    return rev_inv_tree;
};

fp2_elem Fast_Resultant(fp2X_elem const &g, std::vector<std::vector<fp2X_elem>> const &f_tree, std::vector<std::vector<fp2X_elem>> const &f_rev_inv_tree) {
    //Takes in g and a product tree of f
    //computes the product of g mod f_i
    std::vector<fp2X_elem> prev{g};

    for (size_t j = 0 ; j < f_tree.size() ; j++) {
        auto fi_list = f_tree[j];
        auto fi_rev_inv_list = f_rev_inv_tree[j];
        std::vector<fp2X_elem> next;
        for (size_t i = 0 ; i < fi_list.size() ; i++) {
            next.push_back(Poly_mod(prev[i/2], fi_list[i], fi_rev_inv_list[i]));
        }
        prev = next;
    }
    // multiply the leaves
    fp2_elem res = Fp2_one();
    for (auto const& gi : prev) {
        assert (gi.size() == 1);
        res = Fp2_mul(res, gi[0]);
    }

    return res;
}
*/