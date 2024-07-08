#include <pari/pari.h>

extern GEN sqrtm1, p, frob_pol, T;

extern GEN frob_pol_pow[];

GEN frob_i(long);

#define frob_eval(x, i) (FpX_FpXQ_eval((x), frob_i((i)), T, p))

GEN get_i(GEN);
//Let m be the minimal polynomial of a in F_p[x]/(T) | F_p
//y is in F_p[x]/(m)
//return y embedded in F_p[x]/(T)
GEN embed(GEN y, GEN a, GEN T) {
  pari_sp av = avma;
  GEN z = gcopy(pol_1(0)), res = pol_0(0);
  //pari_printf("y: %Ps\n\n", y);
  //pari_printf("a: %Ps\n\n", a);

  for (long i = 2; i < lg(y); i++)
  {
    //pari_printf("gel(y, %ld): %Ps\n\n", i, gel(y, i));
    //pari_printf("z: %Ps\n\n", z);
    //pari_printf("Res before: %Ps\n\n", res);
    res = FpX_add(res, FpX_Fp_mul(z, gel(y, i), p), p);
    //pari_printf("Res after: %Ps\n\n", res);
    z = FpXQ_mul(z, a, T, p);
  }
  return gerepileupto(av, res);
}

//x in F_q^2, where q = p^n, q = 1 mod 4
//return x^((q-1)/4)
//we write (q-1)/4 = ((p-1)/2)*((p+1)/2)*sum_{i=0}^{(n-1}p^i
//and for the p-th power use Frobenius matrix
GEN pow_qm1div4(GEN x, GEN T){
  long k = degree(T)/2;
  GEN res = FpXQ_pow(x, diviuexact(addiu(p, 1), 2), T, p);
  res = FpXQ_pow(res, diviuexact(subiu(p, 1), 2), T, p);
  GEN frobs[k/2];
  frobs[0] = gcopy(res);
  for (long i = 1; i < k/2; i++)
  {
    frobs[i] = frob_eval(frobs[i-1], 2);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  return res;
}
//x in F_q^2, where q = p^n, q = 3 mod 4
//return x^((q-3)/4)
//we write (q-3)/4 = a + p[pa + 3a + 2]*sum_{i=0}^{(n-3)/2}p^{2i}
//where a = (p-3)/4
//and for the p-th power use Frobenius matrix
GEN pow_qm3div4(GEN x, GEN T){
  long k = degree(T)/2;
  //y = x^a
  GEN y = FpXQ_pow(x, diviuexact(subiu(p, 3), 4), T, p);
  //res = x^(3*a + 2)
  GEN res = FpXQ_mul(FpXQ_sqr(x, T, p), FpXQ_powu(y, 3, T, p), T, p);
  //res = x^(p*(p*a + 3*a + 2))
  res = frob_eval(FpXQ_mul(res, frob_eval(y, 1), T, p), 1);
  GEN frobs[(k-1)/2];
  frobs[0] = gcopy(res);
  for (long i = 1; i < (k-1)/2; i++)
  {
    frobs[i] = frob_eval(frobs[i-1], 2);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  res = FpXQ_mul(res, y, T, p);
  return res;
}
//x in F_q, where q = p^n, q = 3 mod 4 (that is, n is odd)
//return x^((q-3)/4)
//we write (q-3)/4 = a + p[pa + 3a + 2]*sum_{i=0}^{(n-3)/2}p^{2i}
//where a = (p-3)/4
//and for the p-th power use Frobenius matrix
GEN pow_qm3div4_odd(GEN x, GEN T){
  long k = degree(T);
  //y = x^a
  GEN y = FpXQ_pow(x, diviuexact(subiu(p, 3), 4), T, p);
  GEN res = FpXQ_mul(FpXQ_sqr(x, T, p), FpXQ_powu(y, 3, T, p), T, p);
  //res = x^(p*(p*a + 3*a + 2))
  res = frob_eval(FpXQ_mul(res, frob_eval(y, 1), T, p), 1);
  GEN frobs[(k-1)/2];
  frobs[0] = gcopy(res);
  for (long i = 1; i < (k-1)/2; i++)
  {
    frobs[i] = frob_eval(frobs[i-1], 2);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  res = FpXQ_mul(res, y, T, p);
  return res;
}
//x in F_q^2, where q = p^n
//return x^((q-1)/2)
GEN pow_qm1div2(GEN x, GEN T){
  long k = degree(T)/2;
  GEN res = FpXQ_pow(x, diviuexact(subiu(p, 1), 2), T, p);
  GEN frobs[k];
  frobs[0] = gcopy(res);
  for (long i = 1; i < k; i++)
  {
    frobs[i] = frob_eval(frobs[i-1], 1);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  return res;
}


GEN sqrt_even(GEN, GEN);
/*Algorithm 9 from https://eprint.iacr.org/2012/685.pdf
a is in F_q^2 with q = 3 mod 4*/
GEN sqrt_even2(GEN a, GEN T) {
  GEN a1, alpha, x0, b;
  //long n = degree(T)/2;
  a1 = pow_qm3div4(a, T);
  x0 = FpXQ_mul(a, a1, T, p);
  alpha = FpXQ_mul(a1, x0, T, p);

  if(gequal1(FpX_neg(alpha, p))) {
    return FpXQ_mul(sqrtm1, x0, T, p);
  }
  b = pow_qm1div2(FpX_add(pol_1(0), alpha, p), T);
  return FpXQ_mul(b, x0, T, p);
}
//a is in F_q, where q = p^n, q = 3 mod 4
GEN sqrt_odd(GEN a, GEN T) {
  GEN a1 = pow_qm3div4_odd(a, T);
  return FpXQ_mul(a1, a, T, p);
}

//a has to be a quadratic residue
GEN better_sqrt2(GEN a, GEN T) {
  pari_sp av = avma;
  if(degree(T) % 4 == 0) return gerepileupto(av, sqrt_even(a, T));
  if(degree(T) % 2 == 0) return sqrt_even2(a, T);
  return gerepileupto(av, sqrt_odd(a, T));
}
/*Algorithm 10 from https://eprint.iacr.org/2012/685.pdf
a is in F_q^2 with q = 1 mod 4*/
GEN sqrt_even(GEN a, GEN T) {
  long n = degree(T)/2;

  GEN b = pow_qm1div4(a, T);
  GEN x;
  GEN bq = frob_eval(b, n);
  GEN b2 = FpXQ_sqr(b, T, p);

  if(gequal1(FpXQ_mul(bq, b, T, p))) {
    GEN new_a = FpXQ_mul(b2, a, T, p);
    GEN minpoly = FpXQ_minpoly(new_a, T, p);

    GEN sqrt_in_smaller = better_sqrt2(pol_x(0), minpoly);
    GEN x0 = embed(sqrt_in_smaller, new_a, T);

    x = FpXQ_mul(x0, bq, T, p);
  }
  else {
    //"Precomputation", which is only used here
    GEN c;
    do {
      c = random_FpX(n*2, 0, p);
    } while (FpXQ_issquare(c, T, p));

    GEN d = pow_qm1div2(c, T);
    GEN dc = FpXQ_mul(c, d, T, p);
    GEN e = FpXQ_inv(dc, T, p);
    GEN f = FpXQ_sqr(dc, T, p);
    
    GEN new_a = FpXQ_mul(f, FpXQ_mul(b2, a, T, p), T, p);
    GEN minpoly = FpXQ_minpoly(new_a, T, p);

    GEN sqrt_in_smaller = better_sqrt2(pol_x(0), minpoly);
    GEN x0 = embed(sqrt_in_smaller, new_a, T);
    
    x = FpXQ_mul(e, FpXQ_mul(x0, bq, T, p), T, p);
  }
  return x;
}


GEN alpha(GEN a, long j, GEN T) {
  GEN x, y, z;

  pari_sp av = avma;
  z = gen_0;
  y = gcopy(a);
  x = gcopy(a);
  
  for (long i = 1; i < j+1; i++)
  {
    x = frob_eval(x, 1);
    y = FpXQ_mul(y, x, T, p);
    z = FpX_add(z, y, p);
    if(gc_needed(av, 1)) {
      gerepileall(av, 3, &x, &y, &z);
    }
  }
  return gerepileupto(av, z);
}
//Algorithm 3 from https://community.ams.org/journals/mcom/2014-83-285/S0025-5718-2013-02715-9/S0025-5718-2013-02715-9.pdf
GEN better_sqrt(GEN a, GEN T) {
  GEN c, b, new_a, lam, beta;
  pari_sp av = avma;
  do
  {
    c = random_FpX(degree(T), 0, p);
    new_a = FpXQ_mul(a, FpXQ_sqr(c, T, p), T, p);
    lam = FpXQ_pow(new_a, divis(subis(p, 1), 2), T, p);
    b = FpX_add(FpX_add(pol_1(0), lam, p), alpha(lam, degree(T)-2, T), p);
  } while (gequal0(b));
  GEN x = FpXQ_mul(new_a, FpXQ_sqr(b, T, p), T, p);
  beta = Fp_sqrt(gel(x, 2), p);
  return gerepilecopy(av, FpX_Fp_mul(FpXQ_inv(FpXQ_mul(b, c, T, p), T, p), beta, p));
}