#include <pari/pari.h>

extern GEN sqrtm1, p;
GEN frob_eval(GEN, GEN);
GEN frobn_eval(GEN, long, GEN);
GEN get_i(GEN);
GEN embed(GEN, GEN, GEN);

//x in F_q^2, where q = p^n, q = 1 mod 4
//return x^((q-1)/4)
//we write (q-1)/4 = ((p-1)/2)*((p+1)/2)*sum_{i=0}^{(n-1}p^i
//and for the p-th power use Frobenius matrix
GEN pow_qm1div4(GEN x, GEN frobm, GEN T){
  long k = degree(T)/2;
  //pari_printf("In pow_qm1div4, q = p^%ld\n", k);
  GEN res = FpXQ_pow(x, diviuexact(addiu(p, 1), 2), T, p);
  res = FpXQ_pow(res, diviuexact(subiu(p, 1), 2), T, p);
  GEN frobs[k/2];
  frobs[0] = gcopy(res);
  GEN frob_mat2 = FpM_mul(frobm,frobm, p);
  for (long i = 1; i < k/2; i++)
  {
    frobs[i] = FpM_FpC_mul_FpX(frob_mat2, RgX_to_RgC(frobs[i-1], degree(T)), p, 0);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  return res;
}
//x in F_q^2, where q = p^n, q = 3 mod 4
//return x^((q-3)/4)
//we write (q-3)/4 = a + p[pa + 3a + 2]*sum_{i=0}^{(n-3)/2}p^{2i}
//where a = (p-3)/4
//and for the p-th power use Frobenius matrix
GEN pow_qm3div4(GEN x, GEN frobm, GEN T){
  long k = degree(T)/2;
  //y = x^a
  GEN y = FpXQ_pow(x, diviuexact(subiu(p, 3), 4), T, p);
  //pari_printf("In pow_qm3div4, q = p^%ld\n", k);
  //res = x^(3*a + 2)
  GEN res = FpXQ_mul(FpXQ_sqr(x, T, p), FpXQ_powu(y, 3, T, p), T, p);
  //res = x^(p*(p*a + 3*a + 2))
  res = frob_eval(FpXQ_mul(res, frob_eval(y, frobm), T, p), frobm);
  GEN frobs[(k-1)/2];
  frobs[0] = gcopy(res);
  GEN frob_mat2 = FpM_mul(frobm,frobm, p);
  for (long i = 1; i < (k-1)/2; i++)
  {
    frobs[i] = FpM_FpC_mul_FpX(frob_mat2, RgX_to_RgC(frobs[i-1], degree(T)), p, 0);
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
GEN pow_qm3div4_odd(GEN x, GEN frobm, GEN T){
  long k = degree(T);
  //y = x^a
  GEN y = FpXQ_pow(x, diviuexact(subiu(p, 3), 4), T, p);
  //pari_printf("In pow_qm3div4, q = p^%ld\n", k);
  //res = x^(3*a + 2)
  GEN res = FpXQ_mul(FpXQ_sqr(x, T, p), FpXQ_powu(y, 3, T, p), T, p);
  //res = x^(p*(p*a + 3*a + 2))
  res = frob_eval(FpXQ_mul(res, frob_eval(y, frobm), T, p), frobm);
  GEN frobs[(k-1)/2];
  frobs[0] = gcopy(res);
  GEN frob_mat2 = FpM_mul(frobm,frobm, p);
  for (long i = 1; i < (k-1)/2; i++)
  {
    frobs[i] = FpM_FpC_mul_FpX(frob_mat2, RgX_to_RgC(frobs[i-1], degree(T)), p, 0);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  res = FpXQ_mul(res, y, T, p);
  return res;
}
//x in F_q^2, where q = p^n
//return x^((q-1)/2)
GEN pow_qm1div2(GEN x, GEN frobm, GEN T){
  long k = degree(T)/2;
  //pari_printf("In pow_qm1div2, q = p^%ld\n", k);
  GEN res = FpXQ_pow(x, diviuexact(subiu(p, 1), 2), T, p);
  GEN frobs[k];
  frobs[0] = gcopy(res);
  for (long i = 1; i < k; i++)
  {
    frobs[i] = frob_eval(frobs[i-1], frobm);
    res = FpXQ_mul(res, frobs[i], T, p);
  }
  return res;
}


GEN sqrt_even(GEN, GEN);
/*Algorithm 9 from https://eprint.iacr.org/2012/685.pdf
a is in F_q^2 with q = 3 mod 4*/
GEN sqrt_even2(GEN a, GEN T) {
  GEN i, a1, alpha, x0, b;
  i = get_i(T);
  GEN frobm = FpX_matFrobenius(T, p);
  //long n = degree(T)/2;
  a1 = pow_qm3div4(a, frobm, T);
  x0 = FpXQ_mul(a, a1, T, p);
  alpha = FpXQ_mul(a1, x0, T, p);

  if(gequal1(FpX_neg(alpha, p))) {
    return FpXQ_mul(i, x0, T, p);
  }
  b = pow_qm1div2(FpX_add(pol_1(0), alpha, p), frobm, T);
  return FpXQ_mul(b, x0, T, p);
}
//a is in F_q, where q = p^n, q = 3 mod 4
GEN sqrt_odd(GEN a, GEN T) {
  //printf("Odd sqrt\n");
  GEN frobm = FpX_matFrobenius(T, p);
  GEN a1 = pow_qm3div4_odd(a, frobm, T);
  return FpXQ_mul(a1, a, T, p);
}

//a has to be a quadratic residue
GEN better_sqrt(GEN a, GEN T) {
  pari_sp av = avma;
  //printf("In better_sqrt with extension: %ld\n", degree(T));
  if(degree(T) % 4 == 0) return gerepileupto(av, sqrt_even(a, T));
  if(degree(T) % 2 == 0) return sqrt_even2(a, T);
  //return FpXQ_sqrt(a, T, p);
  return gerepileupto(av, sqrt_odd(a, T));
}
/*Algorithm 10 from https://eprint.iacr.org/2012/685.pdf
a is in F_q^2 with q = 1 mod 4*/
GEN sqrt_even(GEN a, GEN T) {
  long n = degree(T)/2;
  //GEN q = powiu(p, n);
  GEN frobm = FpX_matFrobenius(T, p);

  //printf("In pari sqrt_even, with q = p^%ld\n", n);
  GEN b = pow_qm1div4(a, frobm, T);
  //GEN b = FpXQ_pow(a, divis(subiu(q, 1), 4), T, p);
  //if(gequal(b,b_other)) printf("bbbbb\n");
  GEN x;
  GEN bq = frobn_eval(b, n, frobm);
  GEN b2 = FpXQ_sqr(b, T, p);

  if(gequal1(FpXQ_mul(bq, b, T, p))) {
    //printf("%ld, d1\n", n);
    GEN new_a = FpXQ_mul(b2, a, T, p);
    GEN minpoly = FpXQ_minpoly(new_a, T, p);

    GEN sqrt_in_smaller = better_sqrt(pol_x(0), minpoly);
    //pari_printf("%Ps\n\n%Ps\n\n", embed(pol_x(0), new_a, T), new_a);
    GEN x0 = embed(sqrt_in_smaller, new_a, T);
    //if(gequal(FpXQ_sqr(x0, T, p), new_a)) printf("Recursion works\n");
    //else printf("Recursion doesn't work, %ld\n", n);

    x = FpXQ_mul(x0, bq, T, p);
  }
  else {
    //printf("%ld, d2\n", n);
    //"Precomputation", which is only used here
    GEN c;
    do {
      c = random_FpX(n*2, 0, p);
    } while (FpXQ_issquare(c, T, p));

    GEN d = pow_qm1div2(c, frobm, T);
    GEN dc = FpXQ_mul(c, d, T, p);
    GEN e = FpXQ_inv(dc, T, p);
    GEN f = FpXQ_sqr(dc, T, p);
    
    GEN new_a = FpXQ_mul(f, FpXQ_mul(b2, a, T, p), T, p);
    GEN minpoly = FpXQ_minpoly(new_a, T, p);

    //pari_printf("a_in_smaller issquare: %d\n", FpXQ_issquare(a_in_smaller, Ts[n], p));
    GEN sqrt_in_smaller = better_sqrt(pol_x(0), minpoly);
    GEN x0 = embed(sqrt_in_smaller, new_a, T);
    //if(gequal(FpXQ_sqr(x0, T, p), new_a)) printf("Recursion works, %ld\n", n);
    //else printf("Recursion doesn't work, %ld\n", n);
    
    x = FpXQ_mul(e, FpXQ_mul(x0, bq, T, p), T, p);
  }
  return x;
}