#include <pari/pari.h>

extern GEN p, sqrtm1;

GEN frob_eval(GEN x, GEN frobm) {
  return FpM_FpC_mul_FpX(frobm, RgX_to_RgC(x, lg(frobm)-1), p, 0);
}

GEN frobn_eval(GEN x, long n, GEN frobm) {
    return frob_eval(x, FpM_powu(frobm, n, p));
}

GEN get_i(GEN T) {
  GEN xff = ffgen(FpX_to_mod(mkpoln(3, gen_1, gen_0, gen_1), p), 0);
  GEN yff = ffgen(FpX_to_mod(T, p), 0);
  gel(xff, 2) = gcopy(sqrtm1);
  GEN m = ffembed(xff, yff);
  return gel(ffmap(m, xff), 2);
}

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