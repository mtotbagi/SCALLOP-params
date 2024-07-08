#include <pari/pari.h>
#include <stdlib.h>
#include <strings.h>
//Compile with
//gcc c_torsion/torsion_basis.c c_torsion/sqrt.c -Wall -Wextra -pedantic -L/usr/local/lib -lpari -o c_torsion/torsion_basis

#define max_frob 1000

GEN better_sqrt(GEN, GEN);

GEN p, T, sqrtm1, frob_pol;
GEN frob_pol_pow[max_frob];

GEN frob_i(long i) {
  if(i == 0) return NULL;
  if(i == 1) return frob_pol;
  if(frob_pol_pow[i] == NULL) {
    frob_pol_pow[i] = FpX_FpXQ_eval(frob_i(i-1), frob_pol, T, p);
  }
  return frob_pol_pow[i];
}

#define frob_eval(x, i) (FpX_FpXQ_eval((x), frob_i((i)), T, p))

void print_pol(GEN f) {
  if(!signe(f)) {
    printf("[0]\n");
    return;
  }

  long len = lg(f);
  printf("[");
  for (long i = 2; i < len; i++)
  {
    pari_printf("%Ps", gel(f,i));
    if(i < len-1) printf(",");
  }
  
  printf("]\n");
  return;
}

GEN dbl(GEN P, GEN a4) {
  pari_sp av = avma;
  GEN x, y, Q, slope;
  if (ell_is_inf(P) || !signe(gel(P,2))) return ellinf();
  x = gel(P,1); y = gel(P,2);
  slope = FpXQ_div(FpX_add(FpX_mulu(FpXQ_sqr(x, T, p), 3, p), a4, p),
                            FpX_mulu(y, 2, p), T, p);
  Q = cgetg(3,t_VEC);
  gel(Q, 1) = FpX_sub(FpXQ_sqr(slope, T, p), FpX_mulu(x, 2, p), p);
  gel(Q, 2) = FpX_sub(FpXQ_mul(slope, FpX_sub(x, gel(Q, 1), p), T, p), y, p);
  return gerepilecopy(av, Q);
}

GEN add(GEN P, GEN Q, GEN a4)
{
  pari_sp av = avma;
  GEN Px, Py, Qx, Qy, R, slope;
  if (ell_is_inf(P)) return Q;
  if (ell_is_inf(Q)) return P;
  Px = gel(P,1); Py = gel(P,2);
  Qx = gel(Q,1); Qy = gel(Q,2);
  if (ZX_equal(Px, Qx))
  {
    if (ZX_equal(Py, Qy))
      return dbl(P, a4);
    else
      return ellinf();
  }
  slope = FpXQ_div(FpX_sub(Py, Qy, p), FpX_sub(Px, Qx, p), T, p);
  R = cgetg(3,t_VEC);
  gel(R, 1) = FpX_sub(FpX_sub(FpXQ_sqr(slope, T, p), Px, p), Qx, p);
  gel(R, 2) = FpX_sub(FpXQ_mul(slope, FpX_sub(Px, gel(R, 1), p), T, p), Py, p);
  return gerepilecopy(av, R);
}

GEN mul_p(GEN P) {

  pari_sp av = avma;
  GEN Q = cgetg(3,t_VEC);
  gel(Q, 1) = frob_eval(gel(P, 1), 2);
  gel(Q, 2) = FpX_neg(frob_eval(gel(P, 2), 2), p);
  return gerepileupto(av, Q);
}

GEN mul_pl(GEN P, long l) {
  pari_sp av = avma;
  GEN Q = cgetg(3,t_VEC);
  gel(Q, 1) = frob_eval(gel(P, 1), 2*l);
  gel(Q, 2) = FpX_Fp_mul(frob_eval(gel(P, 2), 2*l), powiu(gen_m1, l), p);
  return gerepileupto(av, Q);
}

GEN mul_with_pol(GEN P, GEN pol) {
  pari_sp av = avma;
  //pari_printf("in mul with pol\n%Ps\n", pol);
  if(gequal1(pol)) {
    return gcopy(P);
  }
  long k = degree(pol);
  GEN frobs[k+1];
  frobs[0] = gcopy(P);
  for (long i = 1; i < k+1; i++)
  {
    frobs[i] = mul_p(frobs[i-1]);
  }
  GEN Q = ellinf();
  for (long i = 0; i < k+1; i++)
  {
    Q = add(Q, FpXQE_mul(frobs[i], gel(pol, i+2), pol_1(0), T, p), pol_1(0));
  }

  return gerepileupto(av, Q);
}

GEN base_n(GEN x, GEN n) {
  //pari_printf("In base calc, \n%Ps\n%Ps\n", x, n);
  pari_sp av = avma;
  GEN y = gcopy(x);
  long l = 1;
  while(gcmp(y, n) != -1) {
    l++;
    y = divii(y, n);
  }
  y = gerepilecopy(av, x);
  GEN res = cgetg(l+1, t_VEC);
  for (long i = 1; i < l+1; i++)
  {
    gel(res, i) = remii(y, n);
    y = divii(y, n);
  }
  return gerepilecopy(av, res);
}

GEN simul_mul(GEN Ps, GEN cs) {
  //pari_printf("In mul multi\n");
  long k = lg(Ps) - 1;
  long s = itos(powuu(2, k));
  GEN pres[s], A = ellinf();
  pres[0] = ellinf();
  //printf("1 checkpoint\n");
  for (long i = 1; i < s; i++)
  {
    long j = i;
    /*find the largest non-zero bit*/
    long r = 0;
    while(j > 0)
    {
      j = j >> 1;
      r++;
    } 
    //pari_printf("Current point index in base 2: %Ps\n", base_n(stoi(i), gen_2));
    pres[i] = add(pres[i - (1 << (r-1))], gel(Ps, r), pol_1(0));
    //pari_printf("%Ps\n\n", pres[i]);
  }

  pari_sp av = avma;
  /*Find the longest coefficient*/
  long wordlength = 0;
  for (long i = 1; i < lg(cs); i++)
  {
    if(wordlength < lg(gel(cs, i))) wordlength = lg(gel(cs, i));
  }
  wordlength -= 2;
  //printf("3 checkpoint, wordlength: %ld\n", wordlength);
  
  /*In each iteration find the j-th bit of the i-th word of each coef, add the corresponding point from pres*/
  unsigned long words[k+1];
  for (long i = wordlength-1; i >= 0; i--)
  {
    //printf("i is: %ld\n", i);
    for (long j = 1; j < lg(cs); j++)
    {
      //printf("%ld, %ld, %ld\n", i, j, lg(cs));
      //pari_printf("%ld-th coef is: %Ps, its i-th word: %lu\n", j, gel(cs, j), *int_W(gel(cs, j), i));
      words[j] = lg(gel(cs, j)) < i+3 ? 0 : *int_W(gel(cs, j), i);
      //printf("jth word is: %lu\n", words[j]);
    }
    
    for (long j = 63; j >= 0; j--)
    {
      A = dbl(A, pol_1(0));

      /*Which pres we need*/
      long l = 0;
      for (long t = 1; t < k+1; t++)
      {
        //printf("%ld, %ld, %ld, %ld\n", i, j, words[t], ((words[t] & (1UL << j)) >> j) << (t-1));
        l += ((words[t] & (1UL << j)) >> j) << (t-1);
      }
      //printf("%lu, %lu, l is: %lu\n", i, j, l);
      A = add(A, pres[l], pol_1(0));
    }
    
    gerepileall(av, 1, &A);
  }
  return A;
}

GEN mul_using_frob(GEN P, GEN D) {
  //printf("\nIn mul using frob\n");
  pari_sp av = avma;
  GEN Q;
  long extdeg = degree(T) / 2;
  GEN rem, pol = pol_1(0);
  if(extdeg % 4 == 0){
    rem = diviiexact(polcyclo_eval(extdeg, p), D);
    if(!gequal(polcyclo_eval(extdeg, p), mulii(D, rem)))
    {
      fprintf(stderr, "Bad torsion given, can't use frobenius trick\n");
      exit(1);
    }
    GEN degs = divisors(stoi(extdeg));
    for (long i = 1; i < lg(degs); i++)
    {
      if(itos(gel(degs, i)) == extdeg) continue;
      pol = ZX_mul(pol, polcyclo(itos(gel(degs, i)), 0));
    }
  }
  else if(extdeg % 2 == 1){
    rem = diviiexact(polcyclo_eval(2*extdeg, p), D);
    if(!gequal(polcyclo_eval(extdeg*2, p), mulii(D, rem)))
    {
      fprintf(stderr, "Bad torsion given, can't use frobenius trick\n");
      exit(1);
    }
    GEN degs = divisors(stoi(extdeg));
    for (long i = 1; i < lg(degs); i++)
    {
      if(itos(gel(degs, i)) == extdeg) continue;
      pol = ZX_mul(pol, polcyclo(2*itos(gel(degs, i)), 0));
    }
  }
  else {
    rem = diviiexact(polcyclo_eval(extdeg/2, p), D);
    if(!gequal(polcyclo_eval(extdeg/2, p), mulii(D, rem)))
    {
      fprintf(stderr, "Bad torsion given, can't use frobenius trick\n");
      exit(1);
    }
    GEN degs = divisors(stoi(extdeg));

    for (long i = 1; i < lg(degs); i++)
    {
      if(itos(gel(degs, i)) == extdeg / 2) continue;
      pol = ZX_mul(pol, polcyclo(itos(gel(degs, i)), 0));
    }
  }
  Q = mul_with_pol(P, pol);
  //if(gcmp(rem, p) == -1) return gerepilecopy(av, FpXQE_mul(Q, rem, pol_1(0), T, p));
  long rem_size = extdeg - degree(pol);
  long l = rem_size > 10 ? rem_size / 10 + 1 : 1;
  long k = (rem_size - 1) / l + 1;
  GEN cs = base_n(rem, powiu(p, l));
  /*GEN c = gen_0;
  for (long i = 1; i < lg(cs); i++)
  {
    c = addii(c, mulii(gel(cs, i), powiu(p, (i-1)*l)));
  }
  if(!gequal(c, rem)) fprintf(stderr, "base_n not working\n");*/
  GEN Ps = cgetg(k+1, t_VEC);
  gel(Ps, 1) = gcopy(Q);
  for (long i = 2; i < k+1; i++)
  {
    gel(Ps, i) = mul_pl(gel(Ps, i-1), l);
    //if(!gequal(gel(Ps, i), FpXQE_mul(Q, powiu(p, l*(i-1)), pol_1(0), T, p))) printf("Ps not working\n");
  }
  //pari_printf("%ld, %ld, %ld\n", rem_size, l, k);
  //if(!gequal(simul_mul(Ps, cs), FpXQE_mul(Q, rem, pol_1(0), T, p))) printf("noooooooooooooooooooooooooo\n");
  return gerepilecopy(av, simul_mul(Ps, cs));
  //return gerepilecopy(av, FpXQE_mul(Q, rem, pol_1(0), T, p));
}

GEN radical(GEN N) {
  GEN Nfac = gel(Z_factor(N), 1);
  GEN Nrad = cgetg(3, t_MAT);
  long i, l1 = lg(Nfac);
  gel(Nrad,1) = cgetg(l1, t_COL);
  gel(Nrad,2) = cgetg(l1, t_COL);

  for (i = 1; i < l1; ++i)
  {
    gcoeff(Nrad, i, 1) = gel(Nfac, i);
    gcoeff(Nrad, i, 2) = stoi(1);
  }
  setlg(Nrad[1],i); setlg(Nrad[2],i);
  return factorback(Nrad);
}

int has_order(GEN P, GEN N) {

  GEN Nfac = gel(Z_factor(N), 1);
  long i, l1 = lg(Nfac);
  GEN Nrad = radical(N);
  GEN N1 = diviiexact(N, Nrad);
  GEN Psmall = FpXQE_mul(P, N1, pol_1(0), T, p);
  for (i = 1; i < l1; i++)
  {
    if(gequal(FpXQE_mul(Psmall, diviiexact(Nrad, gel(Nfac, i)), pol_1(0), T, p), ellinf())) return 0;
  }
  return gequal(FpXQE_mul(P, N, pol_1(0), T, p), ellinf());
}

int is_independent(GEN P, GEN Q, GEN N) {
  GEN ell5 = mkvecn(5, gen_0, gen_0, gen_0, gen_1, gen_0);
  GEN ff = cgetg(5, t_FFELT);
  ff = ffgen(FpX_to_mod(T, p), 0);
  GEN E = ellinit(ell5, ff, 0);

  GEN Nfac = gel(Z_factor(N), 1);
  long i, l1 = lg(Nfac);
  GEN Nrad = radical(N);
  GEN N1 = diviiexact(N, Nrad);
  GEN Psmall = FpXQE_mul(P, N1, pol_1(0), T, p);
  GEN Qsmall = FpXQE_mul(Q, N1, pol_1(0), T, p);
  for (i = 1; i < l1; i++)
  {
    GEN Pi = FpXQE_mul(Psmall, diviiexact(Nrad, gel(Nfac, i)), pol_1(0), T, p);
    GEN Qi = FpXQE_mul(Qsmall, diviiexact(Nrad, gel(Nfac, i)), pol_1(0), T, p);
    if (gequal1(ellweilpairing(E, Pi, Qi, gel(Nfac, i)))) {
      return 0;
    }
  }
  return 1;
}

GEN endJ(GEN P) {
  GEN Q;
  Q = cgetg(3,t_VEC);
  gel(Q, 1) = FpX_neg(gel(P, 1), p);
  gel(Q, 2) = FpXQ_mul(gel(P, 2), sqrtm1, T, p);
  return Q;
}

GEN endI(GEN P) {
  GEN Q = cgetg(3,t_VEC);
  gel(Q, 1) = frob_eval(gel(P, 1), 1);
  gel(Q, 2) = frob_eval(gel(P, 2), 1);
  return Q;
}

GEN random_point(GEN T) {
  pari_sp ltop = avma;
  GEN x, x2, y, rhs;
  long v = get_FpX_var(T), d = get_FpX_degree(T);
  do
  {
    set_avma(ltop);
    x   = random_FpX(d,v,p);
    x2  = FpXQ_sqr(x, T, p);
    rhs = FpXQ_mul(x, FpX_add(x2, pol_1(0), p), T, p);
  } while ((!signe(rhs) && !signe(FpX_add(FpX_mulu(x2,3,p), pol_1(0), p)))
          || !FpXQ_issquare(rhs, T, p));
  y = better_sqrt(rhs, T);
  if (!y) pari_err_PRIME("random_FpE", p);
  return gerepilecopy(ltop, mkvec2(x, y));
}

GEN torsion_point(GEN D) {
  pari_timer t;
  while(1) {
    timer_delay(&t);
    GEN point = random_point(T);
    GEN P = mul_using_frob(point, D);
    if(has_order(P, D)) return P;
  }
}

int torsion_basis(GEN N, GEN *P, GEN *Q) {
  *P = torsion_point(N);
  //printf("First point found\n");
  //pari_timer t;
  //timer_delay(&t);
  *Q = endJ(*P);
  //timer_printf(&t, "for second point from endj");

  int a = has_order(*Q, N);
  //timer_printf(&t, "for has_order");

  int b = is_independent(*P, *Q, N);
  //timer_printf(&t, "for is_independent");

  if(a && b) return 1;

  *Q = endI(*P);
  //printf("Second point from endi\n");
  if(has_order(*Q, N) && is_independent(*P, *Q, N)) {
    //timer_printf(&t, "for Second point from endi");
    return 1;
  } 
  //timer_printf(&t, "for Second point from endi failed");

  while(1) {
    *Q = torsion_point(N);
    //printf("Second point randomly\n");

    if(is_independent(*P, *Q, N)) return 1;
  }
}

int main(int argc, char *argv[]) {
    
    pari_init(1000000000, 2);
    pari_timer tim;
    timer_start(&tim);
    
    if(argc < 4) {
      /*fprintf(stderr, "Too few parameters\n");
      return 1;*/
      char p_str[] = "596189944777372638172818300422290385024394468583357632698904620718401201446423451780730293885911387310424468466907078162034493044851014467790119468529446035454449745387614406711198662080037561127345507482289795632805587649490909091103733371042789817522706777642319645749180265150818764030351009256231058220919495563789230174870643587414401676910226850987845426748837131252505799030848442553318641506513527604875273309202232111785774057997494359022678091161441638118008960288159517233147661046896608323640602694819387549541102995842204231072991679669787422227224111599861798724874762075043130922161141224785272355703617441872269175696231180105996744064900599544211922383800484570456691509756280863371782883761347598535069684820049133567";
      long k = 42;
      p = gp_read_str(p_str);
      GEN a = stoi(44);
      Fp_sqrt(a, p);
      /*paristack_setsize(k * gsizebyte(p)*5000, k * gsizebyte(p)*5000 * 2);
      p = gp_read_str(p_str);
      GEN P, Q, D = gp_read_str("43");
      T = init_Fq(p, k, 0);
      frob_pol = FpX_Frobenius(T, p);
      for (long i = 0; i < max_frob; i++)
      {
        frob_i(i);
      }
      sqrtm1 = better_sqrt(mkpoln(1, gen_m1), T);

      torsion_basis(D, &P, &Q);*/
      fprintf(stderr, "Time for torsion basis: %lu\n", timer_get(&tim));
      return 0;
    }
    p = gp_read_str(argv[2]);
    long k = atoi(argv[3]);

    for (long i = 0; i < max_frob; i++)
    {
      frob_pol_pow[i] = NULL;
    }
    

    if(!strcmp(argv[1], "get_field")) {
      paristack_setsize(k * gsizebyte(p)*1000, k * gsizebyte(p)*1000 * 2);
      p = gp_read_str(argv[2]);

      T = init_Fq(p, k, 0);
      frob_pol = FpX_Frobenius(T, p);
      sqrtm1 = better_sqrt(mkpoln(1, gen_m1), T);
      print_pol(T);
      print_pol(sqrtm1);
      print_pol(frob_pol);
      fprintf(stderr, "Time for field generation: %lu\n", timer_get(&tim));
      return 0;
    }

    if(!strcmp(argv[1], "get_torsion_basis")) {
      paristack_setsize(k * gsizebyte(p)*10000, k * gsizebyte(p)*10000 * 2);
      p = gp_read_str(argv[2]);

      GEN P, Q, D = gp_read_str(argv[4]);
      T = gp_read_str(argv[5]);
      sqrtm1 = gp_read_str(argv[6]);
      frob_pol = gp_read_str(argv[7]);
      for (long i = 0; i < max_frob; i++)
      {
        frob_i(i);
      }
      torsion_basis(D, &P, &Q);

      print_pol(gel(P, 1));
      print_pol(gel(P, 2));
      print_pol(gel(Q, 1));
      print_pol(gel(Q, 2));
      fprintf(stderr, "Time for torsion basis: %lu\n", timer_get(&tim));

      return 0;
    }

    if(!strcmp(argv[1], "test")) {
      paristack_setsize(k * gsizebyte(p)*5000, k * gsizebyte(p)*5000 * 2);
      p = gp_read_str(argv[2]);

      GEN P, Q, D = gp_read_str(argv[4]);
      T = init_Fq(p, k, 0);
      frob_pol = FpX_Frobenius(T, p);
      for (long i = 0; i < max_frob; i++)
      {
        frob_i(i);
      }
      sqrtm1 = better_sqrt(mkpoln(1, gen_m1), T);
      fprintf(stderr, "Guess: %lu\n", k * gsizebyte(p));

      torsion_basis(D, &P, &Q);
      fprintf(stderr, "Time for torsion basis: %lu\n", timer_get(&tim));

    }
    pari_close();
}