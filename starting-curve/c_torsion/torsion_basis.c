#include <pari/pari.h>
#include <stdlib.h>

//Compile with
//gcc c_torsion/torsion_basis.c c_torsion/sqrt.c c_torsion/field_stuff.c -Wall -Wextra -pedantic -L/usr/local/lib -lpari -o c_torsion/torsion_basis

#define STR_SIZE 1000
#define max_ext 200

GEN frob_eval(GEN, GEN);
GEN frobn_eval(GEN, long, GEN);
GEN get_i(GEN);
GEN embed(GEN, GEN, GEN);
GEN better_sqrt(GEN, GEN);

GEN p;
GEN Ts[max_ext];
GEN frob_mat[max_ext];

FILE *fptr;
GEN sqrtm1;

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

GEN mul_p(GEN P, GEN T) {
  pari_sp av = avma;
  GEN Q = cgetg(3,t_VEC);
  gel(Q, 1) = frobn_eval(gel(P, 1), 2, frob_mat[degree(T)]);
  gel(Q, 2) = FpX_neg(frobn_eval(gel(P, 2), 2, frob_mat[degree(T)]), p);
  return gerepileupto(av, Q);
}

GEN mul_with_pol(GEN P, GEN pol, GEN T) {
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
    frobs[i] = mul_p(frobs[i-1], T);
  }
  GEN Q = ellinf();
  for (long i = 0; i < k+1; i++)
  {
    Q = FpXQE_add(Q, FpXQE_mul(frobs[i], gel(pol, i+2), pol_1(0), T, p), pol_1(0), T, p);
  }

  return gerepileupto(av, Q);
}

GEN mul_using_frob(GEN P, GEN D, GEN T) {
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
  
  Q = mul_with_pol(P, pol, T);
  return gerepilecopy(av, FpXQE_mul(Q, rem, pol_1(0), T, p));
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

int has_order(GEN P, GEN N, GEN T) {

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

int is_independent(GEN P, GEN Q, GEN N, GEN T) {
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

GEN endJ(GEN P, GEN T) {
  GEN i = get_i(T);
  if(i == NULL) return ellinf();
  GEN Q;
  Q = cgetg(3,t_VEC);
  gel(Q, 1) = FpX_neg(gel(P, 1), p);
  gel(Q, 2) = FpXQ_mul(gel(P, 2), i, T, p);
  return Q;
}

GEN endI(GEN P, GEN T) {
  GEN Q = cgetg(3,t_VEC);
  gel(Q, 1) = frob_eval(gel(P, 1), frob_mat[degree(T)]);
  gel(Q, 2) = frob_eval(gel(P, 2), frob_mat[degree(T)]);
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

GEN torsion_point(GEN D, GEN T) {
  pari_timer t;
  while(1) {
    timer_delay(&t);
    GEN point = random_point(T);
    GEN P = mul_using_frob(point, D, T);
    if(has_order(P, D, T)) return P;
  }
}

int torsion_basis(GEN N, GEN T, GEN *P, GEN *Q) {
  *P = torsion_point(N, T);
  //printf("First point found\n");
  //pari_timer t;
  //timer_delay(&t);
  *Q = endJ(*P, T);
  //timer_printf(&t, "for second point from endj");

  int a = has_order(*Q, N, T);
  //timer_printf(&t, "for has_order");

  int b = is_independent(*P, *Q, N, T);
  //timer_printf(&t, "for is_independent");

  if(a && b) return 1;

  *Q = endI(*P, T);
  //printf("Second point from endi\n");
  if(has_order(*Q, N, T) && is_independent(*P, *Q, N, T)) {
    //timer_printf(&t, "for Second point from endi");
    return 1;
  } 
  //timer_printf(&t, "for Second point from endi failed");

  while(1) {
    *Q = torsion_point(N, T);
    //printf("Second point randomly\n");

    if(is_independent(*P, *Q, N, T)) return 1;
  }
}

int main(int argc, char *argv[]) {
    long i = atoi(argv[2]);
    pari_init(i*1000000000, 2);
    pari_timer tim;
    timer_start(&tim);
    char i_filename[STR_SIZE];

    snprintf(i_filename, sizeof(i_filename), "c_torsion/data%s.txt", argv[1]);
    //printf("%s\n", i_filename);
    fptr = fopen(i_filename, "r");
    p = gp_read_stream(fptr);

    sqrtm1 = FpXQ_sqrt(FpX_neg(pol_1(0), p), mkpoln(3, gen_1, gen_0, gen_1), p);
    
    pari_timer t;
    GEN P, Q;
    if(argc >= 5) {
      timer_delay(&t);
      long k = atoi(argv[3]) * 2;
      GEN D = gp_read_str(argv[4]);
      if(argc >= 6) Ts[k] = gp_read_str(argv[5]);
      else Ts[k] = init_Fq(p, k, 0);
      if (k == 2) Ts[k] = mkpoln(3, gen_1, gen_0, gen_1);
      frob_mat[k] = FpX_matFrobenius(Ts[k], p);
      torsion_basis(D, Ts[k], &P, &Q);
      GEN i = get_i(Ts[k]);

      print_pol(Ts[k]);
      print_pol(i);

      print_pol(gel(P, 1));
      print_pol(gel(P, 2));
      print_pol(gel(Q, 1));
      print_pol(gel(Q, 2));
      return 0;
    }


    GEN ntemp = gp_read_stream(fptr); long n = itos(ntemp); cgiv(ntemp);
    long extdeg[n];
    GEN tors[n];

    for (long i = 0; i < n; i++)
    {
      extdeg[i] = itos(gp_read_stream(fptr))*2;
      tors[i] = gp_read_stream(fptr);
    }
    for (long i = 0; i < n; i++)
    {
      timer_delay(&t);
      pari_printf("Finding %Ps torsion in extension %d\n", tors[i], extdeg[i]);
      //pari_printf("Defining polynomial: %Ps\n", Ts[extdeg[i]]);
      Ts[extdeg[i]] = init_Fq(p, extdeg[i], 0);
      frob_mat[extdeg[i]] = FpX_matFrobenius(Ts[extdeg[i]], p);
      torsion_basis(tors[i], Ts[extdeg[i]], &P, &Q);
      timer_printf(&t, "Torsion basis in extension with degree %d", extdeg[i]);
      printf("\n");
    }
    timer_printf(&tim, "overall");
    //Closing stuff
    fclose(fptr);
    pari_close();
}