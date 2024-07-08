#include <pari/pari.h>
#include <stdlib.h>

//Compile with
//gcc c_torsion/testing.c c_torsion/sqrt.c -Wall -Wextra -pedantic -L/usr/local/lib -lpari -o c_torsion/testing

#define STR_SIZE 1000
#define max_ext 200

GEN frob_eval(GEN, GEN);
GEN frobn_eval(GEN, long, GEN);
GEN get_i(GEN);
GEN embed(GEN, GEN, GEN);
GEN better_sqrt(GEN, GEN);
GEN better_sqrt2(GEN, GEN);

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
    
    for (long i = 2; i < 200; i++)
    {
        GEN a, x, y, T;
        T = init_Fq(p, i, 0);
        do
        {
            a = random_FpX(i, 0, p);
        } while (!FpXQ_issquare(a, T, p));
        timer_delay(&t);
        x = better_sqrt(a, T);
        timer_printf(&t, "for better sqrt in extension %ld", i);
        y = better_sqrt2(a, T);
        timer_printf(&t, "for better sqrt 2 in extension %ld", i);
        //y = FpXQ_sqrt(a, T, p);
        //timer_printf(&t, "for builtin sqrt in extension %ld", i);
        if(!gequal(a, FpXQ_sqr(x, T, p))) printf("better_sqrt not working\n");
        if(!gequal(a, FpXQ_sqr(y, T, p))) printf("better_sqrt2 not working\n");

    }
    
    
    fclose(fptr);
    pari_close();
}