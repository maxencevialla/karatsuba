#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "gmp.h"

int determineDegre(mpz_t P);
int* decomposeGrandNombre(int degre, mpz_t P);
void karatsuba(int* P, int degP, int* Q, int degQ);

int main(int argc, char* argv[]){

  //On récupère les 2 emp donnés en entrée et on les stocke dans a et b
  mpz_t a,b;

  mpz_init(a);
  mpz_init(b);

  mpz_set_str(a, argv[1], 10);
  mpz_set_str(b, argv[2], 10);

  //On décompose a et b en polynômes sur la base INT_MAX^i
  int degreA = determineDegre(a);
  int degreB = determineDegre(b);

  int* polyA = malloc(degreA*sizeof(int));
  int* polyB = malloc(degreB*sizeof(int));

  polyA = decomposeGrandNombre(degreA, a);
  polyB = decomposeGrandNombre(degreB, b);

  //Pour vérifier le résultat de notre multiplication
  mpz_t res;
  mpz_init(res);
  mpz_mul(res,a,b);


  printf("a=");
  mpz_out_str(stdout,10,a);
  printf("\nb=");
  mpz_out_str(stdout,10,b);
  printf("\na*b=");
  mpz_out_str(stdout,10,res);
  printf("\n");


  // mpz_t p;
  // mpz_init(p);
  //
  // mpz_out_str(stdout,10,p);
}

//Détermination du dégré du polynôme représentatif
int determineDegre(mpz_t A) {
  int degre = 0;
  mpz_t currentMax;
  mpz_init_set_ui(currentMax, INT_MAX);

  //Tant que A est plus grand que currentMax, c'est qu'on ne peut pas le représenter sous la forme d'un polynôme de degré degre.
  while(mpz_cmp(A, currentMax) > 0) {
    degre++;
    mpz_mul_ui(currentMax, currentMax, INT_MAX);
  }

  return degre;
}

//Décompose le polynôme P en base INT_MAX
int* decomposeGrandNombre(int degre, mpz_t P) {
  mpz_t tmp, currentMax,copyP;
  mpz_init(tmp);
  mpz_init(currentMax);
  mpz_init(copyP);
  mpz_set(copyP, P);

  mpz_ui_pow_ui(currentMax, INT_MAX, degre);

  int* polyTable = malloc(degre*sizeof(int));

  for(int i = degre; i > -1 ; i--) {
    //Détermination du coef de degré i
    mpz_tdiv_qr(tmp, copyP, copyP, currentMax);
    polyTable[i] = mpz_get_ui(tmp);
    mpz_mul_ui(tmp, currentMax, polyTable[i]);

    mpz_divexact_ui(currentMax, currentMax, INT_MAX);

    //printf("Degre %d, coef=%d\n", i, polyTable[i]);
  }

  //Verif
  mpz_set_ui(tmp, 0);
  mpz_t currentPow,tmp2;
  mpz_init_set_ui(currentPow, 1);
  mpz_init(tmp2);

  //On reconstruit le polynôme degré par degré
  for(int j = 0 ; j < degre+1 ; j++) {
    mpz_mul_ui(tmp2,currentPow,polyTable[j]);
    mpz_add(tmp, tmp, tmp2);
    mpz_mul_ui(currentPow, currentPow, INT_MAX);
  }
  printf("P=");
  mpz_out_str(stdout,10,tmp);
  printf("\n");

  return polyTable;
}

void karatsuba(int* P, int degP, int* Q, int degQ) {
}
