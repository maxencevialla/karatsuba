#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "gmp.h"

int determineDegre(mpz_t P);
int* decomposeGrandNombre(int degre, mpz_t P);
int* karatsuba(int* P, int degP, int* Q, int degQ);

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

  karatsuba(polyA, degreA, polyB, degreB);

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
  //Fin vérif

  return polyTable;
}

int* karatsuba(int* P, int degP, int* Q, int degQ) {

  int deg;
  int sizeDeg = deg*sizeof(int);

  if(degP != degQ) {
    printf("Les degrés de A et B sont différents, on sait pas (encore) faire\n");
    return;
  } else {
    deg = degP/2;
  }

  //On découpe A et B pour Karatsuba
  int* A0 = malloc(sizeDeg);
  int* A1 = malloc(sizeDeg);
  int* B0 = malloc(sizeDeg);
  int* B1 = malloc(sizeDeg);
  int* SA = malloc(sizeDeg);
  int* SB = malloc(sizeDeg);

  for(int i = 0 ; i < deg+1 ; i++) {
    A0[i] = P[i];
    B0[i] = Q[i];
    A1[i] = P[i+deg+1];
    B1[i] = Q[i+deg+1];
    SA[i] = A0[i]+A1[i];
    SB[i] = B0[i]+B1[i];
  }

  int* res;

  if(deg == 0) {
    res = malloc(sizeof(int));
    res[0] = A0[0]*B0[0];
  } else {
    int* P2 = malloc(sizeDeg);
    int* P0 = malloc(sizeDeg);
    int* P1 = malloc(sizeDeg);
    P2 = karatsuba(A1, deg, B1, deg);
    P0 = karatsuba(A0, deg, B0, deg);
    P1 = karatsuba(SA, deg, SB, deg);

    free(A0);
    free(A1);
    free(B0);
    free(B1);
    free(SA);
    free(SB);

    for(int i = 0 ; i < deg+1 ; i++) {
      P1[i] = P1[i]-P2[i]-P0[i];
    }

    res = malloc(2*sizeDeg);

    int j;
    for(j=0 ; j<(deg/2) ; j++) {
      res[j] = P0[j];
    }
    for(j=(deg/2) ; j<deg ; j++) {
      res[j] = P1[j-deg/2];
    }
    for(j=deg ; j<2*deg ; j++) {
      res[j] = P2[j-deg];
    }

    free(P0);
    free(P1);
    free(P2);

    return res;
  }
}
