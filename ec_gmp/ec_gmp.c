
/*
 * Courbe elliptique sur Fp : y^2 = x^3 + a x + b (mod p)
 * GMP mpz_t pour la multiprécision.
 *
 * Compilation:
 *   gcc -o ec_gmp ec_gmp.c -lgmp
 *
 * Exécution:
 *   ./ec_gmp
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct {
    mpz_t x;
    mpz_t y;
    int infinity; // 1 => point à l'infini
} ECPoint;

typedef struct {
    mpz_t p;  // p doit être premier
    mpz_t a;
    mpz_t b;
} ECCurve;

static void modp(mpz_t r, const mpz_t x, const mpz_t p) {
    mpz_mod(r, x, p);
    if (mpz_sgn(r) < 0) mpz_add(r, r, p);
}

static void point_init(ECPoint *P) {
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1;
}

static void point_set_infinity(ECPoint *P) {
    P->infinity = 1;
    mpz_set_ui(P->x, 0);
    mpz_set_ui(P->y, 0);
}

static void point_set_xy(ECPoint *P, const mpz_t x, const mpz_t y) {
    P->infinity = 0;
    mpz_set(P->x, x);
    mpz_set(P->y, y);
}

static void point_copy(ECPoint *R, const ECPoint *P) {
    R->infinity = P->infinity;
    mpz_set(R->x, P->x);
    mpz_set(R->y, P->y);
}

static void point_print(const char *name, const ECPoint *P) {
    if (P->infinity) {
        printf("%s = O (point à l'infini)\n", name);
        return;
    }
    gmp_printf("%s = (x=%Zd, y=%Zd)\n", name, P->x, P->y);
}

static void curve_init(ECCurve *E) {
    mpz_init(E->p);
    mpz_init(E->a);
    mpz_init(E->b);
}

static inline void curve_set_str(ECCurve *E, const char *p_str, const char *a_str, const char *b_str)
{
    mpz_set_str(E->p, p_str, 10);
    mpz_set_str(E->a, a_str, 10);
    mpz_set_str(E->b, b_str, 10);
    mpz_mod(E->a, E->a, E->p);
    mpz_mod(E->b, E->b, E->p);
}

static int valid_elliptic_curve(const ECCurve *E) {
    mpz_t t1, t2, s;
    mpz_inits(t1, t2, s, NULL);
    mpz_powm_ui(t1, E->a, 3, E->p);
    mpz_mul_ui(t1, t1, 4);
    modp(t1, t1, E->p);
    mpz_powm_ui(t2, E->b, 2, E->p);
    mpz_mul_ui(t2, t2, 27);
    modp(t2, t2, E->p);
    mpz_add(s, t1, t2);
    modp(s, s, E->p);
    int valid = (mpz_cmp_ui(s, 0) != 0);
    return valid;
}

static void point_add_distinct(const ECCurve *E, ECPoint *R, const ECPoint *P, const ECPoint *Q) {

    mpz_t s, num, den, den_inv, tmp;
    mpz_inits(s, num, den, den_inv, tmp, NULL);

    // Bloc 1 : -------------------------------------
    mpz_sub(num, Q->y, P->y);
    modp(num, num, E->p);
    // ----------------------------------------------

    // Bloc 2 : -------------------------------------
    mpz_sub(den, Q->x, P->x);
    modp(den, den, E->p);
    // ----------------------------------------------
    
    // Bloc 3 : -------------------------------------
    mpz_invert(den_inv, den, E->p);
    mpz_mul(s, num, den_inv);
    modp(s, s, E->p);
    // ----------------------------------------------

    // Bloc 4 : -------------------------------------
    mpz_mul(tmp, s, s);
    modp(tmp, tmp, E->p);
    mpz_sub(tmp, tmp, P->x);
    mpz_sub(tmp, tmp, Q->x);
    modp(tmp, tmp, E->p);
    mpz_set(R->x, tmp);
    // ----------------------------------------------

    // Bloc 5 : -------------------------------------
    mpz_sub(tmp, P->x, R->x);
    modp(tmp, tmp, E->p);
    mpz_mul(tmp, s, tmp);
    modp(tmp, tmp, E->p);
    mpz_sub(tmp, tmp, P->y);
    modp(tmp, tmp, E->p);
    mpz_set(R->y, tmp);
    R->infinity = 0;
    // ----------------------------------------------
    
}

 static void point_double(const ECCurve *E, ECPoint *R, const ECPoint *P) {

     /* Cas particulier : si y == 0 alors le double est le point à l'infini */
     if (mpz_cmp_ui(P->y, 0) == 0) {
         point_set_infinity(R);
         return;
     }

     // TODO

 }

/* ============================================================
 *  Fonction principale d'addition elliptique : R = P + Q
 *      1) P = O
 *      2) Q = O
 *      3) P = -Q
 *      4) P = Q (doublement)
 *      5) Cas général (points distincts)
 * ============================================================ */

static void point_add(const ECCurve *E, ECPoint *R, const ECPoint *P, const ECPoint *Q) {

    /* 1) Si P est le point à l'infini */
    if (P->infinity) {
        point_copy(R, Q);
        return;
    }

    /* 2) Si Q est le point à l'infini */
    if (Q->infinity) {
        point_copy(R, P);
        return;
    }
    
    /* 3) Si P = -Q */
    if (mpz_cmp(P->x, Q->x) == 0) {
        mpz_t yq_neg;
        mpz_init(yq_neg);
        mpz_neg(yq_neg, Q->y);
        modp(yq_neg, yq_neg, E->p);
        if (mpz_cmp(P->y, yq_neg) == 0) {
            point_set_infinity(R);
            return;
        }
    }
    
    /* 4) Si P = Q */
    if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) == 0) {
        point_double(E, R, P);
        return;
    }

    /* 5) Cas général : points distincts */
    point_add_distinct(E, R, P, Q);
}


int main(void) {
    ECCurve E;
    curve_init(&E);
    char p_str[1024], a_str[1024], b_str[1024];
    printf("Donnez p a b (séparés par des espaces) : ");
    scanf("%s %s %s", p_str, a_str, b_str);
    curve_set_str(&E, p_str, a_str, b_str);
    ECPoint P, Q, R;
    point_init(&P);
    point_init(&Q);
    point_init(&R);
    mpz_t x, y;
    mpz_inits(x, y, NULL);
    printf("Donnez x et y de P séparés par des espaces : ");
    gmp_scanf("%Zd %Zd", x, y);
    point_set_xy(&P, x, y);
    printf("Donnez x et y de Q séparés par des espaces : ");
    gmp_scanf("%Zd %Zd", x, y);
    point_set_xy(&Q, x, y);
    point_print("P", &P);
    point_print("Q", &Q);
    point_add(&E, &R, &P, &Q);
    point_print("R = P + Q", &R);
    return 0;
}
