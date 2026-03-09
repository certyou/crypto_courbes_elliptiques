
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
    mpz_t t1, t2, s; // declaration de variables temporaires
    mpz_inits(t1, t2, s, NULL); // initialisation des variables temporaires
    mpz_powm_ui(t1, E->a, 3, E->p); // t1 = a^3
    mpz_mul_ui(t1, t1, 4); // t1 = 4 a^3
    modp(t1, t1, E->p); // t1 = 4 a^3 mod p
    mpz_powm_ui(t2, E->b, 2, E->p); // t2 = b^2
    mpz_mul_ui(t2, t2, 27); // t2 = 27 b^2
    modp(t2, t2, E->p);// t2 = 27 b^2 mod p
    mpz_add(s, t1, t2);// s = 4 a^3 + 27 b^2 mod p
    modp(s, s, E->p); // s = 4 a^3 + 27 b^2 mod p
    int valid = (mpz_cmp_ui(s, 0) != 0); // la courbe est valide si le discriminant est différent de 0
    return valid;
}

static void point_add_distinct(const ECCurve *E, ECPoint *R, const ECPoint *P, const ECPoint *Q) {

    mpz_t s, num, den, den_inv, tmp;
    mpz_inits(s, num, den, den_inv, tmp, NULL);

    // Bloc 1 : -------------------------------------
    mpz_sub(num, Q->y, P->y); // num = y_Q - y_P
    modp(num, num, E->p); // num = (y_Q - y_P) mod p
    // ----------------------------------------------

    // Bloc 2 : -------------------------------------
    mpz_sub(den, Q->x, P->x); // den = x_Q - x_P
    modp(den, den, E->p); // den = (x_Q - x_P) mod p
    // ----------------------------------------------
    
    // Bloc 3 : -------------------------------------
    mpz_invert(den_inv, den, E->p); // den_inv = (x_Q - x_P)^-1 mod p
    mpz_mul(s, num, den_inv); // s = (y_Q - y_P) * (x_Q - x_P)^-1 mod p
    modp(s, s, E->p); // s = (y_Q - y_P) * (x_Q - x_P)^-1 mod p
    // ----------------------------------------------

    // Bloc 4 : -------------------------------------
    mpz_mul(tmp, s, s); // tmp = s^2
    modp(tmp, tmp, E->p); // tmp = s^2 mod p
    mpz_sub(tmp, tmp, P->x); // tmp = s^2 - x_P
    mpz_sub(tmp, tmp, Q->x); // tmp = s^2 - x_P - x_Q
    modp(tmp, tmp, E->p); // tmp = (s^2 - x_P - x_Q) mod p
    mpz_set(R->x, tmp); // R->x = (s^2 - x_P - x_Q) mod p
    // ----------------------------------------------

    // Bloc 5 : -------------------------------------
    mpz_sub(tmp, P->x, R->x); // tmp = x_P - x_R
    modp(tmp, tmp, E->p); // tmp = (x_P - x_R) mod p
    mpz_mul(tmp, s, tmp); // tmp = s * (x_P - x_R)
    modp(tmp, tmp, E->p); // tmp = s * (x_P - x_R) mod p
    mpz_sub(tmp, tmp, P->y); // tmp = s * (x_P - x_R) - y_P
    modp(tmp, tmp, E->p); // tmp = (s * (x_P - x_R) - y_P) mod p
    mpz_set(R->y, tmp); // R->y = (s * (x_P - x_R) - y_P) mod p
    R->infinity = 0; // R n'est pas le point à l'infini
    // ----------------------------------------------
    
}

static void point_double(const ECCurve *E, ECPoint *R, const ECPoint *P) {
    /* Cas particulier : si y == 0 alors la tangente est verticale, 
       le double est le point à l'infini */
    if (mpz_cmp_ui(P->y, 0) == 0) {
        point_set_infinity(R);
        return;
    }

    // Variables temporaires pour les calculs intermédiaires
    mpz_t s, num, den, den_inv, tmp;
    mpz_inits(s, num, den, den_inv, tmp, NULL);

    // Calcul du numérateur de s : num = (3 * x^2 + a) mod p
    mpz_mul(num, P->x, P->x);
    mpz_mul_ui(num, num, 3);
    mpz_add(num, num, E->a);
    modp(num, num, E->p);

    // Calcul du dénominateur de s : den = (2 * y) mod p
    mpz_mul_ui(den, P->y, 2);
    modp(den, den, E->p);

    // Calcul de la pente s : num * den^-1 mod p
    if (mpz_invert(den_inv, den, E->p) == 0) {
        // Si l'inverse n'existe pas, on retourne le point à l'infini
        point_set_infinity(R);
    } else {
        mpz_mul(s, num, den_inv);
        modp(s, s, E->p);

        // Calcul de x_R : (s^2 - 2*x_P) mod p
        mpz_mul(tmp, s, s);
        mpz_submul_ui(tmp, P->x, 2);
        modp(tmp, tmp, E->p);
        mpz_set(R->x, tmp);

        // Calcul de y_R : (s * (x_P - x_R) - y_P) mod p
        mpz_sub(tmp, P->x, R->x);
        mpz_mul(tmp, s, tmp);
        mpz_sub(tmp, tmp, P->y);
        modp(tmp, tmp, E->p);
        mpz_set(R->y, tmp);
        
        R->infinity = 0;
    }
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

static void point_mul(const ECCurve *E, ECPoint *R, const int k, const ECPoint *P) {
    ECPoint res;
    point_init(&res);
    
    int i;
    for (i = 0; i < k; i++) {
        point_add(E, &res, &res, P);
    }
    point_copy(R, &res);
}

int main() {
    ECCurve E;
    mpz_inits(E.p, E.a, E.b, NULL);
    
    printf("Curve (p a b): ");
    gmp_scanf("%Zd %Zd %Zd", E.p, E.a, E.b);

    ECPoint P, Result;
    point_init(&P);
    point_init(&Result);
    int k;
    printf("multiplier k: ");
    scanf("%d", &k);

    printf("Point P (x y): ");
    gmp_scanf("%Zd %Zd", P.x, P.y);
    P.infinity = 0;
    point_mul(&E, &Result, k, &P);
    gmp_printf("Result kP: (%Zd, %Zd)\n", Result.x, Result.y);

    return 0;
}

/*
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
*/