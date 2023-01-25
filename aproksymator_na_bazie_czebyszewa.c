#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(int k, double x)
{
	if (k == 0)
		return 1;
	else if (k == 1)
		return x;
	else
		return 2 * x * fi(k-1, x) - fi(k-2, x); 
}

/* Pierwsza pochodna fi */
double
dfi(int k, double x)
{
	if (k == 0)
		return 0;
	if (k == 1)
		return 1;
	else 
		return 2*fi(k-1, x) + 2*x*dfi(k-1, x) - dfi(k-2, x);
}

/* Druga pochodna fi */
double
d2fi(int k, double x)
{
	if (k == 0)
		return 0;
	if (k == 1)
		return 0;
	else 
		return 2 * x * d2fi(k-1,x) + 4 * dfi(k-1, x) - d2fi(k-2, x);
}

/* Trzecia pochodna fi */
double
d3fi(int k, double x)
{
	if (k == 0)
		return 0;
	if (k == 1)
		return 0;
	else 
		return 2*x*d3fi(k-1, x) + 6*d2fi(k-1, x) - d3fi(k-2, x);
}


void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (k, xx);
				spl->f1[i] += ck * dfi (k, xx);
				spl->f2[i] += ck * d2fi(k, xx);
				spl->f3[i] += ck * d3fi(k, xx);
			}
		}
	}
}
