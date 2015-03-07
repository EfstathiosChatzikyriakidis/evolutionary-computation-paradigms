/*
 *  Copyright (C) 2011 Efstathios Chatzikyriakidis (contact@efxa.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/* include standard C/C++ library headers. */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

/* define unsigned integer type short name. */
typedef unsigned int uint;

/* define generations number & size of colony. */
const uint N = 10, G = 20;

/* random double number generator [0-1]. */
const double
rnd ()
{
  return (double) rand () / RAND_MAX;
}

/* fitness quality calculation. */
const double
fit (const double x)
{
  return sin (x);
}

/* display colony and fitnesses. */
void
display (const double * const c, const double * const f)
{
  for (uint i = 0; i < N; i++) {
    cout << "At. = " << setw (8) << right << i + 1;
    cout << " , x = " << setw (15) << c[i];
    cout << " , y = " << setw (15) << f[i] << endl;
  }
}

/* find the best atom. */
const uint
findBest (const double * const f)
{
  double bestFit = f[0];
  uint best = 0;

  for (uint i = 1; i < N; i++)
    if (f[i] > bestFit) {
      bestFit = f[i];
      best = i;
    }

  return best;
}

/* differential evolution algorithm entry point. */
int
main (void)
{
  /* define parents. */
  double parents[N];

  /* define fitnesses for parents. */
  double pFit[N];

  /* define atom indexes basket. */
  int basket[N];

  /* set the seed for random. */
  srand (time (NULL));

  /* initialize randomly the parents. */
  for (uint i = 0; i < N; i++)
    parents[i] = rnd () * 2 * M_PI;

  /* evaluate the fitness of each parent. */
  for (uint i = 0; i < N; i++)
    pFit[i] = fit (parents[i]);

  /* display the parents with fitnesses. */
  display (parents, pFit);

  /* start evolution process. */
  for (uint g = 0; g < G; g++) {
    cout << "Generation = " << g + 1 << endl;

    /* childs production process. */
    for (uint c = 0; c < N; c++) {
      /* array used for child production. */
      int x[3];

      /* init the atom indexes basket. */
      for (uint i = 0; i < N; i++)
        basket[i] = i;

      /* put the last basket item to the current. */
      basket[c] = basket[N-1];

      /* select three totally different atoms. */
      for (uint m = 0; m < 3; m++) {
        const uint k = random () % (N-m-1);
        x[m] = basket[k];
        basket[k] = basket[N-m-2];
      }

      /* create child. */
      double vi = parents[x[0]] + 0.5 * (parents[x[1]] - parents[x[2]]);

      /* evaluate the fitness of the child. */
      double q = fit (vi);

      /* check if child is better than parent. */
      if (q > pFit[c]) {
        parents[c] = vi;
        pFit[c] = q;
      }
    }

    /* display the childs with fitnesses. */
    display (parents, pFit);

    /* find the best child. */
    uint best = findBest (pFit);

    /* display the best child. */
    cout << "Be. = " << setw (8) << right << best + 1;
    cout << " , x = " << setw (15) << parents[best];
    cout << " , y = " << setw (15) << pFit[best] << endl;
  }

  return EXIT_SUCCESS;
}
