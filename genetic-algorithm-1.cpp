/*
 *  Copyright (C) 2011 Efstathios Chatzikyriakidis (stathis.chatzikyriakidis@gmail.com)
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
const uint N = 10, G = 100;

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

/* roulette wheel selection algorithm. */
const uint
roulette (const double * const f)
{
  uint i; double sum, sumrnd;

  sum = 0;
  for (i = 0; i < N; i++)
    sum += f[i] + 1;

  sumrnd = rnd () * sum;

  sum = 0;
  for (i = 0; i < N; i++) {
    sum += f[i] + 1;
    if (sum > sumrnd)
      break;
  }

  return i;
}

/* genetic algorithm entry point. */
int
main ()
{
  /* define parents and childs. */
  double parents[N], childs[N];

  /* define fitnesses for parents, childs. */
  double pFit[N], cFit[N];

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

    /* elitism process - cloning the best atom. */
    uint best = findBest (pFit);
    childs[0] = parents[best];
    cFit[0] = pFit[best];

    /* childs production process. */
    for (uint c = 1; c < N; c++) {
      /* select two parents. */
      const uint p1 = roulette (pFit);
      const uint p2 = roulette (pFit);

      /* display selected parents. */
      cout << "At. = " << setw (8) << right << p1 + 1;
      cout << " , x = " << setw (15) << parents[p1];
      cout << " , y = " << setw (15) << pFit[p1] << endl;

      cout << "At. = " << setw (8) << right << p2 + 1;
      cout << " , x = " << setw (15) << parents[p2];
      cout << " , y = " << setw (15) << pFit[p2] << endl;

      /* create child (cross-over or cloning). */
      if (rnd () < 0.9)
        childs[c] = (parents[p1] + parents[p2]) / 2;
      else
        childs[c] = rnd () < 0.5 ? parents[p1] : parents[p2];

      /* display child. */
      cout << "Ch. = " << setw (8) << right << c + 1;
      cout << " , x = " << setw (15) << childs[c];
      cout << " , y = " << setw (15) << cFit[c] << endl;

      /* child mutation. */
      if (rnd () < 0.02) {
        childs[c] += rnd () * (M_PI / 4) - (M_PI / 8);

        cout << "Child Mutation = " << childs[c] << endl;

        /* check solution space boundaries. */
        if (childs[c] < 0) childs[c] = 0;
        if (childs[c] > 2 * M_PI) childs[c] = 2 * M_PI;
      }

      /* evaluate the fitness of the child. */
      cFit[c] = fit (childs[c]);
    }

    /* display the childs with fitnesses. */
    display (childs, cFit);

    /* exchange parents with childs. */
    for (uint i = 0; i < N; i++) {
      parents[i] = childs[i];
      pFit[i] = cFit[i];
    }

    /* find the best child. */
    best = findBest (cFit);

    /* display the best child. */
    cout << "Be. = " << setw (8) << right << best + 1;
    cout << " , x = " << setw (15) << childs[best];
    cout << " , y = " << setw (15) << cFit[best] << endl;
  }

  return EXIT_SUCCESS;
}
