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
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

/* define unsigned integer type short name. */
typedef unsigned int uint;

/* define generations number & size of colony. */
const uint N = 20, G = 200;

/* define atom genotype length in bits. */
const uint BITS = 12;

/* random double number generator [0-1]. */
const double
rnd ()
{
  return (double) rand () / RAND_MAX;
}

/* transform an integer to string. */
const string
int2str (const int g)
{
  string x;

  int mask = 1 << (BITS - 1);

  for (int i = 0; i < BITS; i++) {
    if (g & mask)
      x += '1';
    else
      x += '0';

    mask >>= 1;
  }

  return x;
}

/* transform a genotype to phenotype. */
const double
geno2pheno (const int g)
{
  return (double) g / (pow (2, BITS) - 1) * 2 * M_PI;
}

/* fitness quality calculation. */
const double
fit (const int x)
{
  return sin (geno2pheno (x));
}

/* display colony and fitnesses. */
void
display (const int * const c, const double * const f)
{
  for (uint i = 0; i < N; i++) {
    cout << "At. = " << setw (8) << right << i + 1;
    cout << " , x = " << int2str (c[i]) << " (" << setw (15) << right << geno2pheno (c[i]) << ")";
    cout << " , y = " << setw (15) << right << f[i] << endl;
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
  uint i;
  double sum, sumrnd;

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

/* one-point binary crossover. */
const int
crossOver (const int p1, const int p2)
{
  const uint point = random () % (BITS - 1) + 1;
  const int mask = (1 << (BITS - point)) - 1;
  const int child = (p1 & (~mask)) | (p2 & mask);

  cout << "Crossover Point = " << point << endl;

  return child;
}


/* deterministic mutation (probability per bit). */
const int
mutation (int c, const double p)
{
  const double e = BITS * p;
  const int A = (int) e;
  const double rest = e - A;
  const int MP = (rnd () < rest) ? A+1 : A;

  cout << "Mutation , L: " << BITS << " , P: " << p << " , E: " << e
       << " , A: " << A << " , R: " << rest << " , MP: " << MP << endl;

  for (int i = 0; i < MP; i++) {
    cout << "Before = " << int2str (c);

    const uint pnt = random () % BITS;

    int mask = 1 << pnt;

    c ^= mask;

    cout << " , Point: " << pnt;
    cout << " , After = " << int2str (c) << endl;
  }

  return c;
}

/* classic hill climber. */
void
hillClimber (int & x, double & f)
{
  int X = x; double F = f;

  const int steps[4] = { 100, 1, -100, -1 };

  for (int i = 0; i < 4; i++) {
    if ( (x + steps[i]) >= 0 &&
         (x + steps[i]) < pow (2, BITS) ) {
      x += steps[i];
      f = fit (x);

      if (f > F) {
        X = x;
        F = f;
      }
      else {
        x = X;
        f = F;
      }
    }
  }
}

/* genetic algorithm entry point. */
int
main ()
{
  /* define parents and childs. */
  int parents[N], childs[N];

  /* define fitnesses for parents, childs. */
  double pFit[N], cFit[N];

  /* set the seed for random. */
  srand (time (NULL));

  /* initialize randomly the parents. */
  for (uint i = 0; i < N; i++)
    parents[i] = random () % (1 << BITS);

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
      cout << " , x = " << int2str (parents[p1]) << " (" << setw (15) << right << geno2pheno (parents[p1]) << ")";
      cout << " , y = " << setw (15) << right << pFit[p1] << endl;

      cout << "At. = " << setw (8) << right << p2 + 1;
      cout << " , x = " << int2str (parents[p2]) << " (" << setw (15) << right << geno2pheno (parents[p2]) << ")";
      cout << " , y = " << setw (15) << right << pFit[p2] << endl;

      /* create child (cross-over or cloning). */
      if (rnd () < 0.9)
        childs[c] = crossOver (parents[p1], parents[p2]);
      else
        childs[c] = rnd () < 0.5 ? parents[p1] : parents[p2];

      /* display child. */
      cout << "Ch. = " << setw (8) << right << c + 1;
      cout << " , x = " << int2str (childs[c]) << " (" << setw (15) << right << geno2pheno (childs[c]) << ")";
      cout << " , y = " << setw (15) << right << cFit[c]<< endl;

      /* child mutation. */
      childs[c] = mutation (childs[c], 0.01);

      /* evaluate the fitness of the child. */
      cFit[c] = fit (childs[c]);
    }

    /* display the childs with fitnesses. */
    display (childs, cFit);

    /* find the best child. */
    best = findBest (cFit);

    /* perform also hill climber for for a better solution. */
    hillClimber (childs[best], cFit[best]);

    /* exchange parents with childs. */
    for (uint i = 0; i < N; i++) {
      parents[i] = childs[i];
      pFit[i] = cFit[i];
    }

    /* display the best child. */
    cout << "Be. = " << setw (8) << right << best + 1;
    cout << " , x = " << int2str (childs[best]) << " (" << setw (15) << right << geno2pheno (childs[best]) << ")";
    cout << " , y = " << setw (15) << right << cFit[best] << endl;
  }

  return EXIT_SUCCESS;
}
