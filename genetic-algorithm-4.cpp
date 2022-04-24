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
const uint N = 10, G = 100;

/* define atom chromosome length in bits. */
const uint BITS = 16;

/* define number of atom chromosomes. */
const uint VARS = 2;

/* define fitness function offset. */
const int OFFSET = 1;

/* define atom genotype structure. */
typedef struct atom
{
  int geno[VARS];
  double pheno[VARS];
  double fitness;
} atom;

/* random double number generator [0-1]. */
const double
rnd ()
{
  return (double) rand () / RAND_MAX;
}

/* transform a genotype to phenotype. */
void
geno2pheno (atom & g)
{
  for (uint i = 0; i < VARS; i++)
    g.pheno[i] = (double) g.geno[i] / (pow (2, BITS) - 1) * 2 * M_PI;
}

/* fitness quality calculation. */
void
fit (atom & x)
{
  geno2pheno (x);

  const double x1 = x.pheno[0];
  const double x2 = x.pheno[1];

  x.fitness = sin (x1) * sin (x2);
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

/* find the best atom. */
const uint
findBest (const atom * const f)
{
  double bestFit = f[0].fitness;
  uint best = 0;

  for (uint i = 1; i < N; i++)
    if (f[i].fitness > bestFit) {
      bestFit = f[i].fitness;
      best = i;
    }

  return best;
}

/* roulette wheel selection algorithm. */
const uint
roulette (const atom * const f)
{
  uint i;
  double sum, sumrnd;

  sum = 0;
  for (i = 0; i < N; i++)
    sum += f[i].fitness + OFFSET;

  sumrnd = rnd () * sum;

  sum = 0;
  for (i = 0; i < N; i++) {
    sum += f[i].fitness + OFFSET;
    if (sum > sumrnd)
      break;
  }

  return i;
}

/* classic hill climber. */
void
hillClimber (atom & x)
{
  int X; double F;

  const int steps[4] = { 100, 1, -100, -1 };

  for (uint j = 0; j < VARS; j++) {
    X = x.geno[j];
    F = x.fitness;

    for (uint i = 0; i < 4; i++) {
      if ( (x.geno[j] + steps[i]) >= 0 &&
           (x.geno[j] + steps[i]) < pow (2, BITS) ) {
        x.geno[j] += steps[i];

        fit (x);

        if (x.fitness > F) {
          X = x.geno[j];
          F = x.fitness;
        }
        else {
          x.geno[j] = X;
          x.fitness = F;
        }
      }
    }
  }
}

/* one-point binary crossover. */
const atom
crossOver (const atom p1, const atom p2)
{
  atom child;

  const int point = random () % (BITS * VARS - 1) + 1;
  const int chrome = point / BITS;

  for (uint i = 0; i < chrome; i++)
    child.geno[i] = p1.geno[i];

  for (uint i = chrome + 1; i < VARS; i++)
    child.geno[i] = p2.geno[i];

  const int offset = point % BITS;
  const int mask = (1 << (BITS - offset)) - 1;

  child.geno[chrome] = ((p1.geno[chrome] & (~mask)) | (p2.geno[chrome] & mask));

  cout << "Crossover: Point = " << point << " , Chrome = " << chrome << " , Offset = " << offset << endl;

  return child;
}

/* copy a chromosome of a genotype. */
void
copyChrome (atom & x)
{
  int i = random () % VARS;
  int j = random () % VARS;
  x.geno[i] = x.geno[j];
}

/* deterministic mutation (probability per bit). */
void
mutation (atom & c, const double p)
{
  const double e = BITS * p;
  const int A = (int) e;
  const double rest = e - A;
  const int MP = (rnd () < rest) ? A+1 : A;

  cout << "Mutation , L: " << BITS*VARS << " , P: " << p << " , E: " << e
       << " , A: " << A << " , R: " << rest << " , MP: " << MP << endl;

  for (uint i = 0; i < MP; i++) {
    const int pnt = random () % (BITS * VARS);
    const int chrome = pnt / BITS;
    const int offset = pnt % BITS;
    const int mask = 1 << offset;

    cout << "Before = ";
    for (uint j = 0; j < VARS; j++)
      cout << int2str (c.geno[j]) << " ";

    cout << ", Point: " << pnt;

    c.geno[chrome] ^= mask;

    cout << " , After = ";
    for (uint j = 0; j < VARS; j++)
      cout << int2str (c.geno[j]) << " ";

    cout << endl;
  }
}

/* print atom information. */
void
printAtom (const atom x, const string what, const uint index)
{
  cout << what << setw (8) << right << index + 1 << " , ";

  cout << "g = ";
  for (uint j = 0; j < VARS; j++)
    cout << int2str (x.geno[j]) << " ";

  cout << ", p = ";
  for (uint j = 0; j < VARS; j++)
    cout << setw (15) << x.pheno[j] << " ";

  cout << ", f = " << setw (15) << x.fitness << endl;
}

/* display colony and fitnesses. */
void
display (const atom * const c)
{
  for (uint i = 0; i < N; i++)
    printAtom (c[i], "At. = ", i);
}

/* genetic algorithm entry point. */
int
main (void)
{
  /* define parents and childs. */
  atom parents[N], childs[N];

  /* set the seed for random. */
  srand (time (NULL));

  /* initialize randomly the parents. */
  for (uint i = 0; i < N; i++)
    for (uint j = 0; j < VARS; j++)
      parents[i].geno[j] = random () % (1 << BITS);

  /* evaluate the fitness of each parent. */
  for (uint i = 0; i < N; i++)
    fit (parents[i]);

  /* display the parents with fitnesses. */
  display (parents);

  /* start evolution process. */
  for (uint g = 0; g < G; g++) {
    cout << "Generation = " << g + 1 << endl;

    /* elitism process - cloning the best atom. */
    uint best = findBest (parents);
    childs[0] = parents[best];

    /* childs production process. */
    for (uint c = 1; c < N; c++) {
      /* select two parents. */
      const uint p1 = roulette (parents);
      const uint p2 = roulette (parents);

      /* display selected parents. */
      printAtom (parents[p1], "At. = ", p1);
      printAtom (parents[p2], "At. = ", p2);

      /* create child (cross-over or cloning). */
      if (rnd () < 0.9)
        childs[c] = crossOver (parents[p1], parents[p2]);
      else
        childs[c] = rnd () < 0.5 ? parents[p1] : parents[p2];

      /* display child. */
      printAtom (childs[c], "Ch. = ", c);

      /* child mutation. */
      mutation (childs[c], 0.01);

      /* copy chrome. */
      if (rnd () < 0.1)
        copyChrome (childs[c]);

      /* evaluate the fitness of the child. */
      fit (childs[c]);
    }

    /* display the childs with fitnesses. */
    display (childs);

    /* find the best child. */
    best = findBest (parents);

    /* perform also hill climber for for a better solution. */
    hillClimber (childs[best]);

    /* exchange parents with childs. */
    for (int i = 0; i < N; i++)
      parents[i] = childs[i];

    /* display the best child. */
    printAtom (childs[best], "Be. = ", best);
  }

  return EXIT_SUCCESS;
}
