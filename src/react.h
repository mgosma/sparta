/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_REACT_H
#define SPARTA_REACT_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class React : protected Pointers {
 public:
  char *style;
  int nlist;                 // # of reactions read from file

  int recombflag;            // 1 if any recombination reactions defined
  int recombflag_user;       // 0 if user has turned off recomb reactions
  int recomb_species;        // species of 3rd particle in recomb reaction
  int computeChemRates;      // 1 if only computing a TCE rate without
                             // actually doing reaction
  int partition;         // 0 using normal chem model, 1 if using partition function to compute reverse rates
  int thermo;         // 0 using normal chem model, 1 if using thermodynamic NASA7 Polynomial data to compute reverse rates
  int TCELL;              // 0 if using normal TCE model, 1 if using cell-temperature based method
  int isothermal;          // 0 if post-reaction energies are distributed as normal
  int partialEnergy;         // 0 if using rDOF model, 1 if using all energy

  int tempwhich,tempindex,invoked_per_grid;
  char *id_temp;
  class Compute *ctemp;
  class Fix *ftemp;
  double *temp;

  double recomb_density;     // num density of particles in collision grid cell
  double recomb_boost;       // rate boost param for recombination reactions
  double recomb_boost_inverse;   // inverse of boost parameter
  double recomb_Af;
  Particle::OnePart *recomb_part3;  // ptr to 3rd particle in recomb reaction

  int copy,copymode;         // 1 if class copy

  React(class SPARTA *, int, char **);
  React(class SPARTA *sparta) : Pointers(sparta) { style = NULL; random = NULL; }
  virtual ~React();
  virtual void init() {}
  virtual int recomb_exist(int, int) = 0;
  virtual void ambi_check() = 0;
  virtual int attempt(Particle::OnePart *, Particle::OnePart *,
                      double, double, double, double &, int &) = 0;
  virtual char *reactionID(int) = 0;
  virtual double extract_tally(int) = 0;

  void modify_params(int, char **);
  virtual void compute_per_grid();
  RanKnuth* get_random() { return random; }

 protected:
  int nglocal;
  class RanKnuth *random;
};

}

#endif
