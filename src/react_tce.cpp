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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "random_knuth.h"
#include "error.h"
#include "modify.h"
#include "compute.h"

using namespace SPARTA_NS;
//#include <iostream>
//using namespace std;
enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION,REVERSE_EXCHANGE};   // other files

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess,z,zi,zj,iTvib,jTvib;
  int inmode,jnmode;
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int icell = ip->icell;
  double ievib = ip->evib;
  double jevib = jp->evib;
  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;
  //cout << "react called" << endl;
  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();
  zi = 0.0;
  zj = 0.0;
  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    }
    else {
       ecc = pre_etotal;
       if (pre_etotal+r->coeff[4] <= 0.0) continue; // Cover cases where coeff[1].neq.coeff[4]
       z = pre_ave_rotdof;
       if (collide->vibstyle == SMOOTH) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
       else if (collide->vibstyle == DISCRETE) {
	   if (ievib > 1e-26) {
              inmode = species[isp].nvibmode;
	      iTvib = newtonTvib(inmode,ievib,particle->species[isp].vibtemp,3000,1e-4,1000); //Instantaneous T
	      zi = (2 * ievib)/(update->boltz * iTvib); //Instantaneous z
	    }

	    if (jevib > 1e-26) {
              jnmode = species[jsp].nvibmode;
	      jTvib = newtonTvib(jnmode,jevib,particle->species[jsp].vibtemp,3000,1e-4,1000); //Instantaneous T
	      zj = (2 * jevib)/(update->boltz * jTvib); //Instantaneous z
	    }  

            //cout << zi << " " << ievib << " " << zj << " " << jevib << endl;
	    if (isnan(zi) || isnan(zj) || zi<0 || zj<0) {
              //cout << zi << " " << zj << " " << ievib << " " << jevib << endl;
              error->all(FLERR,"Root-Finding Error");
	    }
            z = pre_ave_rotdof + 0.5 * (zi+zj);
       }
       
    }

    if (temp[icell] < 1 || temp[icell] > 1000000) compute_per_grid();
    if (temp[icell] < 1 || temp[icell] > 1000000) continue;
    //if (temp[icell] < 1) cout << "grid cell temp still zero" << endl;
    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
        else e_excess = ecc + r->coeff[4];
        if (e_excess <= 0.0) continue; 
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) * 
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
        if (isnan(react_prob)) { 
           //cout << zi << " " << zj << endl;
           error->all(FLERR,"Reaction Test Error");
        }
        break;
      }

    case REVERSE_EXCHANGE:
      {
        double Kb, Qi, Qj, Qk, Ql, Qreact;
        if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
        else e_excess = ecc + r->coeff[4];
        if (e_excess <= 0.0) continue; 
        if (partition && partialEnergy==0) {
		int ksp = r->products[0]; 
		int lsp = r->products[1]; 
		Qi = Partition(species[isp].nvibmode, species[isp].rotdof, particle->species[isp].symnum, particle->species[isp].elecdegen, particle->species[isp].vibtemp, particle->species[isp].RotTemp, species[isp].mass, temp[icell]);
		Qj = Partition(species[jsp].nvibmode, species[jsp].rotdof, particle->species[jsp].symnum, particle->species[jsp].elecdegen, particle->species[jsp].vibtemp, particle->species[jsp].RotTemp, species[jsp].mass, temp[icell]);
		Qk = Partition(species[ksp].nvibmode, species[ksp].rotdof, particle->species[ksp].symnum, particle->species[ksp].elecdegen, particle->species[ksp].vibtemp, particle->species[ksp].RotTemp, species[ksp].mass, temp[icell]);
		Ql = Partition(species[lsp].nvibmode, species[lsp].rotdof, particle->species[lsp].symnum, particle->species[lsp].elecdegen, particle->species[lsp].vibtemp, particle->species[lsp].RotTemp, species[lsp].mass, temp[icell]);
		double EactivBackward = r->coeff[1] + (-1 * r->coeff[4]);
		if (EactivBackward < 0) EactivBackward = 0;
		Qreact = (Qi * Qj)/(Qk * Ql)*exp((-1) * (r->coeff[1]-EactivBackward)  / (temp[icell] * update->boltz));
		react_prob += (1/Qreact) * r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
		  pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
		  pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
		//cout << react_prob << endl;
                //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
		//cout << Qj << " " << Qk << " " << Ql << " " << Qi << " " << Qreact << " " << EactivBackward << " " << r->coeff[1]-EactivBackward << " " << temp[icell] << " " << react_prob << endl;
        }
        else {
		react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
		  pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) * 
		  pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
        }


        if (isnan(react_prob)) {   
           //cout << Qi << " " << Qj << " " << Qk << " " << Ql << endl;
           //cout << zi << " " << zj << endl;
           //cout << icell << " " << temp[icell] << " " << sizeof(temp)/sizeof(temp[0]) << endl;
           error->all(FLERR,"Reaction Test Error");
        }
        break;
      }

    case RECOMBINATION:
      {
        // skip if no 3rd particle chosen by Collide::collisions()
        //   this includes effect of boost factor to skip recomb reactions
        // check if this recomb reaction is the same one
        //   that the 3rd particle species maps to, else skip it
        // this effectively skips all recombinations reactions
        //   if selected a 3rd particle species that matches none of them
        // scale probability by boost factor to restore correct stats

        if (recomb_species < 0) continue;
        int *sp2recomb = reactions[isp][jsp].sp2recomb;
        if (sp2recomb[recomb_species] != list[i]) continue;
        if (partialEnergy) {
            react_prob += recomb_boost * recomb_density * r->coeff[2] *
              pow(ecc,r->coeff[3]) *
              pow(1.0-r->coeff[1]/ecc,r->coeff[5]);
        } else {

		double Kb;
		double diam = collide->extract(isp,jsp,"diam");
		double omega = collide->extract(isp,jsp,"omega");
		double tref = collide->extract(isp,jsp,"tref");
		static const double MY_PI  = 3.14159265358979323846; // pi
		static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
                double Qi,Qj,Qk,Qreact;

		double epsilon = 1.0;         
		if (isp == jsp) epsilon = 2.0;

		double mr = species[isp].mass * species[jsp].mass /
		     (species[isp].mass + species[jsp].mass);
		double sigma = MY_PI*diam*diam;
		if (partition) {
			int ksp = r->products[0]; //Begin modifications from Nizenkov et al.
			Qi = Partition(species[isp].nvibmode, species[isp].rotdof, particle->species[isp].symnum, 
		                    particle->species[isp].elecdegen, particle->species[isp].vibtemp, 
		                    particle->species[isp].RotTemp, species[isp].mass, temp[icell]);
			Qj = Partition(species[jsp].nvibmode, species[jsp].rotdof, particle->species[jsp].symnum,
		                    particle->species[jsp].elecdegen, particle->species[jsp].vibtemp, 
		                    particle->species[jsp].RotTemp, species[jsp].mass, temp[icell]);
			Qk = Partition(species[ksp].nvibmode, species[ksp].rotdof, particle->species[ksp].symnum,
		                    particle->species[ksp].elecdegen, particle->species[ksp].vibtemp, 
		                    particle->species[ksp].RotTemp, species[ksp].mass, temp[icell]);
			Qreact = (Qi * Qj)/Qk;
			Kb = (r->coeff[7] * pow(temp[icell],r->coeff[3])) / Qreact;
		} else {
			Kb = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
                }
		double *vi = ip->v;
		double *vj = jp->v;

		double Tcoll = (mr * (pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2))) / (update->boltz * (5 - 2*omega));  //(Nizenkov);
		double Rcoll = (2 * MY_PIS / (epsilon)) * pow(diam,2) * pow(Tcoll/tref,1-omega) * pow(2*update->boltz*tref/mr,.5);
		Rcoll *= pow(2.5-omega,1-omega) * tgamma(2.5-omega) / tgamma(3.5-(2*omega));

		react_prob += recomb_density * Kb / Rcoll; //Third order recombination     
		//cout << react_prob << endl;
                //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
		//cout << Kb << " " << Rcoll << " " << recomb_density * Kb / Rcoll << " " << random_prob << " " << Qreact << " "  << pow(temp[icell],r->coeff[3]) << " " << temp[icell] <<endl;
		if (isnan(react_prob)) { 
		   //cout << "Partition error: " << Qi << " " << Qj << " " << Qk << endl;
		   error->all(FLERR,"Reaction Test Error");
		}
        }
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    // test against random number to see if this reaction occurs
    // if it does, reset species of I,J and optional K to product species
    // J particle is destroyed in recombination reaction, set species = -1
    // K particle can be created in a dissociation or ionization reaction,
    //   set its kspecies, parent will create it
    // important NOTE:
    //   does not matter what order I,J reactants are in compared
    //     to order the reactants are listed in the reaction file
    //   for two reasons:
    //   a) list of N possible reactions above includes all reactions
    //      that I,J species are in, regardless of order
    //   b) properties of pre-reaction state are stored in precoln:
    //      computed by setup_collision()
    //      used by perform_collision() after reaction has taken place
    //      precoln only stores combined properties of I,J
    //      nothing that is I-specific or J-specific

    if (react_prob > random_prob) {
      tally_reactions[list[i]]++;

      if (!computeChemRates) {
          ip->ispecies = r->products[0];

          switch (r->type) {
          case DISSOCIATION:
          case IONIZATION:
          case EXCHANGE:
          case REVERSE_EXCHANGE:
            {
              jp->ispecies = r->products[1];
              break;
            }
          case RECOMBINATION:
            {
              // always destroy 2nd reactant species

              jp->ispecies = -1;
              break;
            }
          }

          if (r->nproduct > 2) kspecies = r->products[2];
          else kspecies = -1;

          post_etotal = pre_etotal + r->coeff[4];
      }

      return 1;
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_Evib(int nmode, double Tvib,
	                    double vibtemp[],
	                    double Evib)
{

  // COMPUTES f FOR NEWTON'S SEARCH METHOD OUTLINED IN "newtonTvib".

  double f = -Evib;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    f += (((kb*vibtemp[i])/(exp(vibtemp[i]/Tvib)-1)));
  }

  return f;

}


/* ---------------------------------------------------------------------- */

double ReactTCE::bird_dEvib(int nmode, double Tvib, double vibtemp[])
{
  
  // COMPUTES df FOR NEWTON'S SEARCH METHOD

  double df = 0.0;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    df += ((pow(vibtemp[i],2)*kb*exp(vibtemp[i]/Tvib))/(pow(Tvib,2)*pow(exp(vibtemp[i]/Tvib)-1,2)));
  }

  return df;

}

/* ---------------------------------------------------------------------- */

double ReactTCE::newtonTvib(int nmode, double Evib, double vibTemp[],
               double Tvib0,
               double tol,
               int nmax)
{


  // FUNCTION FOR CONVERTING VIBRATIONAL ENERGY TO VIBRATIONAL TEMPERATURE
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO)
  // spacing.
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol".


  double f;
  double df;
  double Tvib, Tvib_prev;
  double err;
  int i;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels.

  // f and df are computed for Newton's search
  f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error
  Tvib = Tvib0 - (f/df);
  err = fabs(Tvib-Tvib0);

  i=2;

  // Continue to search for Tvib until the error is greater than the tolerance:
  while((err >= tol) && (i <= nmax))
  {

    Tvib_prev = Tvib;

    f = bird_Evib(nmode,Tvib,vibTemp,Evib);
    df = bird_dEvib(nmode,Tvib,vibTemp);

    Tvib = Tvib_prev-(f/df);
    err = fabs(Tvib-Tvib_prev);

    i=i+1;

  }
  
  return Tvib;

}


/* ---------------------------------------------------------------------- */

double ReactTCE::Partition(int nvibmode, int nrotmode, int sym, int edegen, double VibTemp[], double RotTemp[], double mass, double CellTemp)
{

  // FUNCTION FOR CALCULATING PARTITION FUNCTION FOR MOLECULES

  double q, qtrans, qrot, qvib, qel;
  double h = 6.62607015e-34;
  double kb = 1.38064852e-23;
  static const double MY_PI  = 3.14159265358979323846; // pi

  qtrans = pow(((2*MY_PI*mass*kb*CellTemp)/pow(h,2)),1.5);
  
  if (nrotmode == 0) {
    qrot = 1.0;
  }else if (nrotmode == 2) {
    qrot = CellTemp/((double)sym*RotTemp[0]);
  } else {
    qrot = (1/(double)sym)*pow(MY_PI,.5)*pow((pow(CellTemp,3)/(RotTemp[0]*RotTemp[1]*RotTemp[2])),.5);
  }

  qvib = 1.0;
  if (nvibmode != 0) {
    for (int i = 0; i < nvibmode; i++) {
      qvib *= (1/(1-exp(-VibTemp[i]/CellTemp)));
    }
  } 

  qel = (double)edegen;
  //std::cout << "Rot Temps: " << RotTemp[0] << " " << RotTemp[1] << " " << RotTemp[2] << std::endl;
  //std::cout << qtrans << " " << qrot << " " << qvib << " " << qel << std::endl;
  q = qtrans*qrot*qvib*qel; 
  return q;

}


