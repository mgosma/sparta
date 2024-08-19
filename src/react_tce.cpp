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
#include <iostream>
using namespace std;
enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION,REVERSE_EXCHANGE,REVERSE_DISSOCIATION,REVERSE_RECOMBINATION};   // other files

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
  double pre_etotal,ecc,e_excess,z,zi,zj,iTvib,jTvib,EReact;
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
  double Rcoll;
  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();
  zi = 0.0;
  zj = 0.0;
  iTvib = 0.0;
  jTvib = 0.0;
  
  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    //temp[icell] = 5000.0;
    //cout << temp[icell] << endl;
    r = &rlist[list[i]];

    pre_etotal = pre_etrans + pre_erot + pre_evib;
    EReact = r->coeff[4];

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    }
    else if (TCELL) {
       if (temp[icell] < 1 || temp[icell] > 10000000) compute_per_grid();
       if (temp[icell] < 1 || temp[icell] > 10000000 || isnan(temp[icell])) continue;
       if (thermo) {
         ecc = pre_etotal;
         //cout << "about to enter T loop" << endl;
         //EReact = Ereact(r, temp[icell], recomb_density);
         //e_excess = ecc + EReact;
         //if (e_excess <= 0.0) continue;
       }
       else {
         ecc = pre_etotal;
         if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
         else e_excess = ecc + r->coeff[4];
         //if (e_excess <= 0.0) continue;
       }
 
       double Kb;
       double diam = collide->extract(isp,jsp,"diam");
       double omega = collide->extract(isp,jsp,"omega");
       double tref = collide->extract(isp,jsp,"tref");
       static const double MY_PI  = 3.14159265358979323846; // pi
       static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)

       double epsilon = 1.0;         
       if (isp == jsp) epsilon = 2.0;
       
       double mr = species[isp].mass * species[jsp].mass /
            (species[isp].mass + species[jsp].mass);
       double sigma = MY_PI*diam*diam;
       double *vi = ip->v;
       double *vj = jp->v;
       double Tcoll = (mr * (pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2))) / (update->boltz * (5 - 2*omega));  //(Nizenkov);
       Rcoll = (2 * MY_PIS / (epsilon)) * pow(diam,2) * pow(Tcoll/tref,1-omega) * pow(2*update->boltz*tref/mr,.5);
       Rcoll *= pow(2.5-omega,1-omega) * tgamma(2.5-omega) / tgamma(3.5-(2*omega));

       if (isnan(Rcoll)) { 
         cout << Rcoll << " " << Tcoll << " " << vi[0] << " " << vi[1] << " " << vi[2] << " " << vj[0] << endl;
       }  
    }
    else {
       if (temp[icell] < 1 || temp[icell] > 10000000) compute_per_grid();
       if (temp[icell] < 1 || temp[icell] > 10000000 || isnan(temp[icell])) continue;
       if (thermo) {
         ecc = pre_etotal;
         //cout << (pre_etotal/(3 * update->boltz))  << endl;
         //if (!isothermal) EReact = Ereact(r, temp[icell], recomb_density);
         if (isothermal) EReact = (-1) * (2.5-0.7223) * update->boltz * (r->nreactant - r->nproduct) * temp[icell];
         else EReact = 0.0;
         e_excess = ecc + EReact;
         //cout << e_excess << endl;
         if (e_excess <= 0.0) continue;
         //cout << ecc << " " << EReact << endl;
       }
       else {
         ecc = pre_etotal;
         if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
         else e_excess = ecc + r->coeff[4];
         if (e_excess <= 0.0) continue;
       }
       if (collide->rotstyle == NONE) z = 0;
       else z = pre_ave_rotdof;
       if (collide->vibstyle == SMOOTH) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
       else if (collide->vibstyle == DISCRETE) {             
            inmode = species[isp].nvibmode;
            jnmode = species[jsp].nvibmode;

	    //Cell-Averaged z for diatomic molecules (note, this should probably be Tvib instead of Tcell)
            if (inmode == 0) zi=0.0;
            else if (inmode == 1 && temp[icell] > 300.0) zi = 2. * (1 / (exp(particle->species[isp].vibtemp[0] / temp[icell]) - 1)) * log(1.0 / (1 / (exp(particle->species[isp].vibtemp[0] / temp[icell]) - 1)) + 1.0 ); 
            else if (inmode > 1) {
	      if (ievib < 1e-26 ) zi = 0.0; //Low Energy Cut-Off to prevent nan solutions to newtonTvib
              //Instantaneous T for polyatomic
              else {
                iTvib = newtonTvib(inmode,ievib,particle->species[isp].vibtemp,3000,1e-4,1000); 
	        zi = (2 * ievib)/(update->boltz * iTvib); 
              }
	    }

            if (jnmode == 0) zj=0.0;
            if (jnmode == 1 && temp[icell] > 300.0) zj = 2. * (1 / (exp(particle->species[jsp].vibtemp[0] / temp[icell]) - 1)) * log(1.0 / (1 / (exp(particle->species[jsp].vibtemp[0] / temp[icell]) - 1)) + 1.0 ); 
            else if (jnmode > 1) {
	      if (jevib < 1e-26) zj = 0.0;
              else {
                jTvib = newtonTvib(jnmode,jevib,particle->species[jsp].vibtemp,3000,1e-4,1000);
	        zj = (2 * jevib)/(update->boltz * jTvib);  
              }         
	    }  
	    if (isnan(zi) || isnan(zj) || zi<0 || zj<0) {
              error->all(FLERR,"Root-Finding Error");
	    }
            z += 0.5 * (zi+zj);
            
       }
       //cout << z << " " << zi << " " << zj << " " << inmode << " " << jnmode << " " << ievib << " " <<jevib << endl;
    }

    //if (temp[icell] < 1) cout << "grid cell temp still zero" << endl;
    // compute probability of reaction
    //cout << react_prob << endl;
    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        //if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
        //else e_excess = ecc + r->coeff[4];
        //if (e_excess <= 0.0) continue; 
        if (TCELL) {
          double Kf = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
          
          if (r->coeff[0] == 0.0) react_prob += Kf / (Rcoll * recomb_density); //First order dissociation
	  else react_prob += Kf / (Rcoll); //Standard Second order reactions
          //cout << react_prob << " " << Kf << " " << Rcoll << " " << recomb_density << endl;         
        }
        else if (thermo) {
          double EaNew = MAX(-EReact,0.0);
          double Kf = r->coeff[2]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz));
          double n = 0.0;
          if (r->coeff[0] == 0.0) Kf *= 1 / (recomb_density); //First order dissociation
          react_prob += Kf * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+n+1.5)) *
            pow(ecc-EaNew,n-1+r->coeff[5]) * 
            pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,n-1.0+r->coeff[5]);
          //cout << r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz)) << " " << react_prob << endl;
          //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
        }
        else {
          react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
            pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) * 
            pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,r->coeff[3]-1.0+r->coeff[5]);
          //cout << z << endl;
        }
        if (isnan(react_prob)) { 
           //cout << zi << " " << zj << endl;
           error->all(FLERR,"Reaction Test Error");
        }
        break;
      }

    case REVERSE_DISSOCIATION:
      {
        double Qreact,EactivBackward;
	if (partition) {
                double Qi,Ql,Qk;
		int disssp = r->reactants[0];		
		int ksp = r->products[0]; 
		int lsp = r->products[1];
                
		Qi = Partition(species[disssp].nvibmode, species[disssp].rotdof, species[disssp].nelecmode, particle->species[disssp].symnum, 
		                    particle->species[disssp].elecdeg, particle->species[disssp].electemp, particle->species[disssp].vibtemp, 
		                    particle->species[disssp].rottemp, species[disssp].mass, temp[icell]);
		Qk = Partition(species[ksp].nvibmode, species[ksp].rotdof, species[ksp].nelecmode, particle->species[ksp].symnum, 
		                    particle->species[ksp].elecdeg, particle->species[ksp].electemp, particle->species[ksp].vibtemp, 
		                    particle->species[ksp].rottemp, species[ksp].mass, temp[icell]);
		Ql = Partition(species[lsp].nvibmode, species[lsp].rotdof, species[lsp].nelecmode, particle->species[lsp].symnum, 
		                    particle->species[lsp].elecdeg, particle->species[lsp].electemp, particle->species[lsp].vibtemp, 
		                    particle->species[lsp].rottemp, species[lsp].mass, temp[icell]);
	        EactivBackward = (-1 * r->coeff[4]);
	        Qreact = (Qi)/(Qk * Ql) * exp((EactivBackward) / (temp[icell] * update->boltz));
        }
        else if (thermo && partialEnergy==0) {
                double Gi,Gl,Gk;
                int disssp = r->reactants[0];
		int ksp = r->products[0]; 
		int lsp = r->products[1];
		Gi = Gibbs(species[disssp].thermo_T, species[disssp].thermo_Low_Poly, species[disssp].thermo_High_Poly, temp[icell], recomb_density);
		Gk = Gibbs(species[ksp].thermo_T, species[ksp].thermo_Low_Poly, species[ksp].thermo_High_Poly, temp[icell], recomb_density);
		Gl = Gibbs(species[lsp].thermo_T, species[lsp].thermo_Low_Poly, species[lsp].thermo_High_Poly, temp[icell], recomb_density);
                //cout << Gi << " " << Gl << " " << Gk << " " << temp[icell] << endl;
                Qreact = exp((Gk + Gl) - (Gi) - log(recomb_density));
        }
        else {
		error->all(FLERR,"Partition Function or thermo calculation currently required for reverse reactions.");
        }

	if (TCELL) {
	  double Kf = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
	  double Kb = Kf / Qreact;
          //cout << r->coeff[7] << " " << pow(temp[icell],r->coeff[3]) << " " << exp((-1) * (r->coeff[1]-EactivBackward)  / (temp[icell] * update->boltz)) << " " << Qreact << " " << EactivBackward << " "  << (-1) * (EactivBackward)  / (temp[icell] * update->boltz) << endl;
          //cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << endl;
          if (r->coeff[0] == 0.0) react_prob += Kb / (Rcoll * recomb_density); //First order dissociation
	  else react_prob += Kb / (Rcoll); //Second order dissociation
	  if (isnan(react_prob)) { 
	     //cout << "Partition error: " << Qi << " " << Qj << " " << Qk << " " << temp[icell] << endl;
             cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << " " << Rcoll << endl;
	     //error->all(FLERR,"Reaction Test Error");
	  }  
        }

	else if (thermo) {
          double EaNew = MAX(-EReact,0.0);
	  //double Kf = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
          double Kf = r->coeff[2]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz));
	  double Kb = Kf / Qreact;

          if (r->coeff[0] == 0.0) Kb *= 1 / (recomb_density); //First order dissociation

          double n = 0.0;
          //cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << endl;
          react_prob += Kb * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+n+1.5)) *
            pow(ecc-EaNew,n-1+r->coeff[5]) * 
            pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,n-1.0+r->coeff[5]);
          //cout << r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz)) << " " << react_prob << endl;
          //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
          //cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << " " << react_prob << " " << EaNew << endl;
	  if (isnan(react_prob)) { 
             cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << " " << react_prob << " " << EaNew << endl;
             cout << tgamma(z+2.5-r->coeff[5]) << " " << MAX(1.0e-6,tgamma(z+n+1.5)) << " " << pow(ecc-EaNew,n-1+r->coeff[5]) << " " << pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) << " " << pow(update->boltz,n-1.0+r->coeff[5]) << " " << ecc << " " << EaNew << endl;
	     error->all(FLERR,"Reaction Test Error");
	  }  
        }

        else {
	  react_prob += (1/Qreact) * r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
	    pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
	    pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,r->coeff[3]-1.0+r->coeff[5]);
        }
   
	//cout << react_prob << endl;
        //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
	//cout << Kb << " " << Rcoll << " " << recomb_density << " " << react_prob << " " << Qreact << " " << (r->coeff[7] * pow(temp[icell],r->coeff[3]) * exp((-1) * (r->coeff[1])  / (temp[icell] * update->boltz))) << " " << exp((-1) * (EactivBackward)  / (temp[icell] * update->boltz)) << " " << temp[icell] << " " << Qi << " " << Ql << " " << Qk << " " << (r->coeff[7] * pow(temp[icell],r->coeff[3])) * exp((-1) * (EactivBackward)  / (temp[icell] * update->boltz)) / Qreact << endl;
	if (isnan(react_prob)) { 
	   //cout << "Partition error: " << Qi << " " << Qj << " " << Qk << " " << temp[icell] << endl;
           //cout << Kf << " " << Qreact << " " << Kb << " " << recomb_density << " " << temp[icell] << endl;
	   error->all(FLERR,"Reaction Test Error");
	}
        
        break;
      }

    case REVERSE_EXCHANGE:
      {
        double Qreact,EactivBackward;
        if (partition && partialEnergy==0) {
                double Qi, Qj, Qk, Ql;
		int ksp = r->products[0]; 
		int lsp = r->products[1]; 
		Qi = Partition(species[isp].nvibmode, species[isp].rotdof, species[isp].nelecmode, particle->species[isp].symnum, 
		                    particle->species[isp].elecdeg, particle->species[isp].electemp, particle->species[isp].vibtemp, 
		                    particle->species[isp].rottemp, species[isp].mass, temp[icell]);
		Qj = Partition(species[jsp].nvibmode, species[jsp].rotdof, species[jsp].nelecmode, particle->species[jsp].symnum, 
		                    particle->species[jsp].elecdeg, particle->species[jsp].electemp, particle->species[jsp].vibtemp, 
		                    particle->species[jsp].rottemp, species[jsp].mass, temp[icell]);
		Qk = Partition(species[ksp].nvibmode, species[ksp].rotdof, species[ksp].nelecmode, particle->species[ksp].symnum, 
		                    particle->species[ksp].elecdeg, particle->species[ksp].electemp, particle->species[ksp].vibtemp, 
		                    particle->species[ksp].rottemp, species[ksp].mass, temp[icell]);
		Ql = Partition(species[lsp].nvibmode, species[lsp].rotdof, species[lsp].nelecmode, particle->species[lsp].symnum, 
		                    particle->species[lsp].elecdeg, particle->species[lsp].electemp, particle->species[lsp].vibtemp, 
		                    particle->species[lsp].rottemp, species[lsp].mass, temp[icell]);
	        EactivBackward = r->coeff[1] + (-1 * r->coeff[4]);
	        if (EactivBackward < 0) EactivBackward = 0;
	        Qreact = (Qi * Qj)/(Qk * Ql)*exp((-1) * (r->coeff[1]-EactivBackward)  / (temp[icell] * update->boltz));
        }
        else if (thermo && partialEnergy==0) {
                double Gi, Gj, Gk, Gl;
		int ksp = r->products[0]; 
		int lsp = r->products[1];
		Gi = Gibbs(species[isp].thermo_T, species[isp].thermo_Low_Poly, species[isp].thermo_High_Poly, temp[icell], recomb_density);
		Gj = Gibbs(species[jsp].thermo_T, species[jsp].thermo_Low_Poly, species[jsp].thermo_High_Poly, temp[icell], recomb_density);
		Gk = Gibbs(species[ksp].thermo_T, species[ksp].thermo_Low_Poly, species[ksp].thermo_High_Poly, temp[icell], recomb_density);
		Gl = Gibbs(species[lsp].thermo_T, species[lsp].thermo_Low_Poly, species[lsp].thermo_High_Poly, temp[icell], recomb_density);
                Qreact = exp((Gk + Gl) - (Gi + Gj));
                //cout << isp << " " << jsp << " " << ksp << " " << lsp << endl;
                //cout << Gi << " " << Gj << " " << Gk << " " << Gl << endl;
        }
        else {
		error->all(FLERR,"Partition Function or thermo calculation currently required for reverse reactions.");
        }
        //cout << react_prob << endl;

	if (TCELL) {
	  double Kf = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
	  double Kb = Kf / Qreact;
	  react_prob += Kb / (Rcoll); //Second order reaction  
          //cout << Kb << " " << Kf << " " << Qreact << endl;
          //double dm = (species[isp].mass + species[jsp].mass - (species[r->products[0]].mass + species[r->products[1]].mass));
          //if (dm > 1e-29) {
	  //  cout << "Reverse Exchange" << endl;
	  //  cout << ip->ispecies << " " << jp->ispecies << endl;
	  //  cout << r->reactants[0] << " " << r->reactants[1] << " -> " << r->products[0] << " " << r->products[1] << endl;
	  //  cout << (species[isp].mass + species[jsp].mass - (species[r->products[0]].mass + species[r->products[1]].mass)) << endl;
          //}
      //cout << " -> " << (species[r->products[0]].mass + species[r->products[1]].mass + species[r->products[2]].mass) << endl;
        }
        else if (thermo) {
          double EaNew = MAX(-EReact,0.0);
          double Kf = r->coeff[2]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz));
          double n = 0.0;
          double Kb = Kf / Qreact;
          react_prob += Kb * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+n+1.5)) *
            pow(ecc-EaNew,n-1+r->coeff[5]) * 
            pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,n-1.0+r->coeff[5]);
          //cout << r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz)) << " " << react_prob << endl;
          //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
        }
        else {
	  react_prob += (1/Qreact) * r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
	    pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
	    pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,r->coeff[3]-1.0+r->coeff[5]);
        }   
	//cout << react_prob << endl;
        //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
	//cout << Kb << " " << Rcoll << " " << recomb_density << " " << react_prob << " " << Qreact << " " << (r->coeff[7] * pow(temp[icell],r->coeff[3])) << " " << temp[icell] <<endl;
	if (isnan(react_prob)) { 
	   //cout << "Partition error: " << Qi << " " << Qj << " " << Qk << " " << temp[icell] << endl;
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
              pow(1.0-r->coeff[1]/ecc,r->coeff[5]) / pow(update->boltz,r->coeff[3]-1.0+r->coeff[5]);
        } 
        else {
		double Kb;
                Kb = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz));
		if (TCELL) {
			if (r->coeff[0] == 0.0) react_prob += Kb / (Rcoll); //Second Order Recombination
			else react_prob += recomb_density * Kb / Rcoll; //Third order recombination 
		} 
		else if (thermo) {
		  double EaNew = MAX(-EReact,0.0);
		  double Kf = r->coeff[2]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz));
		  double n = 0.0;
		  if (r->coeff[0] == 0.0) react_prob *= recomb_density; //Second Order Recombination
		  react_prob += Kf * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+n+1.5)) *
		    pow(ecc-EaNew,n-1+r->coeff[5]) * 
		    pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,n-1.0+r->coeff[5]);
		  //cout << r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz)) << " " << react_prob << endl;
		  //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
		}
                else {
                        //Still Currently need to use the Tcell Method for Recombination
			double diam = collide->extract(isp,jsp,"diam");
			double omega = collide->extract(isp,jsp,"omega");
			double tref = collide->extract(isp,jsp,"tref");
			static const double MY_PI  = 3.14159265358979323846; // pi
			static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
			double epsilon = 1.0;         
			if (isp == jsp) epsilon = 2.0;

			double mr = species[isp].mass * species[jsp].mass /
			     (species[isp].mass + species[jsp].mass);
			double sigma = MY_PI*diam*diam;

			double *vi = ip->v;
			double *vj = jp->v;

			double Tcoll = (mr * (pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2))) / (update->boltz * (5 - 2*omega));  //(Nizenkov);
			double Rcoll = (2 * MY_PIS / (epsilon)) * pow(diam,2) * pow(Tcoll/tref,1-omega) * pow(2*update->boltz*tref/mr,.5);
			Rcoll *= pow(2.5-omega,1-omega) * tgamma(2.5-omega) / tgamma(3.5-(2*omega));
	  
			if (r->coeff[0] == 0.0) react_prob += Kb / (Rcoll); //Second Order Recombination
			else react_prob += recomb_density * Kb / Rcoll; //Third order recombination 
                }

		if (isnan(react_prob)) { 
		   error->all(FLERR,"Reaction Test Error");
		}
        }
        break;
      }

    case REVERSE_RECOMBINATION:
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
		double Kb, Qreact;
		if (partition) {
                        double Qi,Qj,Qk;
			int ksp = r->products[0];
			Qi = Partition(species[isp].nvibmode, species[isp].rotdof, species[isp].nelecmode, particle->species[isp].symnum, 
		                    particle->species[isp].elecdeg, particle->species[isp].electemp, particle->species[isp].vibtemp, 
		                    particle->species[isp].rottemp, species[isp].mass, temp[icell]);
			Qj = Partition(species[jsp].nvibmode, species[jsp].rotdof, species[jsp].nelecmode, particle->species[jsp].symnum, 
		                    particle->species[jsp].elecdeg, particle->species[jsp].electemp, particle->species[jsp].vibtemp, 
		                    particle->species[jsp].rottemp, species[jsp].mass, temp[icell]);
			Qk = Partition(species[ksp].nvibmode, species[ksp].rotdof, species[ksp].nelecmode, particle->species[ksp].symnum, 
		                    particle->species[ksp].elecdeg, particle->species[ksp].electemp, particle->species[ksp].vibtemp, 
		                    particle->species[ksp].rottemp, species[ksp].mass, temp[icell]);
			Qreact = (Qi * Qj)/Qk;
			Kb = (r->coeff[7] * pow(temp[icell],r->coeff[3])) / Qreact;
		} 
		else if (thermo) {
                        double Gi,Gj,Gk;
			int ksp = r->products[0]; 
			Gi = Gibbs(species[isp].thermo_T, species[isp].thermo_Low_Poly, species[isp].thermo_High_Poly, temp[icell], recomb_density);
			Gj = Gibbs(species[jsp].thermo_T, species[jsp].thermo_Low_Poly, species[jsp].thermo_High_Poly, temp[icell], recomb_density);
			Gk = Gibbs(species[ksp].thermo_T, species[ksp].thermo_Low_Poly, species[ksp].thermo_High_Poly, temp[icell], recomb_density);
		        Qreact = exp((Gk) - (Gi + Gj) + log(recomb_density));
			//Kb = r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-r->coeff[1]/(temp[icell]*update->boltz)) / Qreact;
		}
                else {
			error->all(FLERR,"Partition Function or thermo calculation currently required for reverse reactions.");
                }
		if (TCELL) {
			if (r->coeff[0] == 0.0) react_prob += Kb / (Rcoll); //Second Order Recombination
			else react_prob += recomb_density * Kb / Rcoll; //Third order recombination 
		}
		else if (thermo) {
		  double EaNew = MAX(-EReact,0.0);
		  double Kf = r->coeff[2]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz));
		  double n = 0.0;
		  double Kb = Kf / Qreact;
		  if (r->coeff[0] == 0.0) react_prob *= recomb_density; //Second Order Recombination
		  react_prob += Kb * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+n+1.5)) *
		    pow(ecc-EaNew,n-1+r->coeff[5]) * 
		    pow(1.0-EaNew/ecc,z+1.5-r->coeff[5]) / pow(update->boltz,n-1.0+r->coeff[5]);
		  //cout << r->coeff[7]*pow(temp[icell],r->coeff[3])*exp(-(r->coeff[1]-EaNew)/(temp[icell]*update->boltz)) << " " << react_prob << endl;
		  //cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;
		}
                else {
                        //Still Currently need to use the Tcell Method for Recombination
			double diam = collide->extract(isp,jsp,"diam");
			double omega = collide->extract(isp,jsp,"omega");
			double tref = collide->extract(isp,jsp,"tref");
			static const double MY_PI  = 3.14159265358979323846; // pi
			static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
			double epsilon = 1.0;         
			if (isp == jsp) epsilon = 2.0;

			double mr = species[isp].mass * species[jsp].mass /
			     (species[isp].mass + species[jsp].mass);
			double sigma = MY_PI*diam*diam;

			double *vi = ip->v;
			double *vj = jp->v;

			double Tcoll = (mr * (pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2))) / (update->boltz * (5 - 2*omega));  //(Nizenkov);
			double Rcoll = (2 * MY_PIS / (epsilon)) * pow(diam,2) * pow(Tcoll/tref,1-omega) * pow(2*update->boltz*tref/mr,.5);
			Rcoll *= pow(2.5-omega,1-omega) * tgamma(2.5-omega) / tgamma(3.5-(2*omega));
	  
			if (r->coeff[0] == 0.0) react_prob += Kb / (Rcoll); //Second Order Recombination
			else react_prob += recomb_density * Kb / Rcoll; //Third order recombination 
                }

		if (isnan(react_prob)) { 
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
      //if (react_prob > 1.0) cout << "React Prob Exceeds Unity: " << react_prob << " Reactants: " << r->reactants[0] << " " << r->reactants[1] << " Products: " << r->products[0] << endl;
      //if (thermo) {
        //ecc = pre_etrans;
        //cout << "about to enter T loop" << endl;
        //EReact = Ereact(r, temp[icell], recomb_density);
        //cout << EReact << endl;
        //e_excess = ecc + EReact;
        //if (e_excess <= 0.0) {
          //cout << "energy check failed" << endl;
          //continue;
        //}
        //else cout << "pass" << endl;
      //}

      if (!computeChemRates) {
          ip->ispecies = r->products[0];

          switch (r->type) {
          case DISSOCIATION:
          case REVERSE_DISSOCIATION:
          case IONIZATION:
          case EXCHANGE:
          case REVERSE_EXCHANGE:
            {
              jp->ispecies = r->products[1];
              break;
            }
          case RECOMBINATION:
          case REVERSE_RECOMBINATION:
            {
              jp->ispecies = -1;
              break;
            }
          }

          if (r->nproduct > 2) kspecies = r->products[2];
          else kspecies = -1; 

	  //if (isothermal) post_etotal = pre_etotal;
          post_etotal = pre_etotal + EReact;
          //cout << pre_etotal << " " << post_etotal << endl;
          //cout << "react 1: " << r->reactants[0] << " " << "react 2: " << r->reactants[1] << " product 1: " << r->products[0] << " " << "product 2: " << r->products[1] << " Energy: " << ecc << " " << post_etotal << " " << EReact << " " << e_excess << endl;
          if (isnan(post_etotal)) cout << r->coeff[0] << " " << r->coeff[1] << " " << r->coeff[2] << " " << r->coeff[3] << " " << r->coeff[4] << " " << r->coeff[7] << " " << endl;;
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

double ReactTCE::Partition(int nvibmode, int nrotmode, int nelmode, int sym, int edegen[], double ElecTemp[], double VibTemp[], double RotTemp[], double mass, double CellTemp)
{

  // FUNCTION FOR CALCULATING PARTITION FUNCTION FOR MOLECULES

  double q, qtrans, qrot, qvib, qel;
  double h = 6.62607015e-34;
  double kb = 1.38064852e-23;
  static const double MY_PI  = 3.14159265358979323846; // pi

  qtrans = pow(((2*MY_PI*mass*kb*CellTemp)/pow(h,2)),1.5);

  if (nrotmode>0 && RotTemp[0]==0) {
    error->one(FLERR,"Species rotational data not read");
  }

  if (nrotmode == 0) {
    qrot = 1.0;
  } else if (nrotmode == 2) {
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

  //qel = (double)edegen;
  qel = 0.0;
  if (edegen[0]==0) {
    error->one(FLERR,"Species electronic data not read");
  }
  for (int i = 0; i < nelmode; i++) {
    if ((ElecTemp[i]) > (10*CellTemp)) break; //Check to improve computational time    
    qel += ((double)edegen[i])*exp(-(ElecTemp[i])/(CellTemp));
  }
  //std::cout << "Rot Temps: " << RotTemp[0] << " " << RotTemp[1] << " " << RotTemp[2] << std::endl;
  //std::cout << qtrans << " " << qrot << " " << qvib << " " << qel << " " <<  qtrans*qrot*qvib*qel  << std::endl;
  q = qtrans*qrot*qvib*qel; 
  return q;

}

/* ---------------------------------------------------------------------- */

double ReactTCE::Gibbs(double ThermoTemp[], double LowPoly[], double HighPoly[], double CellTemp, double Number_Density)
{

  // FUNCTION FOR CALCULATING GIBBS FREE ENERGY FOR MOLECULES

  double H, S, G;
  double Poly[7];

  if (ThermoTemp[1]==0) {
    error->one(FLERR,"Species thermo data not read");
  }

  if (ThermoTemp[1] > CellTemp) {
    for (int i = 0; i < 7; i++) {
      Poly[i] = LowPoly[i];
    }
  } else {
    for (int i = 0; i < 7; i++) {
      Poly[i] = HighPoly[i];
    }
  }
  
  H = Poly[0];
  S = Poly[0] * log(CellTemp);

  for (int i = 1; i < 5; i++) {
    H += Poly[i] * pow(CellTemp,i) / (i+1);
    S +=  Poly[i] * pow(CellTemp,i) / (i);
  }

  H += Poly[5] / CellTemp;
  S += Poly[6];

  // Pressure-Dependent Adjustment of Entropy from Ref Values (Assume Reference pressure is 1 atm)
  double Pref = 101325.0;
  double kb = 1.38064852e-23;
  double P = Number_Density * kb * CellTemp;
  S -= log(P/Pref);

  G = H - S;
  return G;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::Ereact(OneReaction *r, double CellTemp, double Number_Density)
{
  //cout << "entered loop" << endl;
  Particle::Species *species = particle->species;
  double f;
  double df;
  double Treact, Treact_prev, Ereact, Hi;
  double err;
  int i,j;

  double kb = 1.38064852e-23;
  double tol = 1e-4;
  int nmax = 1000;
  int nreactant = r->nreactant;
  int nproduct = r->nproduct;
  double Treact0 = CellTemp;
  double Greact;
  int sp;

  //Removes Third Body from Recombination reactions
  switch (r->type) {
  case DISSOCIATION:
  case REVERSE_DISSOCIATION:
  case IONIZATION:
  case EXCHANGE:
  case REVERSE_EXCHANGE:
    {
      break;
    }
  case RECOMBINATION:
  case REVERSE_RECOMBINATION:
    {
      nreactant = 3;
      break;
    }
  }

  Greact = 0.0;
  Hi = 0.0;
  //if (Greact != 0.0) cout << Greact << endl;
  // Uses Newton's method to solve for a post-reaction temp given species thermodynamic data
  // Determine enthalpy of reaction
  for (int j = 0; j < nreactant; j++) {
    if (j==2) sp = r->products[1];
    else sp = r->reactants[j];
    Greact -= HSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, CellTemp, i);
    //cout << Greact << " " << sp << " " << CellTemp << endl;
  }

  for (int j = 0; j < nproduct; j++) {
    sp = r->products[j];
    Hi += HSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, CellTemp, i);
    //cout << Hi << " " << sp << " " << CellTemp << endl;
  }

  Greact += Hi;

  // f and df are computed for Newton's search
  // f is attempting to match Greact
  f = -Greact - Hi;
  df = 0.0;
  for (int j = 0; j < nproduct; j++) {
    sp = r->products[j];
    f += HSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, Treact0, i);
    df += dHSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, Treact0);
  }

  // Update guess for Tvib and compute error
  Treact = Treact0 - (f/df);
  if (Treact < 100) Treact = 100;
  err = fabs(Treact-Treact0);
  //cout << i << " " << Treact << " " << err << " " << CellTemp << endl;
  i=1;

  // Continue to search for Tvib until the error is greater than the tolerance:
  while((err >= tol) && (i <= nmax))
  {
    Treact_prev = Treact;
    f = -Greact - Hi;
    df = 0.0;
    for (int j = 0; j < nproduct; j++) {
      sp = r->products[j];
      f += HSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, Treact, i);
      df += dHSpec(species[sp].thermo_T, species[sp].thermo_Low_Poly, species[sp].thermo_High_Poly, Treact);
    }

    Treact = Treact_prev-(f/df);
    if (Treact < 100) Treact = 100;
    //if (Treact > 6000) Treact = 6000;
    err = fabs(Treact-Treact_prev);
    //cout << i << " " << Treact << " " << err << " " << f << " " << df << endl;
    

    i=i+1;
  }
  
  Ereact = (3 / 2) * nproduct * kb * (Treact - Treact0);
  //cout << i << " " << Treact0 << " " << Treact << Greact << " " << Number_Density << " " << Ereact << " " << nreactant << " " << nproduct << " " << r->reactants[0] << " " << r->reactants[1] << endl;
  if (isnan(Treact)) cout << i << " " << Treact0 << " " << Treact << " " << Greact << " " << Number_Density << " " << Ereact << " " << nreactant << " " << nproduct << endl;
  if (isnan(Treact)) error->all(FLERR,"Reaction Test Error");
  return Ereact;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::HSpec(double ThermoTemp[], double LowPoly[], double HighPoly[], double CellTemp, int index)
{

  // FUNCTION FOR CALCULATING GIBBS FREE ENERGY FOR MOLECULES

  double H;
  double Poly[7];

  if (ThermoTemp[1]==0) {
    error->one(FLERR,"Species thermo data not read");
  }

  if (ThermoTemp[1] > CellTemp) {
    for (int i = 0; i < 7; i++) {
      Poly[i] = LowPoly[i];
    }
  } else {
    for (int i = 0; i < 7; i++) {
      Poly[i] = HighPoly[i];
      //cout << Poly[i] << endl;
    }
  }
  
  H = Poly[0];

  for (int i = 1; i < 5; i++) {
    H += Poly[i] * pow(CellTemp,i) / (i+1);
    
  }

  H += Poly[5] / CellTemp;
  //if (index == 0) cout << H << " "  << " " << CellTemp << endl;
  return H;
}

double ReactTCE::dHSpec(double ThermoTemp[], double LowPoly[], double HighPoly[], double CellTemp)
{

  // FUNCTION FOR CALCULATING GIBBS FREE ENERGY FOR MOLECULES

  double H;
  double Poly[7];

  if (ThermoTemp[1]==0) {
    error->one(FLERR,"Species thermo data not read");
  }

  if (ThermoTemp[1] > CellTemp) {
    for (int i = 0; i < 7; i++) {
      Poly[i] = LowPoly[i];
    }
  } else {
    for (int i = 0; i < 7; i++) {
      Poly[i] = HighPoly[i];
    }
  }
  
  H = 0;

  for (int i = 1; i < 5; i++) {
    H += Poly[i] * i * pow(CellTemp,i-1) / (i+1);
  }

  H -= Poly[5] / pow(CellTemp,2);
  return H;
}

double ReactTCE::dGibbs(double ThermoTemp[], double LowPoly[], double HighPoly[], double CellTemp, double Number_Density)
{

  // FUNCTION FOR CALCULATING GIBBS FREE ENERGY FOR MOLECULES

  double H, S, G;
  double Poly[7];

  if (ThermoTemp[1]==0) {
    error->one(FLERR,"Species thermo data not read");
  }

  if (ThermoTemp[1] > CellTemp) {
    for (int i = 0; i < 7; i++) {
      Poly[i] = LowPoly[i];
    }
  } else {
    for (int i = 0; i < 7; i++) {
      Poly[i] = HighPoly[i];
    }
  }
  
  H = 0;
  S = Poly[0] / CellTemp;

  for (int i = 1; i < 5; i++) {
    H += Poly[i] * i * pow(CellTemp,i-1) / (i+1);
    S +=  Poly[i] * i * pow(CellTemp,i-1) / (i);
  }

  H -= Poly[5] / pow(CellTemp,2);
  S += 0;

  // Pressure-Dependent Adjustment of Entropy from Ref Values (Assume Reference pressure is 1 atm)
  double Pref = 101325.0;
  double kb = 1.38064852e-23;
  double P = Number_Density * kb * CellTemp;
  S -= 1 / CellTemp;

  G = H - S;
  return G;
}


