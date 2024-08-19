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
#include "collide_vss.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "fix_vibmode.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <iostream>
using namespace std;
using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(SPARTA *sparta, int narg, char **arg) :
  Collide(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal collide command");

  // optional args

  relaxflag = CONSTANT;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"relax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
      if (strcmp(arg[iarg+1],"constant") == 0) relaxflag = CONSTANT;
      else if (strcmp(arg[iarg+1],"variable") == 0) relaxflag = VARIABLE;
      else error->all(FLERR,"Illegal collide command");
      iarg += 2;
    } else error->all(FLERR,"Illegal collide command");
  }

  // proc 0 reads file to extract params for current species
  // broadcasts params to all procs

  nparams = particle->nspecies;
  if (nparams == 0)
    error->all(FLERR,"Cannot use collide command with no species defined");

  memory->create(params,nparams,nparams,"collide:params");
  if (comm->me == 0) read_param_file(arg[2]);
  MPI_Bcast(params[0],nparams*nparams*sizeof(Params),MPI_BYTE,0,world);

  // allocate per-species prefactor array

  memory->create(prefactor,nparams,nparams,"collide:prefactor");
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  if (copymode) return;

  memory->destroy(params);
  memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  Collide::init();
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()
------------------------------------------------------------------------- */

double CollideVSS::vremax_init(int igroup, int jgroup)
{
  // parent has set mixture ptr

  Particle::Species *species = particle->species;
  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double cxs = params[isp][jsp].diam*params[isp][jsp].diam*MY_PI;
      prefactor[isp][jsp] = cxs * pow(2.0*update->boltz*params[isp][jsp].tref/
	params[isp][jsp].mr,params[isp][jsp].omega-0.5) /
	tgamma(2.5-params[isp][jsp].omega);
      double beta = MAX(vscale[isp],vscale[jsp]);
      double vrm = 2.0 * cxs * beta;
      vrmgroup = MAX(vrmgroup,vrm);
    }
  }

  return vrmgroup;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int np, double volume)
{
  double fnum = update->fnum;
  double dt = update->dt;

  double nattempt;

  if (remainflag) {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
    remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
  } else {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + random->uniform();
  }

  return nattempt;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int igroup, int jgroup,
				     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 // return 2x the value for igroup != jgroup, since no J,I pairing

 double npairs;
 if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 else npairs = ngroup[igroup] * (ngroup[jgroup]);
 //else npairs = 0.5 * ngroup[igroup] * (ngroup[jgroup]);

 nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

 if (remainflag) {
   nattempt += remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else nattempt += random->uniform();

 return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
			       Particle::OnePart *ip, Particle::OnePart *jp)
{
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  //cout << ispecies << endl;
  //cout << jspecies << endl;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double vro  = pow(vr2,1.0-params[ispecies][jspecies].omega);

  // although the vremax is calculated for the group,
  // the individual collisions calculated species dependent vre
  //if (vro == 0) return 0;
  double vre = vro*prefactor[ispecies][jspecies];
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if ((vre/vremax[icell][igroup][jgroup]) > .05) {
    //cout << vre << " " << vro << " " << ispecies << " " << jspecies << endl;
    //cout << vre/vremax[icell][igroup][jgroup] << endl;
  }
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;


/*
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  int nlocal = particle->nlocal;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  //double omega1 = params[ispecies].omega;
  //double omega2 = params[jspecies].omega;
  //double omega = 0.5 * (omega1+omega2);
  double omega = CSomega[ispecies][jspecies];
  //cout << omega << endl;
  double vro  = pow(vr2,1.0-omega);

  // although the vremax is calcualted for the group,
  // the individual collisions calculated species dependent vre
  //int kspecies, reactflag;
  //precoln.etrans = 0.5 * precoln.mr * precoln.vr2;
  //precoln.erot = ip->erot + jp->erot;
  //ievib = ip->erot;
  //jevib = jp->erot;
  //precoln.evib = ip->evib + jp->evib;
  //if (react) 
   // reactflag = react->attempt(ip,jp,
    //                           precoln.etrans,precoln.erot,
    //                           precoln.evib,postcoln.etotal,kspecies);
  
  double epsilon = 2.0;         
  double mr = 1.67356E-27 * 1.67356E-27 / (1.67356E-27 + 1.67356E-27);
  double Tcoll = (mr * (pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2))) / (update->boltz * (5 - 2*omega));  //Paper Instructions (Nizenkov);
  double Rcoll = (2 * MY_PIS / (epsilon)) * pow(2.97189e-10,2) * pow(Tcoll/273,1-omega) * pow(2*update->boltz*273/mr,.5);
  Rcoll *= pow(2.5-omega,1-omega) * tgamma(2.5-omega) / tgamma(3.5-(2*omega));
  double sigmaR = 1.66058E-15 / Rcoll;
 
 
  //cout << ispecies << " " << jspecies << endl;
  double vre = vro*prefactor[ispecies][jspecies];
 // cout << react->prob << endl;
  //if (react->prob > 1) vre = react->prob*vro*prefactor[ispecies][jspecies]; //In input file, species 2 is H; Added to correct for hydrogen recomb underprediction 
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;
*/
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  double imass = precoln.imass = species[isp].mass;
  double jmass = precoln.jmass = species[jsp].mass;

  precoln.etrans = 0.5 * params[isp][jsp].mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (imass+jmass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((imass*vi[0])+(jmass*vj[0])) * divisor;
  precoln.vcmf = ((imass*vi[1])+(jmass*vj[1])) * divisor;
  precoln.wcmf = ((imass*vi[2])+(jmass*vj[2])) * divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

int CollideVSS::perform_collision(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp)
{
  int reactflag,kspecies;
  double x[3],v[3];
  Particle::OnePart *p3;

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1
  //cout << ip->ispecies << " " << jp->ispecies << endl;
  if (react) {
    //cout << "Performing reaction" << endl;
    reactflag = react->attempt(ip,jp,
                               precoln.etrans,precoln.erot,
                               precoln.evib,postcoln.etotal,kspecies);
    //cout << "Reaction complete" << endl;
  }
  
  else reactflag = 0;

  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  if (reactflag && !react->computeChemRates) {
    //if (isnan(postcoln.etotal)) cout << "nan occurs in react loop" << endl;
    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct

    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();
      
      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int reallocflag =
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        ip = particle->particles + (ip - particles);
        jp = particle->particles + (jp - particles);
      }

      kp = &particle->particles[particle->nlocal-1];
      //cout << "Diss energy exchange" << endl;
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      //cout << "Diss exchange complete" << endl;
      SCATTER_ThreeBodyScattering(ip,jp,kp);

    // remove 2nd J particle if recombination reaction removed it
    // p3 is 3rd particle participating in energy exchange

    } else if (jp->ispecies < 0) {
      double *vi = ip->v;
      double *vj = jp->v;

      double divisor = 1.0 / (precoln.imass + precoln.jmass);
      double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
      double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

      vi[0] = ucmf;
      vi[1] = vcmf;
      vi[2] = wcmf;

      jp = NULL;
      p3 = react->recomb_part3;
      //cout << "React third body" << p3->ispecies << endl; 
      // properly account for 3rd body energy with another call to setup_collision()
      // it needs relative velocity of recombined species and 3rd body

      double *vp3 = p3->v;
      double du  = vi[0] - vp3[0];
      double dv  = vi[1] - vp3[1];
      double dw  = vi[2] - vp3[2];
      double vr2 = du*du + dv*dv + dw*dw;
      precoln.vr2 = vr2;

      // internal energy of ip particle is already included
      //   in postcoln.etotal returned from react->attempt()
      // but still need to add 3rd body internal energy

      double partial_energy =  postcoln.etotal + p3->erot + p3->evib;

      ip->erot = 0;
      ip->evib = 0;
      p3->erot = 0;
      p3->evib = 0;
	
      // returned postcoln.etotal will increment only the
      //   relative translational energy between recombined species and 3rd body
      // add back partial_energy to get full total energy

      setup_collision(ip,p3);
      postcoln.etotal += partial_energy;

      if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp);
      SCATTER_TwoBodyScattering(ip,p3);

    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }

  } else {
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp);
    SCATTER_TwoBodyScattering(ip,jp);
  }

  return reactflag;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip,
					   Particle::OnePart *jp)
{
  double ua,vb,wc;
  double vrc[3];

  Particle::Species *species = particle->species;
  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;

  double alpha_r = 1.0 / params[isp][jsp].alpha;

  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2.0 * postcoln.etrans / params[isp][jsp].mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (params[isp][jsp].mr * precoln.vr2));
    double cosX = 2.0*pow(random->uniform(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.0e-6) {
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] );
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}

/* ---------------------------------------------------------------------- */
/*
void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip,
						Particle::OnePart *jp)
{

  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib,irot;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip;
      else p = jp;

      int sp = p->ispecies;
      rotdof = species[sp].rotdof;
      double rotn_phi = species[sp].rotrel;

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= random->uniform()) {
          if (rotstyle == NONE) {
            p->erot = 0.0;
          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot =
              1- pow(random->uniform(),
		     (1/(2.5-params[ip->ispecies][jp->ispecies].omega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose *
              sample_bl(random,0.5*species[sp].rotdof-1.0,
                        1.5-params[ip->ispecies][jp->ispecies].omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;

      vibdof = species[sp].vibdof;
      double vibn_phi = species[sp].vibrel[0];

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= random->uniform()) {
          if (vibstyle == NONE) {
            p->evib = 0.0;

          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              Fraction_Vib =
                1.0 - pow(random->uniform(),
			  (1.0/(2.5-params[ip->ispecies][jp->ispecies].omega)));
              p->evib= Fraction_Vib * E_Dispose;
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int>
                  (random->uniform()*(max_level+AdjustFactor));
                p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - params[ip->ispecies][jp->ispecies].omega));
              } while (State_prob < random->uniform());
              E_Dispose -= p->evib;
            }

          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              p->evib = E_Dispose *
                sample_bl(random,0.5*species[sp].vibdof-1.0,
                          1.5-params[ip->ispecies][jp->ispecies].omega);
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              p->evib = 0.0;

              int nmode = particle->species[sp].nvibmode;
              int **vibmode =
                particle->eiarray[particle->ewhich[index_vibmode]];
              int pindex = p - particle->particles;

              for (int imode = 0; imode < nmode; imode++) {
                ivib = vibmode[pindex][imode];
                E_Dispose += ivib * update->boltz *
                  particle->species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));

                do {
                  ivib = static_cast<int>
                    (random->uniform()*(max_level+AdjustFactor));
                  pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - params[ip->ispecies][jp->ispecies].omega));
                } while (State_prob < random->uniform());

                vibmode[pindex][imode] = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
            }
          } // end of vibstyle/vibdof if
        }
        postcoln.evib += p->evib;
      } // end of vibdof if
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}
*/

/* ---------------------------------------------------------------------- */

// Prohibits Double Relaxation
void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{

  double State_prob,Fraction_Vib,E_Dispose;
  double phi,factor,transdof,pevib,Tt,vibdof_dis,A,omega;
  int i,sp,rotdof,vibdof,max_level,ivib,imode,relaxflag1;

  Particle::OnePart *p1,*p2,*p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  relaxflag1 = 0;
  phi = 0.0;
  factor= 1.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    do {
      if (0.5 < random->uniform()) {
          p1 = ip;
          p2 = jp;
      }   else {
          p1 = jp;
          p2 = ip;
      }

      for (i = 0; i < 2; i++) {
        if (i == 0) p = p1;
        else if (i == 1) p = p2;

          sp = p->ispecies;
          vibdof = species[sp].vibdof;
          //transdof = 5.0-2.0*params[sp].omega;
          omega = params[ip->ispecies][jp->ispecies].omega;
          transdof = 5.0-2.0*omega;
          Tt = (E_Dispose+p->evib) / (update->boltz * (transdof+vibdof));

          if (vibdof) {
              if (vibstyle == NONE) {
                p->evib = 0.0;

              } else if (vibstyle == SMOOTH) {

                  factor *= 1/(1-phi);
                  if (relaxflag == VARIABLE) phi = factor*(1.0 + vibdof/transdof)/species[sp].vibrel[0];
//phi = factor*(1.0 + vibdof/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof+transdof);
                  //else phi = factor*(1.0 + vibdof/transdof)/species[sp].vibrel[0];
                  else phi = 1.0/species[sp].vibrel[0];
                  if (phi >= random->uniform()) {
                      E_Dispose += p->evib;
                      if (vibdof == 2) {
                          Fraction_Vib = 1.0 - pow(random->uniform(),(1.0/(2.5-omega)));
                          p->evib= Fraction_Vib * E_Dispose;
                      }
                      else if (vibdof > 2) p->evib = E_Dispose * sample_bl(random,0.5*vibdof-1.0,1.5-omega);
                      E_Dispose -= p->evib;
                      postcoln.evib += p->evib;
                      relaxflag1 = 1;
                      break;
                  }

              } else if (vibstyle == DISCRETE) {
                  if (vibdof == 2) {

                      vibdof_dis = 2.0*(species[sp].vibtemp[0]/Tt) / (exp(species[sp].vibtemp[0]/Tt)-1);
                      A = pow(vibdof_dis,2)*exp(species[sp].vibtemp[0]/Tt)/2.0;

                      factor *= 1/(1-phi);
                      if (relaxflag == VARIABLE) phi = factor*(1.0 + A/transdof)/species[sp].vibrel[0];
//phi = factor*(1.0 + A/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof_dis+transdof);
                      //else phi = factor*(1.0 + A/transdof)/species[sp].vibrel[0];
                      else phi = 1.0/species[sp].vibrel[0];
                      if (phi >= random->uniform()) {
                          int **vibmode = particle->eiarray[particle->ewhich[index_vibmode]];
                          int pindex = p - particle->particles;
                          E_Dispose += p->evib;
                          max_level = static_cast<int>
                            (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
                          do {
                            ivib = static_cast<int>
                              (random->uniform()*(max_level+AdjustFactor));
                            p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                            State_prob = pow((1.0 - p->evib / E_Dispose),
                                             (1.5 - omega));
                          } while (State_prob < random->uniform());
                          E_Dispose -= p->evib;
                          postcoln.evib += p->evib;
                          //vibmode[pindex][0] = ivib;
                          relaxflag1 = 1;
                          break;
                      }

                  } else if (vibdof > 2) {

                      int nmode = particle->species[sp].nvibmode;
                      int **vibmode =
                        particle->eiarray[particle->ewhich[index_vibmode]];
                      int pindex = p - particle->particles;

                      imode = 0;
                      while (imode < nmode) {

                        vibdof_dis = 2.0*(species[sp].vibtemp[imode]/Tt) / (exp(species[sp].vibtemp[imode]/Tt)-1);
                        A = pow(vibdof_dis,2)*exp(species[sp].vibtemp[imode]/Tt)/2.0;

                        factor *= 1/(1-phi);
                        if (relaxflag == VARIABLE) phi = factor*(1.0 + A/transdof)/species[sp].vibrel[imode];
//phi = factor*(1.0 + A/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof_dis+transdof);
                        //else phi = factor*(1.0 + A/transdof)/species[sp].vibrel[imode];
                        else phi = 1.0/species[sp].vibrel[0];
                        if (phi >= random->uniform()) {
                            ivib = vibmode[pindex][imode];
                            E_Dispose += ivib * update->boltz *
                                    particle->species[sp].vibtemp[imode];
                            p->evib -= ivib * update->boltz *
                                    particle->species[sp].vibtemp[imode];
                            max_level = static_cast<int>
                              (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));

                            do {
                              ivib = static_cast<int>
                                (random->uniform()*(max_level+AdjustFactor));
                              pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                              State_prob = pow((1.0 - pevib / E_Dispose),
                                               (1.5 - omega));
                            } while (State_prob < random->uniform());

                            vibmode[pindex][imode] = ivib;
                            p->evib += pevib;
                            E_Dispose -= pevib;
                            postcoln.evib += pevib;
                            relaxflag1 = 1;
                            break;
                        }
                        imode++;
                      }
                      if (relaxflag1) break;
                    }
                } // end of vibstyle if
            } // end of vibdof if
          if (relaxflag1) break;
      }

      if (relaxflag1) break;

      for (i = 0; i < 2; i++) {
        if (i == 0) p = p1;
        else if (i == 1) p = p2;

          sp = p->ispecies;
          rotdof = species[sp].rotdof;
          transdof = 5.0-2.0*omega;

          if (rotdof) {
            if (rotstyle == NONE) {
                p->erot = 0.0;
            } else {
                factor *= 1/(1-phi);
                if (relaxflag == VARIABLE)  phi = factor*(1.0 + rotdof/transdof)/species[sp].rotrel;
//phi = factor*(1.0 + rotdof/transdof)/rotrel_prohibdouble(sp,E_Dispose+p->erot);
                //else phi = factor*(1.0 + rotdof/transdof)/species[sp].rotrel;
                else phi = 1.0/species[sp].rotrel;
                if (phi >= random->uniform()) {
                    E_Dispose += p->erot;
                    if (rotdof == 2) {
                        Fraction_Vib = 1.0 - pow(random->uniform(),(1.0/(2.5-omega)));
                        p->erot= Fraction_Vib * E_Dispose;
                    }
                    else if (rotdof > 2) p->erot = E_Dispose * sample_bl(random,0.5*rotdof-1.0,1.5-omega);
                    E_Dispose -= p->erot;
                    postcoln.erot += p->erot;
                    break;
                }
            }
          }
      }
    } while (0 > 1);
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_ThreeBodyScattering(Particle::OnePart *ip,
			  		     Particle::OnePart *jp,
			  		     Particle::OnePart *kp)
{
  double vrc[3],ua,vb,wc;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;
  double mass_k = species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha_r = 1.0 / params[isp][jsp].alpha;
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(random->uniform(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) {
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0];
      vb = scale * sinX*vrc[0]*cos(eps);
      wc = scale * sinX*vrc[0]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  vi[0] = precoln.ucmf + (mass_k*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_k*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_k*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_ij*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_ij*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_ij*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip,
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib,irot;
  double aveomega,pevib;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    //cout << "1 " << ip->ispecies << endl;
    //cout << jp->ispecies << endl;
    //cout << ip->erot << " " << jp->erot << endl;
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega =  params[ip->ispecies][jp->ispecies].omega; 
    //cout << "2 species" << endl;
  } else {
    //cout << "1 " << ip->ispecies << endl;
    //cout << kp->ispecies << endl;
    //cout << jp->ispecies << endl;
    //cout << ip->erot << " " << jp->erot << " " << kp->erot << endl;
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (1/3)*(params[ip->ispecies][jp->ispecies].omega + params[ip->ispecies][kp->ispecies].omega + params[jp->ispecies][kp->ispecies].omega); 
    //cout << "3 species" << endl;
  }

  // handle each kind of energy disposal for reacting reactants
  // clean up memory for the products
  
  double E_Dispose = postcoln.etotal;
  double E_PreLB = postcoln.etotal;
  //cout << E_Dispose << endl;
  double zc, Tcol, zki;
  if (E_Dispose < 0)
    error->all(FLERR,"Impossible Energy Post-collision");

  if (numspecies == 3) { 
    int isp = ip->ispecies;
    int jsp = jp->ispecies;
    int ksp = kp->ispecies;
    double aveomega12 = params[ip->ispecies][jp->ispecies].omega; 
    double aveomega123 = aveomega;
    int nvibmode[] = {species[isp].nvibmode, species[jsp].nvibmode, species[ksp].nvibmode};

    
    //Need to calculate a collisional tempertature Tcol between the particles; in addition, need to find post-collisional vibrational degrees of 
    //freedom for colliding particles; values are then used for computation of the vibrational energies, assuming each vib mode is a harmonic oscillator
    double omega[] = {aveomega12, aveomega123};
    double nrotmode[] = {(double)species[isp].rotdof, (double)species[jsp].rotdof, (double)species[ksp].rotdof};
    double xguess[] = {3000.0 , 2.0 , 2.0 , 2.0};
    //cout << " newtonTcol4 Function call (494)" << endl;
    double * newtinfo = newtonTcol4(4, nvibmode, postcoln.etotal, particle->species[isp].vibtemp, particle->species[jsp].vibtemp, particle->species[ksp].vibtemp, nrotmode, omega, xguess,
               1e-4,
               100);
    //cout << "newtonTcol4 successfully executed" << endl;
    Tcol = newtinfo[0];
    double zvibi = newtinfo[1];
    double zvibj = newtinfo[2];
    double zvibk = newtinfo[3];
   
    //cout << "Returned values to collide_vss: " << Tcol << " " << zvibi << " " << zvibj << " " << zvibk << endl;
    if ((Tcol < 0) || (zvibi < 0) || (zvibj < 0) || (zvibk < 0)) error->all(FLERR,"Negative returns from root solver");
    zc = zvibi + zvibj + zvibk + species[isp].rotdof + species[jsp].rotdof + species[ksp].rotdof + 10 - 2*aveomega12 - 2*aveomega123;

  } else if (numspecies == 2) {
    int isp = ip->ispecies;
    int jsp = jp->ispecies;


    //double aveomega =  CSomega[ip->ispecies][jp->ispecies]; 
    int nvibmode[] = {species[isp].nvibmode, species[jsp].nvibmode};

    //Need to calculate a collisional tempertature Tcol between the particles; in addition, need to find post-collisional vibrational degrees of 
    //freedom for colliding particles; values are then used for computation of the vibrational energies, assuming each vib mode is a harmonic oscillator
    //cout << "New Material Start (490)" << endl;

    double nrotmode[] = {(double)species[isp].rotdof, (double)species[jsp].rotdof};
    double xguess[] = {3000.0 , 2.0 , 2.0};
    //cout << " newtonTcol3 Function call (494)" << endl;
    double * newtinfo = newtonTcol3(3, nvibmode, postcoln.etotal, particle->species[isp].vibtemp, particle->species[jsp].vibtemp, nrotmode, aveomega, xguess,
               1e-4,
               100);
    //cout << "newtonTcol3 successfully executed" << endl;
    Tcol = newtinfo[0];
    double zvibi = newtinfo[1];
    double zvibj = newtinfo[2];
    //cout << "Returned values to collide_vss: " << Tcol << " " << zvibi << " " << zvibj << endl;
    if ((Tcol < 0) || (zvibi < 0) || (zvibj < 0)) error->all(FLERR,"Negative returns from root solver");
    zc = zvibi + zvibj + species[isp].rotdof + species[jsp].rotdof + 5 - 2*aveomega;
  }

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip; 
    else if (i == 1) p = jp; 
    else p = kp;

    int sp = p->ispecies;

    vibdof = species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        int **vibmode = particle->eiarray[particle->ewhich[index_vibmode]];
        int pindex = p - particle->particles;
        max_level = static_cast<int> 
          (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int> 
               (random->uniform()*(max_level+AdjustFactor));
          p->evib = ivib * update->boltz * species[sp].vibtemp[0];
          zki = (2 * species[sp].vibtemp[0] / Tcol) * (1 / (exp(species[sp].vibtemp[0]/Tcol) - 1));
          State_prob = pow(((E_Dispose - p->evib) / (E_Dispose)),
                                 (.5 * (zc - zki) - 1.0));
        } while (State_prob < random->uniform());
        E_Dispose -= p->evib;
        zc -= zki;
        //vibmode[pindex][0] = ivib;
      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(random->uniform(),(1.0 / (2.5-aveomega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
          p->evib = E_Dispose * 
          sample_bl(random,0.5*species[sp].vibdof-1.0,
                   1.5-aveomega);
          E_Dispose -= p->evib;
      } else if (vibdof > 2 && vibstyle == DISCRETE) {
          p->evib = 0.0;

          int nmode = particle->species[sp].nvibmode;
          int **vibmode = particle->eiarray[particle->ewhich[index_vibmode]];
          int pindex = p - particle->particles;
          double zpe = 0.0;
          for (int imode = 0; imode < nmode; imode++) {
            zpe += 0.5 * update->boltz * species[sp].vibtemp[imode];
          }
          for (int imode = 0; imode < nmode; imode++) {
            //ivib = vibmode[pindex][imode];
            //E_Dispose += ivib * update->boltz * 
            //particle->species[sp].vibtemp[imode];
            //double pevibwithzpe;
            max_level = static_cast<int>
            (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));
            do {
              ivib = static_cast<int> 
              (random->uniform()*(max_level+AdjustFactor));
              pevib = ivib * update->boltz * species[sp].vibtemp[imode];
              //pevibwithzpe = (ivib + .5) * update->boltz * species[sp].vibtemp[imode];
              zki = (2 * species[sp].vibtemp[0] / Tcol) * (1 / (exp(species[sp].vibtemp[0]/Tcol) - 1));
              State_prob = pow(((E_Dispose - pevib) / (E_Dispose)),
                                 (.5 * (zc - zki) - 1.0));
              //State_prob = pow(((E_Dispose - pevibwithzpe) / (E_Dispose - 0.5 * update->boltz * species[sp].vibtemp[imode])),
                                 //(.5 * (zc - zki) - 1.0));
              //cout << E_Dispose << " " << ivib << " " << State_prob << endl;
            } while (State_prob < random->uniform());

            vibmode[pindex][imode] = ivib;
            p->evib += pevib;
            E_Dispose -= pevib;
            //E_Dispose -= pevibwithzpe;
            zc -= zki;
            //cout << pevib << " " << zki << endl;
          }        
      }
    }
  } 


  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip; 
    else if (i == 1) p = jp; 
    else p = kp;

    int sp = p->ispecies;
    rotdof = species[sp].rotdof;

    if (rotdof) {
      if (rotstyle == NONE) {
        p->erot = 0.0;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(random->uniform(),(1/(2.5-aveomega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;
        
      } else if (rotdof > 2) {
        p->erot = E_Dispose * 
          sample_bl(random,0.5*species[sp].rotdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->erot;
      }
    }
  }
  
  // compute post-collision internal energies
  
  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;
  
  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }
  
  // compute portion of energy left over for scattering
  
  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
  //cout << "Pre React: " << precoln.etotal << " Pre_Energy: " << E_PreLB << " Vib E: " << postcoln.evib << " Rot E: " << postcoln.erot << " Trans E: " << postcoln.etrans << endl;
  if (E_Dispose < 0) error->all(FLERR,"Impossible Energy Post-collision");
}

/* ---------------------------------------------------------------------- */

double CollideVSS::sample_bl(RanKnuth *random, double Exp_1, double Exp_2)
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = random->uniform();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < random->uniform());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::rotrel(int isp, double Ec)
{
  // Because we are only relaxing one of the particles in each call, we only
  //  include its DoF, consistent with Bird 2013 (3.32)

  double Tr = Ec /(update->boltz * (2.5-params[isp][isp].omega + particle->species[isp].rotdof/2.0));
  double rotphi = (1.0+params[isp][isp].rotc2/sqrt(Tr) + params[isp][isp].rotc3/Tr)
                / params[isp][isp].rotc1;
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (3.5-params[isp][isp].omega));
  double vibphi = 1.0 / (params[isp][isp].vibc1/pow(Tr,params[isp][isp].omega) *
                         exp(params[isp][isp].vibc2/pow(Tr,1.0/3.0)));
  return vibphi;
}


/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideVSS::read_param_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open VSS parameter file %s",fname);
    error->one(FLERR,str);
  }

  // set all species diameters to -1, so can detect if not read
  // set all cross-species parameters to -1 to catch no-reads, as
  // well as user-selected average

  for (int i = 0; i < nparams; i++) {
    params[i][i].diam = -1.0;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].diam = params[i][j].omega = params[i][j].tref = -1.0;
      params[i][j].alpha = params[i][j].rotc1 = params[i][j].rotc2 = -1.0;
      params[i][j].rotc3 = params[i][j].vibc1 = params[i][j].vibc2 = -1.0;
    }
  }

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least REQWORDS, which depends on VARIABLE flag

  int REQWORDS = 5;
  if (relaxflag == VARIABLE) REQWORDS = 9;
  char **words = new char*[REQWORDS+1]; // one extra word in cross-species lines
  char line[MAXLINE];
  int isp,jsp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    int nwords = wordparse(REQWORDS+1,line,words);
    if (nwords < REQWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    jsp = particle->find_species(words[1]);
    
    // if we don't match a species with second word, but it's not a number,
    // skip the line (it involves a species we aren't using)

    //cout << isp << " " << jsp << endl;
    if ( jsp < 0 &&  !(atof(words[1]) > 0) ) continue;

    if (jsp < 0 ) {
      params[isp][isp].diam = atof(words[1]);
      params[isp][isp].omega = atof(words[2]);
      params[isp][isp].tref = atof(words[3]);
      params[isp][isp].alpha = atof(words[4]);
      if (relaxflag == VARIABLE) {
        params[isp][isp].rotc1 = atof(words[5]);
        params[isp][isp].rotc2 = atof(words[6]);
        params[isp][isp].rotc3 = (MY_PI+MY_PI2*MY_PI2)*params[isp][isp].rotc2;
        params[isp][isp].rotc2 = (MY_PI*MY_PIS/2.)*sqrt(params[isp][isp].rotc2);
        params[isp][isp].vibc1 = atof(words[7]);
        params[isp][isp].vibc2 = atof(words[8]);
      }
    }else {
      if (nwords < REQWORDS+1)  // one extra word in cross-species lines
        error->one(FLERR,"Incorrect line format in VSS parameter file");
      params[isp][jsp].diam = params[jsp][isp].diam = atof(words[2]);
      params[isp][jsp].omega = params[jsp][isp].omega = atof(words[3]);
      params[isp][jsp].tref = params[jsp][isp].tref = atof(words[4]);
      params[isp][jsp].alpha = params[jsp][isp].alpha = atof(words[5]);
      if (relaxflag == VARIABLE) {
        params[isp][jsp].rotc1 = params[jsp][isp].rotc1 = atof(words[6]);
        params[isp][jsp].rotc2 = atof(words[7]);
        params[isp][jsp].rotc3 = params[jsp][isp].rotc3 =
        		(MY_PI+MY_PI2*MY_PI2)*params[isp][jsp].rotc2;
        if(params[isp][jsp].rotc2 > 0)
        	params[isp][jsp].rotc2 = params[jsp][isp].rotc2 =
        			(MY_PI*MY_PIS/2.)*sqrt(params[isp][jsp].rotc2);
        params[isp][jsp].vibc1 = params[jsp][isp].vibc1= atof(words[8]);
        params[isp][jsp].vibc2 = params[jsp][isp].vibc2= atof(words[9]);
      }
    }
  }

  delete [] words;
  fclose(fp);

  // check that params were read for all species
  for (int i = 0; i < nparams; i++) {

    if (params[i][i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
	      particle->species[i].id);
      error->one(FLERR,str);
    }
  }

  for ( int i = 0; i<nparams; i++) {
    params[i][i].mr = particle->species[i].mass / 2;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].mr = params[j][i].mr = particle->species[i].mass *
	particle->species[j].mass / (particle->species[i].mass + particle->species[j].mass);

      if(params[i][j].diam < 0) params[i][j].diam = params[j][i].diam =
				  0.5*(params[i][i].diam + params[j][j].diam);
      if(params[i][j].omega < 0) params[i][j].omega = params[j][i].omega =
				   0.5*(params[i][i].omega + params[j][j].omega);
      if(params[i][j].tref < 0) params[i][j].tref = params[j][i].tref =
				  0.5*(params[i][i].tref + params[j][j].tref);
      if(params[i][j].alpha < 0) params[i][j].alpha = params[j][i].alpha =
				   0.5*(params[i][i].alpha + params[j][j].alpha);

      if (relaxflag == VARIABLE) {
	if(params[i][j].rotc1 < 0) params[i][j].rotc1 = params[j][i].rotc1 =
				     0.5*(params[i][i].rotc1 + params[j][j].rotc1);
	if(params[i][j].rotc2 < 0) params[i][j].rotc2 = params[j][i].rotc2 =
				     0.5*(params[i][i].rotc2 + params[j][j].rotc2);
	if(params[i][j].rotc3 < 0) params[i][j].rotc3 = params[j][i].rotc3 =
				     0.5*(params[i][i].rotc3 + params[j][j].rotc3);
	if(params[i][j].vibc1 < 0) params[i][j].vibc1 = params[j][i].vibc1 =
				     0.5*(params[i][i].vibc1 + params[j][j].vibc1);
	if(params[i][j].vibc2 < 0) params[i][j].vibc2 = params[j][i].vibc2 =
				     0.5*(params[i][i].vibc2 + params[j][j].vibc2);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
------------------------------------------------------------------------- */

int CollideVSS::wordparse(int maxwords, char *line, char **words)
{
  int nwords = 1;
  char * word;

  words[0] = strtok(line," \t\n");
  while ((word = strtok(NULL," \t\n")) != NULL && nwords < maxwords) {
    words[nwords++] = word;
  }
  return nwords;
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, int jsp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp][jsp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp][jsp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp][jsp].tref;
  else error->all(FLERR,"Request for unknown parameter from collide");
  return 0.0;
}

/* ged.cpp -- Gaussian elimination linear equation solvers.
 *
 *  (C) 2001, C. Bond. All rights reserved.
 *
 *  Simple pivoting on zero diagonal element supported.
 *   Eliminates unnecessary  zeroing of lower triangle.
 *   Does not scale rows to unity pivot value.
 *   Swaps b[] as well as a[][], so a pivot ID vector
 *   is not required.
 */

double * CollideVSS::gelimd3(double mat[3][4])
{
    static double res[3];
    int i,j,k;
    int n = 3;

/*    for(i=0;i<n;i++) 
    {                   
        for(j=i+1;j<n;j++)
        {
            if(abs(mat[i][i]) < abs(mat[j][i]))
            {
                for(k=0;k<n+1;k++)
                {
                    /* swapping mat[i][k] and mat[j][k] 
        mat[i][k]=mat[i][k]+mat[j][k];
                    mat[j][k]=mat[i][k]-mat[j][k];
                    mat[i][k]=mat[i][k]-mat[j][k];
                }
            }
      }
    }
 */  
     /* performing Gaussian elimination */
    for(i=0;i<n-1;i++)
    {
        for(j=i+1;j<n;j++)
        {
            float f=mat[j][i]/mat[i][i];
            for(k=0;k<n+1;k++)
            {
              mat[j][k]=mat[j][k]-f*mat[i][k];
      }
        }
    }
    /* Backward substitution for discovering values of unknowns */
    for(i=n-1;i>=0;i--)          
    {                     
        res[i]=mat[i][n];
                    
        for(j=i+1;j<n;j++)
        {
          if(i!=j)
          {
              res[i]=res[i]-mat[i][j]*res[j];
          }          
        }
        res[i]=res[i]/mat[i][i];
    }   
    //std::cout << res[0] << " " << res[1] << " " << res[2] << std::endl;
    return res;
}

double * CollideVSS::gelimd4(double mat[4][5])
{
    static double res[4];
    int i,j,k;
    int n = 4;

 /*   for(i=0;i<n;i++) 
    {                   
        for(j=i+1;j<n;j++)
        {
            if(abs(mat[i][i]) < abs(mat[j][i]))
            {
                for(k=0;k<n+1;k++)
                {
                    /* swapping mat[i][k] and mat[j][k] 
        mat[i][k]=mat[i][k]+mat[j][k];
                    mat[j][k]=mat[i][k]-mat[j][k];
                    mat[i][k]=mat[i][k]-mat[j][k];
                }
            }
      }
    } */
   
     /* performing Gaussian elimination */
    for(i=0;i<n-1;i++)
    {
        for(j=i+1;j<n;j++)
        {
            float f=mat[j][i]/mat[i][i];
            for(k=0;k<n+1;k++)
            {
              mat[j][k]=mat[j][k]-f*mat[i][k];
      }
        }
    }
    /* Backward substitution for discovering values of unknowns */
    for(i=n-1;i>=0;i--)          
    {                     
        res[i]=mat[i][n];
                    
        for(j=i+1;j<n;j++)
        {
          if(i!=j)
          {
              res[i]=res[i]-mat[i][j]*res[j];
          }          
        }
        res[i]=res[i]/mat[i][i];
    }   
    //std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << std::endl;
    return res;
}


/* ---------------------------------------------------------------------- */

double CollideVSS::nizenkov_zvib(int nmode, double Tcol, double zeta, double VibT[])
{
  double f;
  f = -zeta;
  if (nmode != 0) {
    for (int i = 0; i < nmode; i++) {
      f += (2 * VibT[i] / Tcol)*(1 / (exp(VibT[i] / Tcol) - 1));
     // std::cout << f << std::endl;
    }
  } 
  return f;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::nizenkov_dzvib(int nmode, double Tcol, double zeta, double VibT[])
{
  double f;
  f = 0.0;
  if (nmode != 0) {
    for (int i = 0; i < nmode; i++) {
      f += (2 * VibT[i] / pow(Tcol,3)) * (VibT[i]*exp(VibT[i] / Tcol) - Tcol*exp(VibT[i] / Tcol) + Tcol) * pow((1 / (exp(VibT[i] / Tcol) - 1)),2);
    }
  } 
  return f;
}

/* ----------------------------------------------------------------------
   compute post-reaction energy information
------------------------------------------------------------------------- */
double * CollideVSS::newtonTcol3(int n, int nmode[], double Ecol, double vibTempi[], double vibTempj[], double zrot[], double omega,
               double x0[],
               double tol,
               int nmax)
{
  // FUNCTION FOR CONVERTING COLLISIONAL ENERGY TO COLLISIONAL TEMPERATURE AND POST-COLLISION VIBRATIONAL DOFS
  // Computes Tcol and zeta_vib for particles using Newtons Search method and system of Equations from Nizenkov et al.
  // Search for values begins at some initial values "x0" until the search reaches a tolerance level "tol".

  //std::cout << "Now inside newtonTcol function" << std::endl;
  double f[3];
  double df1dx1, df1dx2, df1dx3, df1dx4, df2dx1, df3dx1, df4dx1;
  double * x, x_prev[3];
  double err[3];
  int i;

  double boltz = 1.38064852e-23;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels.

  // f and df are computed for Newton's search
  //std::cout << "Calling zvib and dzvib functions" << std::endl;
  //std::cout << Ecol << " " << x0[0] << " " << x0[1] << " " << x0[2] << " " << zrot[0] << " " << zrot[1] << " " << omega << std::endl;
  f[0] = -Ecol + .5 * boltz * (x0[1]+x0[2]+zrot[0]+zrot[1]+5-(2*omega)) * x0[0];
  //std::cout << "f0 calculated" << std::endl;
  f[1] = nizenkov_zvib(nmode[0],x0[0],x0[1],vibTempi);
  //std::cout << "f1 calculated" << std::endl;
  f[2] = nizenkov_zvib(nmode[1],x0[0],x0[2],vibTempj);
  //std::cout << f[0] << " " << f[1] << " " << f[2] << std::endl;

  df1dx1 = .5*(x0[1]+x0[2])*boltz;
  df1dx2 = df1dx3 = df1dx4 = .5*x0[0]*boltz;
  df2dx1 = nizenkov_dzvib(nmode[0],x0[0],x0[1],vibTempi);
  df3dx1 = nizenkov_dzvib(nmode[1],x0[0],x0[2],vibTempj);
  //std::cout << df1dx1 << " " << df1dx2 << " " << df2dx1 << " " << df3dx1 << std::endl;
  //std::cout << "zvib and dzvib functions completed" << std::endl;
  // Create Jacobian, then sends to solver
  double jac[3][4];
  jac[0][0] = df1dx1;
  jac[0][1] = df1dx2;
  jac[0][2] = df1dx3;

  jac[1][0] = df2dx1;
  jac[1][1] = -1.0;
  jac[1][2] = 0.0;

  jac[2][0] = df3dx1;
  jac[2][1] = 0.0;
  jac[2][2] = -1.0;

  jac[0][3] = f[0];
  jac[1][3] = f[1];
  jac[2][3] = f[2];
  //std::cout << "sending jac to gelimd3" << std::endl;
  x = gelimd3(jac);
  //std::cout << "gelimd3 completed" << std::endl;
  // Update guesses for variables and compute error
  for (int j = 0; j < n; j++) {
     //std::cout << *(x) << " " << *(x+1) << " " << *(x+2) << std::endl; 
     x[j] = x0[j] - x[j];
     err[j] = fabs(x[j]-x0[j]);
  }
  if (x[0] < 0.0) x[0] = 500;  //These checks occur to correct for any negative solutions returned, which occassionally occur when the difference between the root and guess are large.
  if (x[1] < 0.0) x[1] = 1.0;
  if (x[2] < 0.0) x[2] = 1.0;
 // std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  i=2;

  // Continue to search for Tvib until the error is less than the tolerance:
  while(((err[0] >= tol) || (err[1] >= tol) || (err[2] >= tol)) && (i <= nmax))
  {
    for (int j = 0; j < n; j++) {
       x_prev[j] = x[j];
    }

    f[0] = -Ecol + .5 * boltz * (x[1]+x[2]+zrot[0]+zrot[1]+5-(2*omega)) * x[0];
    f[1] = nizenkov_zvib(nmode[0],x[0],x[1],vibTempi);
    f[2] = nizenkov_zvib(nmode[1],x[0],x[2],vibTempj);
    //std::cout << f[0] << " " << f[1] << " " << f[2] << std::endl;

    df1dx1 = .5*(x[1]+x[2]+x[3])*boltz;
    df1dx2 = df1dx3 = df1dx4 = .5*x[0]*boltz;
    df2dx1 = nizenkov_dzvib(nmode[0],x0[0],x0[1],vibTempi);
    df3dx1 = nizenkov_dzvib(nmode[1],x0[0],x0[2],vibTempj);

    jac[0][0] = df1dx1;
    jac[0][1] = df1dx2;
    jac[0][2] = df1dx3;

    jac[1][0] = df2dx1;
    jac[1][1] = -1.0;
    jac[1][2] = 0.0;

    jac[2][0] = df3dx1;
    jac[2][1] = 0.0;
    jac[2][2] = -1.0;

    jac[0][3] = f[0];
    jac[1][3] = f[1];
    jac[2][3] = f[2];

    x = gelimd3(jac);

    for (int j = 0; j < n; j++) {
       x[j] = x_prev[j] - x[j];
       err[j] = fabs(x[j]-x_prev[j]);
    }
    if (x[0] < 0.0) x[0] = 500;  //These checks occur to correct for any negative solutions returned, which occassionally occur.
    if (x[1] < 0.0) x[1] = 1.0;
    if (x[2] < 0.0) x[2] = 1.0;
    i=i+1;

  }
  //std::cout << "returned vaules to collide vss: " << x[0] << " " << x[1] << " " << x[2] << std::endl;
  return x;

}

double * CollideVSS::newtonTcol4(int n, int nmode[], double Ecol, double vibTempi[], double vibTempj[], double vibTempk[], double zrot[], double omega[],
               double x0[],
               double tol,
               int nmax)
{
  // FUNCTION FOR CONVERTING COLLISIONAL ENERGY TO COLLISIONAL TEMPERATURE AND POST-COLLISION VIBRATIONAL DOFS
  // Computes Tcol and zeta_vib for particles using Newtons Search method and system of Equations from Nizenkov et al.
  // Search for values begins at some initial values "x0" until the search reaches a tolerance level "tol".


  double f[4];
  double df1dx1, df1dx2, df1dx3, df1dx4, df2dx1, df3dx1, df4dx1;
  double * x, x_prev[4];
  double err[4];
  int i;
  double boltz = 1.38064852e-23;
  //std::cout << "Now inside newtonTcol function" << std::endl;
  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels.

  // f and df are computed for Newton's search
  //std::cout << Ecol << " " << x0[0] << " " << x0[1] << " " << x0[2] << " " << x0[3] << " " << zrot[0] << " " << zrot[1] << " " << zrot[2] << " " << omega[0] << " " << omega[1] << " " << std::endl;
  f[0] = -Ecol + .5 * boltz * (x0[1]+x0[2]+x0[3]+zrot[0]+zrot[1]+zrot[2]+10 - 2*(omega[0]+omega[1])) * x0[0];
  //std::cout << f[0] << std::endl;
  f[1] = nizenkov_zvib(nmode[0],x0[0],x0[1],vibTempi);
  f[2] = nizenkov_zvib(nmode[1],x0[0],x0[2],vibTempj);
  f[3] = nizenkov_zvib(nmode[2],x0[0],x0[3],vibTempk);
  //std::cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;

  df1dx1 = .5*(x0[1]+x0[2]+x0[3])*boltz;
  df1dx2 = df1dx3 = df1dx4 = .5*x0[0]*boltz;
  df2dx1 = nizenkov_dzvib(nmode[0],x0[0],x0[1],vibTempi);
  df3dx1 = nizenkov_dzvib(nmode[1],x0[0],x0[2],vibTempj);
  df4dx1 = nizenkov_dzvib(nmode[2],x0[0],x0[3],vibTempk);
  //std::cout << "f and df functions called" << std::endl;

  // Create Jacobian, then sends to solver

  double jac[4][5];
  jac[0][0] = df1dx1;
  jac[0][1] = df1dx2;
  jac[0][2] = df1dx3;
  jac[0][3] = df1dx4;

  jac[1][0] = df2dx1;
  jac[1][1] = -1.0;
  jac[1][2] = 0.0;
  jac[1][3] = 0.0;

  jac[2][0] = df3dx1;
  jac[2][1] = 0.0;
  jac[2][2] = -1.0;  
  jac[2][3] = 0.0; 

  jac[3][0] = df4dx1;
  jac[3][1] = 0.0;
  jac[3][2] = 0.0;  
  jac[3][3] = -1.0; 

  jac[0][4] = f[0];
  jac[1][4] = f[1];
  jac[2][4] = f[2];
  jac[3][4] = f[3];
  //std::cout << "jac solver call" << std::endl;
  x = gelimd4(jac);
  //std::cout << "jac solver returned" << std::endl;
  //std::cout << *(x) << " " << *(x+1) << " " << *(x+2) << " " << *(x+3) << std::endl; 
  // Update guesses for variables and compute error
  for (int j = 0; j < n; j++) {
     //std::cout << *(x) << " " << *(x+1) << " " << *(x+2) << " " << *(x+3) << std::endl; 
     x[j] = x0[j] - x[j];
     err[j] = fabs(x[j]-x0[j]);
  }
  if (x[0] < 0.0) x[0] = 500;  //These checks occur to correct for any negative solutions returned, which occassionally occur.
  if (x[1] < 0.0) x[1] = 1.0;
  if (x[2] < 0.0) x[2] = 1.0;
  if (x[3] < 0.0) x[3] = 1.0;
  i=2;

  // Continue to search for Tvib until the error is less than the tolerance:
  while(((err[0] >= tol) || (err[1] >= tol) || (err[2] >= tol) || (err[3] >= tol)) && (i <= nmax))
  {
    for (int j = 0; j < n; j++) {
       x_prev[j] = x[j];
    }
    //std::cout << "loop f and df call" << std::endl;
    //std::cout << "loop x value: " << x[0] << " " << x[1] << " " << x[2] <<  " " << x[3] << std::endl;
    f[0] = -Ecol + .5 * boltz * (x[1]+x[2]+x[3]+zrot[0]+zrot[1]+zrot[2]+10 - 2*omega[0] - 2*omega[1]) * x[0];
    f[1] = nizenkov_zvib(nmode[0],x[0],x[1],vibTempi);
    f[2] = nizenkov_zvib(nmode[1],x[0],x[2],vibTempj);
    f[3] = nizenkov_zvib(nmode[2],x[0],x[3],vibTempk);

    df1dx1 = .5*(x[1]+x[2]+x[3])*boltz;
    df1dx2 = df1dx3 = df1dx4 = .5*x[0]*boltz;
    df2dx1 = nizenkov_dzvib(nmode[0],x0[0],x0[1],vibTempi);
    df3dx1 = nizenkov_dzvib(nmode[1],x0[0],x0[2],vibTempj);
    df4dx1 = nizenkov_dzvib(nmode[2],x0[0],x0[3],vibTempk);
    //std::cout << "loop f and df call complete" << std::endl;
    //std::cout << "loop f vals: " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;
    jac[0][0] = df1dx1;
    jac[0][1] = df1dx2;
    jac[0][2] = df1dx3;
    jac[0][3] = df1dx4;

    jac[1][0] = df2dx1;
    jac[1][1] = -1.0;
    jac[1][2] = 0.0;
    jac[1][3] = 0.0;

    jac[2][0] = df3dx1;
    jac[2][1] = 0.0;
    jac[2][2] = -1.0;  
    jac[2][3] = 0.0; 

    jac[3][0] = df4dx1;
    jac[3][1] = 0.0;
    jac[3][2] = 0.0;  
    jac[3][3] = -1.0; 

    jac[0][4] = f[0];
    jac[1][4] = f[1];
    jac[2][4] = f[2];
    jac[3][4] = f[3];

    x = gelimd4(jac);

    for (int j = 0; j < n; j++) {
       x[j] = x_prev[j] - x[j];
       err[j] = fabs(x[j]-x_prev[j]);
    }
    if (x[0] < 0.0) x[0] = 500;  //These checks occur to correct for any negative solutions returned, which occassionally occur.
    if (x[1] < 0.0) x[1] = 1.0;
    if (x[2] < 0.0) x[2] = 1.0;
    if (x[3] < 0.0) x[3] = 1.0;
    i=i+1;

  }
  //std::cout << "returned vaules to collide vss: " << x[0] << " " << x[1] << " " << x[2] <<  " " << x[3] << std::endl;
  return x;

}


