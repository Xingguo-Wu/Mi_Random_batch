/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Xingguo Wu (SJTU)
------------------------------------------------------------------------- */

#include "mi_rbe.h"
#include <mpi.h>
#include "atom.h"
#include "comm.h"
#include "force.h" 
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "utils.h"
#include <cstring>
#include <iostream>
#include <ctime>
 
#include <unistd.h>
#include <stdexcept>
#include <cassert>
#include <cstddef> 
#include <vector>
#include <array>
#include <algorithm>
#include <random>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;
#define SMALL 0.00001
#define LARGE 1000000
inline double ProbMiRbeIso(double m, double alpha, double L, double S1, double pi)
{
	double k = 2 * pi * m / L;
	return exp(-(k) * (k) / (4 * alpha)) / S1;
}

static double gaussian_mirbe_iso(double k2, double alpha)
{
    return exp(-k2 / (4 * alpha));
}

static double mi_MH_D_iso(double m, double alpha, double L, double pi)
{
	double a;
	if (m == 0 || fabs(m) < 1.0e-13)
	{
		a = erf(0.5 / sqrt(alpha * L * L / (pi * pi)));
	}
	else
	{
		a = 0.5 * (erf((abs(m) + 0.5) / sqrt(alpha * L * L / (pi * pi))) - erf((abs(m) - 0.5) / sqrt(alpha * L * L / (pi * pi))));
	}
	return a;
}

/* ---------------------------------------------------------------------- */
MiRbe::MiRbe(LAMMPS* lmp) : KSpace(lmp)
{
	TimeSet = new double[1000];
	K_All = nullptr;
	K_Sample = nullptr;
	P = 0;
	RankID = 0;
	mirbe_rho = 0.0;
	mirbe_alpha = 0.0;
	mirbe_eps = 1.0e-8;
	random_seed_base = 0U;
	mirbe_updates = 0;
	mirbe_state_index = -1;
	write_timing = 0;
	force_s0 = 1.0;
	mirbe_gain_index = -1;
	ewaldflag = 1;
}
 
void MiRbe::settings(int narg, char** arg)
{
	MPI_Comm_rank(world, &me);
	int nprocs;
	MPI_Comm_size(world, &nprocs);

	if ((narg != 7) && (narg != 8) && (narg != 9))
		error->all(FLERR, "Illegal kspace_style mi/rbe command: expected accuracy P rho alpha eps S0 seed [RankID] [write_timing]");
	accuracy_relative = fabs(utils::numeric(FLERR,arg[0],false,lmp));
	if (accuracy_relative > 1.0)
		error->all(FLERR, "Invalid relative accuracy {:g} for kspace_style {}",
				accuracy_relative, force->kspace_style);
	P = int(utils::numeric(FLERR, arg[1], false, lmp));
	mirbe_rho = utils::numeric(FLERR, arg[2], false, lmp);
	mirbe_alpha = utils::numeric(FLERR, arg[3], false, lmp);
	mirbe_eps = utils::numeric(FLERR, arg[4], false, lmp);
	force_s0 = utils::numeric(FLERR, arg[5], false, lmp);
	if (force_s0 <= 0.0) error->all(FLERR, "Mi-RBE S0 must be > 0");
	random_seed_base = static_cast<unsigned int>(utils::inumeric(FLERR, arg[6], false, lmp));
	RankID = (narg >= 8) ? int(utils::numeric(FLERR, arg[7], false, lmp)) : nprocs;
	if (RankID <= 0 || RankID > nprocs)
		error->all(FLERR, "Kspace style mi/rbe requires 0 < RankID <= MPI size");
	write_timing = (narg == 9) ? int(utils::numeric(FLERR, arg[8], false, lmp)) : 0;
	mirbe_updates = 0;

	pi = 3.141592653589793;
	pi_sqrt = 1.772453850905516;

	K_All = new int[100][3];
	K_Sample = new int* [P];
	for (int i = 0; i < P; i++) K_Sample[i] = new int[3];
}

/* ---------------------------------------------------------------------- */

void MiRbe::init()
{
	if (comm->me == 0) utils::logmesg(lmp,"Ewald initialization ...\n");

	// error check

	triclinic_check();
	if (domain->dimension == 2)
		error->all(FLERR,"Cannot use Ewald with 2d simulation");

	if (!atom->q_flag) error->all(FLERR,"Kspace style requires atom attribute q");

	if (mirbe_alpha != 0.0) {
		int flag, cols;
		mirbe_state_index = atom->find_custom("mirbe_state", flag, cols);
		if (mirbe_state_index < 0 || flag != 1 || cols != 0)
			error->all(FLERR, "Kspace style mi/rbe requires 'fix <id> all property/atom d_mirbe_state'");
	}
	int flag, cols;
	mirbe_gain_index = atom->find_custom("mirbe_gain", flag, cols);

	if (slabflag == 0 && domain->nonperiodic > 0)
		error->all(FLERR,"Cannot use non-periodic boundaries with Ewald");
	if (slabflag) {
		if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
			domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
		error->all(FLERR,"Incorrect boundaries with slab Ewald");
		if (domain->triclinic)
		error->all(FLERR,"Cannot (yet) use Ewald with triclinic box "
					"and slab correction");
	}

	// compute two charge force

	two_charge();

	// extract short-range Coulombic cutoff from pair style

	triclinic = domain->triclinic;
	pair_check();

	int itmp;
	auto p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
	if (p_cutoff == nullptr)
		error->all(FLERR,"KSpace style is incompatible with Pair style");
	double cutoff = *p_cutoff;

	// compute qsum & qsqsum and warn if not charge-neutral

	scale = 1.0;
	qqrd2e = force->qqrd2e;
	qsum_qsq();
	natoms_original = atom->natoms;

	// set accuracy (force units) from accuracy_relative or accuracy_absolute

	if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
	else accuracy = accuracy_relative * two_charge_force;

	// setup K-space resolution

	bigint natoms = atom->natoms;

	// use xprd,yprd,zprd even if triclinic so grid size is the same
	// adjust z dimension for 2d slab Ewald
	// 3d Ewald just uses zprd since slab_volfactor = 1.0

	double xprd = domain->xprd;
	double yprd = domain->yprd;
	double zprd = domain->zprd;
	double zprd_slab = zprd*slab_volfactor;

	// make initial g_ewald estimate
	// based on desired accuracy and real space cutoff
	// fluid-occupied volume used to estimate real-space error
	// zprd used rather than zprd_slab

	if (!gewaldflag) {
		if (accuracy <= 0.0)
		error->all(FLERR,"KSpace accuracy must be > 0");
		if (q2 == 0.0)
		error->all(FLERR,"Must use 'kspace_modify gewald' for uncharged system");
		g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
		if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
		else g_ewald = sqrt(-log(g_ewald)) / cutoff;
	}

	// setup Ewald coefficients so can print stats

	setup();

	if (me == 0) {
    	std::string mesg = fmt::format("  G vector (1/distance) = {:.8g}\n",g_ewald);
    	utils::logmesg(lmp,mesg);
  	}
}


/* ---------------------------------------------------------------------- */

void MiRbe::setup()
{

	int me;
	MPI_Comm_rank(world, &me);
	alpha = g_ewald * g_ewald;
	double xprd = domain->xprd;
	double yprd = domain->yprd;
	double zprd = domain->zprd * slab_volfactor;
	volume = xprd * yprd * zprd;

	S_X = sqrt(alpha * xprd * xprd / MY_PI) * (1 + 2 * exp(-(alpha * xprd * xprd)));
	S_Y = sqrt(alpha * yprd * yprd / MY_PI) * (1 + 2 * exp(-(alpha * yprd * yprd)));
	S_Z = sqrt(alpha * zprd * zprd / MY_PI) * (1 + 2 * exp(-(alpha * zprd * zprd)));
	S = S_X * S_Y * S_Z - 1;
    Step = 0;
}

/* ---------------------------------------------------------------------- */

void MiRbe::compute(int eflag, int vflag)
{

	ev_init(eflag, vflag);
	alpha = g_ewald * g_ewald;
	// if(me==0){
	// 	cout << "g_ewald: " << g_ewald << ' ' << "alpha: " << alpha << endl;
	// 	cout << "S: " << S << endl;
	// }

	// if atom count has changed, update qsum and qsqsum
	 
	if (atom->natoms != natoms_original) {
		qsum_qsq();
		natoms_original = atom->natoms;
	}

	double xprd = domain->xprd;
	double yprd = domain->yprd;
	double zprd = domain->zprd * slab_volfactor;
	scale = 1;

	double V = xprd * yprd * zprd;
	std::vector<std::array<double,3>> K(P);
	Step = update->ntimestep;

	const int this_rank = Step % RankID;

	if (((Step % RankID) == 0) && (me < RankID)) {
		S1_X = sqrt(alpha * xprd * xprd / MY_PI) * (1 + 2 * exp(-(alpha * xprd * xprd)));
		S1_Y = sqrt(alpha * yprd * yprd / MY_PI) * (1 + 2 * exp(-(alpha * yprd * yprd)));
		S1_Z = sqrt(alpha * zprd * zprd / MY_PI) * (1 + 2 * exp(-(alpha * zprd * zprd)));
		S0 = S1_X * S1_Y * S1_Z - 1;
		xprd0 = xprd;
		yprd0 = yprd;
		zprd0 = zprd;

		std::mt19937 rng(random_seed_base + 104729U * static_cast<unsigned int>(me) +
		                 1000003U * static_cast<unsigned int>(Step));
		std::uniform_real_distribution<double> uniform01(0.0, 1.0);
		auto normal_x = [&]() {
			std::normal_distribution<double> dist(0.0, sqrt(alpha * xprd * xprd / (2 * pi * pi)));
			return dist(rng);
		};
		auto normal_y = [&]() {
			std::normal_distribution<double> dist(0.0, sqrt(alpha * yprd * yprd / (2 * pi * pi)));
			return dist(rng);
		};
		auto normal_z = [&]() {
			std::normal_distribution<double> dist(0.0, sqrt(alpha * zprd * zprd / (2 * pi * pi)));
			return dist(rng);
		};

		std::vector<int> mx(5 * P), my(5 * P), mz(5 * P);
		do {
			mx[0] = static_cast<int>(round(normal_x()));
			my[0] = static_cast<int>(round(normal_y()));
			mz[0] = static_cast<int>(round(normal_z()));
		} while (mx[0] == 0 && my[0] == 0 && mz[0] == 0);

		for (int i = 0; i < 5 * P - 1; i++) {
			double x = normal_x();
			double mold = mx[i];
			double mnew = round(x);
			double pup = ProbMiRbeIso(mnew, alpha, xprd, S1_X, pi);
			double qup = mi_MH_D_iso(mold, alpha, xprd, pi);
			double pdown = ProbMiRbeIso(mold, alpha, xprd, S1_X, pi);
			double qdown = mi_MH_D_iso(mnew, alpha, xprd, pi);
			double acce = pup * qup / (pdown * qdown) > 1.0 ? 1.0 : pup * qup / (pdown * qdown);
			const bool accept_trial_x = uniform01(rng) < acce;
			mx[i + 1] = accept_trial_x ? static_cast<int>(mnew) : static_cast<int>(mold);

			double xx = normal_y();
			mold = my[i];
			mnew = round(xx);
			pup = ProbMiRbeIso(mnew, alpha, yprd, S1_Y, pi);
			qup = mi_MH_D_iso(mold, alpha, yprd, pi);
			pdown = ProbMiRbeIso(mold, alpha, yprd, S1_Y, pi);
			qdown = mi_MH_D_iso(mnew, alpha, yprd, pi);
			acce = pup * qup / (pdown * qdown) > 1.0 ? 1.0 : pup * qup / (pdown * qdown);
			const bool accept_trial_y = uniform01(rng) < acce;
			my[i + 1] = accept_trial_y ? static_cast<int>(mnew) : static_cast<int>(mold);

			double xxx = normal_z();
			mold = mz[i];
			mnew = round(xxx);
			pup = ProbMiRbeIso(mnew, alpha, zprd, S1_Z, pi);
			qup = mi_MH_D_iso(mold, alpha, zprd, pi);
			pdown = ProbMiRbeIso(mold, alpha, zprd, S1_Z, pi);
			qdown = mi_MH_D_iso(mnew, alpha, zprd, pi);
			acce = pup * qup / (pdown * qdown) > 1.0 ? 1.0 : pup * qup / (pdown * qdown);
			const bool accept_trial_z = uniform01(rng) < acce;
			mz[i + 1] = accept_trial_z ? static_cast<int>(mnew) : static_cast<int>(mold);

			if (mx[i + 1] == 0 && my[i + 1] == 0 && mz[i + 1] == 0) {
				i = i - 1;
			}
		}
		for (int i = 0; i < P; i++) {
			K_Sample[i][0] = mx[5*i + 4];
			K_Sample[i][1] = my[5*i + 4];
			K_Sample[i][2] = mz[5*i + 4];
		}
	}

	if (comm->me == this_rank) {
		for (int i = 0; i < P; i++) {
			K[i][0] = static_cast<double>(K_Sample[i][0]);
			K[i][1] = static_cast<double>(K_Sample[i][1]);
			K[i][2] = static_cast<double>(K_Sample[i][2]);
		}
	}

	MPI_Bcast(reinterpret_cast<double *>(K.data()), 3 * P, MPI_DOUBLE, this_rank, MPI_COMM_WORLD);

	double** x = atom->x;
	double** f = atom->f;
	double* mirbe_state = (mirbe_alpha != 0.0) ? atom->dvector[mirbe_state_index] : nullptr;
	double* mirbe_gain = (mirbe_gain_index >= 0) ? atom->dvector[mirbe_gain_index] : nullptr;
	double* q = atom->q;
	int* type = atom->type;
	int nlocal = atom->nlocal;
	double qqrd2e = force->qqrd2e;
	double dielectric = force->dielectric;

	std::vector<std::array<double,3>> X(nlocal), F(nlocal);
	std::vector<double> Q(nlocal);
	for (int i = 0; i < nlocal; i++) {
		X[i][0] = x[i][0];
		X[i][1] = x[i][1];
		X[i][2] = x[i][2];
		F[i] = {0.0, 0.0, 0.0};
		Q[i] = q[i];
	}

	const std::array<double,3> pxyz = {2 * pi / xprd, 2 * pi / yprd, 2 * pi / zprd};
	std::vector<std::array<double,2>> Rho(P), Rho_All(P);
	std::vector<std::array<double,3>> midterm(P), Kx(P);
	std::vector<double> fac(P);

	for (int i = 0; i < P; i++) {
		Kx[i][0] = K[i][0] * pxyz[0];
		Kx[i][1] = K[i][1] * pxyz[1];
		Kx[i][2] = K[i][2] * pxyz[2];
		const double Kx0 = Kx[i][0] * xprd / xprd0;
		const double Ky0 = Kx[i][1] * yprd / yprd0;
		const double Kz0 = Kx[i][2] * zprd / zprd0;
		const double K2 = Kx[i][0] * Kx[i][0] + Kx[i][1] * Kx[i][1] + Kx[i][2] * Kx[i][2];
		const double K2_0 = Kx0 * Kx0 + Ky0 * Ky0 + Kz0 * Kz0;
		fac[i] = S0 / S * gaussian_mirbe_iso(K2, alpha) / gaussian_mirbe_iso(K2_0, alpha);
		midterm[i][0] = Kx[i][0] / K2;
		midterm[i][1] = Kx[i][1] / K2;
		midterm[i][2] = Kx[i][2] / K2;
	}

	for (int i = 0; i < P; i++) {
		double real = 0.0;
		double imag = 0.0;
		for (int j = 0; j < nlocal; j++) {
			const double moment = Kx[i][0] * X[j][0] + Kx[i][1] * X[j][1] + Kx[i][2] * X[j][2];
			real += Q[j] * cos(moment);
			imag += Q[j] * sin(moment);
		}
		Rho[i][0] = real;
		Rho[i][1] = imag;
	}

	MPI_Allreduce(reinterpret_cast<double *>(Rho.data()), reinterpret_cast<double *>(Rho_All.data()),
	              2 * P, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	const double MIDTERM = -4 * pi * (S0 / (P + 0.00)) * qqrd2e / V;
	for (int i = 0; i < nlocal; i++) {
		for (int j = 0; j < P; j++) {
			const double moment = -(Kx[j][0] * X[i][0] + Kx[j][1] * X[i][1] + Kx[j][2] * X[i][2]);
			const double c = cos(moment);
			const double s = sin(moment);
			const double imag = (c * Rho_All[j][1] + s * Rho_All[j][0]) * MIDTERM * fac[j];
			F[i][0] += Q[i] * imag * midterm[j][0];
			F[i][1] += Q[i] * imag * midterm[j][1];
			F[i][2] += Q[i] * imag * midterm[j][2];
		}
	}

	// mirbe_alpha == 0 is the unmodified RBE baseline.  Keep it as a
	// separate path so rho and the second-moment state cannot affect the
	// force, even through otherwise unnecessary floating-point operations.
	if (mirbe_alpha == 0.0) {
		for (int i = 0; i < nlocal; i++) {
			if (mirbe_gain) mirbe_gain[i] = 1.0;
			for (int dim = 0; dim < 3; dim++) f[i][dim] += F[i][dim];
		}
	} else {
		const double ema_weight = 1.0 - mirbe_rho;
		const double whitening_power = 0.5 * mirbe_alpha;
		const bool has_mirbe_history = mirbe_updates > 0;
		double safe_bias_correction = 1.0;
		if (has_mirbe_history) {
			const double bias_correction = 1.0 - pow(mirbe_rho, static_cast<double>(mirbe_updates));
			safe_bias_correction = (bias_correction > 1.0e-16) ? bias_correction : 1.0e-16;
		}

		for (int i = 0; i < nlocal; i++) {
			const double force2_over3 =
				(F[i][0] * F[i][0] + F[i][1] * F[i][1] + F[i][2] * F[i][2]) / 3.0;

			double gain = 1.0;
			if (has_mirbe_history) {
				const double corrected_state = mirbe_state[i] / safe_bias_correction;
				gain = pow((corrected_state + mirbe_eps) / force_s0, -whitening_power);
			}
			mirbe_state[i] = mirbe_rho * mirbe_state[i] + ema_weight * force2_over3;

			if (mirbe_gain) mirbe_gain[i] = gain;
			for (int dim = 0; dim < 3; dim++) f[i][dim] += gain * F[i][dim];
		}
		mirbe_updates++;
	}

	// sum global energy across Kspace vevs and add in volume-dependent term
		const double qscale = qqrd2e * scale;

        if (eflag_global) {
			double KXX[3];
			for (int i = 0; i < P; i++)
			{
				KXX[0] = K[i][0] * pxyz[0];
				KXX[1] = K[i][1] * pxyz[1];
				KXX[2] = K[i][2] * pxyz[2];

				energy += fac[i] * (2 * pi * S0 / ((P + 0.00) * V)) * (Rho_All[i][0] * Rho_All[i][0] + Rho_All[i][1] * Rho_All[i][1]) / (KXX[0] * KXX[0] + KXX[1] * KXX[1] + KXX[2] * KXX[2]);
			}
			energy -= g_ewald * qsqsum / MY_PIS +
				MY_PI2 * qsum * qsum / (g_ewald * g_ewald * volume);// self energy
			energy *= qscale;
			// if(me == 0){cout << "energy: " << energy << endl;}
		}


// global virial
		if (vflag_global) {
			double KKXX[3],Moment_Term= 2 * pi * (S0 / (P + 0.00)) * qqrd2e / (V+0.00);

#if defined(LMP_SIMD_COMPILER)
#pragma vector aligned
#pragma simd
#endif
			for (int i = 0; i < P; i++)
			{
				KKXX[0] = K[i][0] * pxyz[0];
				KKXX[1] = K[i][1] * pxyz[1]; 
				KKXX[2] = K[i][2] * pxyz[2];

				double coef1 = fac[i] * (Rho_All[i][0]* Rho_All[i][0]+ Rho_All[i][1]* Rho_All[i][1])/ (KKXX[0] * KKXX[0] + KKXX[1] * KKXX[1] + KKXX[2] * KKXX[2]);
				double coef2 = (1.0 / (4 * g_ewald * g_ewald)) + 1.0 / (KKXX[0] * KKXX[0] + KKXX[1] * KKXX[1] + KKXX[2] * KKXX[2]);
				virial[0] += coef1 * Moment_Term * (1 - coef2 * 2 * KKXX[0] * KKXX[0]);
				virial[1] += coef1 * Moment_Term * (1 - coef2 * 2 * KKXX[1] * KKXX[1]);
				virial[2] += coef1 * Moment_Term * (1 - coef2 * 2 * KKXX[2] * KKXX[2]);

				virial[3] += coef1 * Moment_Term * (- coef2 * 2 * KKXX[0] * KKXX[1]);
				virial[4] += coef1 * Moment_Term * (- coef2 * 2 * KKXX[0] * KKXX[2]);
				virial[5] += coef1 * Moment_Term * (- coef2 * 2 * KKXX[1] * KKXX[2]);
			}
			/*if(comm->me == 0){
			cout << "virial: " << virial[0] << " " << virial[1] << " " << virial[2] << " " << virial[3] << " " 
			<< virial[4] << " " << virial[5] << endl;
			}*/
		}

// per-atom energy/virial
// energy includes self-energy correction
		if (evflag_atom) {
			if (eflag_atom) {
#if defined(LMP_SIMD_COMPILER)
#pragma vector aligned
#pragma simd
#endif
				for (int i = 0; i < nlocal; i++) {
					eatom[i] -= g_ewald * q[i] * q[i] / sqrt(pi) + pi*pi * q[i] * qsum /
						(g_ewald * g_ewald * volume);
					eatom[i] *= qscale;
				}
			}

		if (vflag_atom)
#if defined(LMP_SIMD_COMPILER)
#pragma vector aligned
#pragma simd
#endif
				for (int i = 0; i < nlocal; i++)
					for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
		}
    
	if (slabflag == 1) slabcorr();
}


/* ---------------------------------------------------------------------- */
MiRbe::~MiRbe()
{
	if (K_Sample) {
		for (int i = 0; i < P; i++) delete[] K_Sample[i];
		delete[] K_Sample;
	}
	delete[] K_All;
	delete[] TimeSet;
}

/* ---------------------------------------------------------------------- */

double MiRbe::memory_usage()
{
	double bytes = 0.0;
	bytes += sizeof(double) * 1000;
	bytes += sizeof(int) * 100 * 3;
	bytes += sizeof(int) * P * 3;
	return bytes;
}

/* ---------------------------------------------------------------------- */

void MiRbe::slabcorr()
{
	// compute local contribution to global dipole moment
	qqrd2e = force->qqrd2e;
	scale = 1;

	double* q = atom->q;
	double** x = atom->x;
	double zprd = domain->zprd;
	int nlocal = atom->nlocal;

	double dipole = 0.0;
	for (int i = 0; i < nlocal; i++) dipole += q[i] * x[i][2];

	// sum local contributions to get global dipole moment

	double dipole_all;
	MPI_Allreduce(&dipole, &dipole_all, 1, MPI_DOUBLE, MPI_SUM, world);

	// need to make non-neutral systems and/or
	//  per-atom energy translationally invariant

	double dipole_r2 = 0.0;
	if (eflag_atom || fabs(qsum) > SMALL) {
		for (int i = 0; i < nlocal; i++)
			dipole_r2 += q[i] * x[i][2] * x[i][2];

		// sum local contributions

		double tmp;
		MPI_Allreduce(&dipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
		dipole_r2 = tmp;
	}

	// compute corrections

	const double e_slabcorr = MY_2PI * (dipole_all * dipole_all -
		qsum * dipole_r2 - qsum * qsum * zprd * zprd / 12.0) / volume;
	const double qscale = qqrd2e * scale;

	if (eflag_global) energy += qscale * e_slabcorr;


	// per-atom energy

	// add on force corrections

	double ffact = qscale * (-4.0 * MY_PI / volume);
	double** f = atom->f;

	for (int i = 0; i < nlocal; i++) f[i][2] += ffact * q[i] * (dipole_all - qsum * x[i][2]);

}

