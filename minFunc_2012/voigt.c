/* voigt.c: computes exact Voigt profiles in terms of the complex
   error function (requires libcerf) */

#include "mex.h"
#include <cerf.h>
#include <math.h>

#define LAMBDAS_ARG    prhs[0]
#define Z_ARG          prhs[1]
#define N_ARG          prhs[2]
#define NUM_LINES_ARG  prhs[3]

#define PROFILE_ARG    plhs[0]

/* number of lines in the Lyman series to consider */
#define NUM_LINES 31

/* note: all units are CGS */

/* physical constants */

static const double c   = 2.99792458e+10;              /* speed of light          cm s⁻¹        */
/* static const double k   = 1.38064852e-16; */        /* Boltzmann constant      erg K⁻¹       */
/* static const double m_p = 1.672621898e-24;*/        /* proton mass             g             */
/* static const double m_e = 9.10938356e-28; */        /* electron mass           g             */
/* e = 1.6021766208e-19 * c / 10; */
/* static const double e   = 4.803204672997660e-10; */ /* elementary charge       statC         */

/* CIV doublet */

static const double transition_wavelengths[] =         /* transition wavelengths  cm            */
  {
    1.54820490e-05,
    1.55077845e-05
  };

static const double oscillator_strengths[] =           /* oscillator strengths    dimensionless */
  {
    0.189900,
    0.094750
  };

static const double Gammas[] =                         /* transition rates        s^-1          */
  {
    2.642e+08,
    2.628e+08
  };

/* assumed constant */
/* static const double T = 1e+04; */                   /* gas temperature         K             */

/* derived constants */

/* b = sqrt(2 * k * T / m_p); */
/* static const double b =
     1.28486551932562422e+06; */                       /* Doppler parameter       cm s⁻¹        */

/* sigma = b / M_SQRT2; */
static const double sigma = 9.08537121627923800e+05;   /* Gaussian width          cm s⁻¹        */

/* leading_constants[i] =
       M_PI * e * e * oscillator_strengths[i] * transition_wavelengths[i] / (m_e * c);
*/
static const double leading_constants[] =              /* leading constants       cm²           */
  {
    1.34347262962625339e-07,
    2.15386482180851912e-08,
    7.48525170087141461e-09,
    3.51375347286007472e-09,
    1.94112336271172934e-09,
    1.18916112899713152e-09,
    7.82448627128742997e-10,
    5.42930932279390593e-10,
    3.92301197282493829e-10,
    2.92796010451409027e-10,
    2.24422239410389782e-10,
    1.75895684469038289e-10,
    1.40338556137474778e-10,
    1.13995374637743197e-10,
    9.37706429662300083e-11,
    7.79453203101192392e-11,
    6.55369055970184901e-11,
    5.58100321584169051e-11,
    4.77895916635794548e-11,
    4.12301389852588843e-11,
    3.58872072638707592e-11,
    3.12745536798214080e-11,
    2.76337116167110415e-11,
    2.44791750078032772e-11,
    2.15681362798480253e-11,
    1.93850080479346101e-11,
    1.72025364178111889e-11,
    1.55051698336865945e-11,
    1.40504672409331934e-11,
    1.28383057589411395e-11,
    1.16264059622218997e-11
  };

/* gammas[i] = Gammas[i] * transition_wavelengths[i] / (4 * M_PI); */
static const double gammas[] =                         /* Lorentzian widths       cm s⁻¹        */
  {
    6.06075804241938613e+02,
    1.54841462408931704e+02,
    6.28964942715328164e+01,
    3.17730561586147395e+01,
    1.82838676775503330e+01,
    9.15463131005758157e+00,
    6.08448802613156925e+00,
    4.24977523573725779e+00,
    3.08542121666345803e+00,
    2.31184525202557767e+00,
    1.77687796208123139e+00,
    1.39477990932179852e+00,
    1.11505539984541979e+00,
    9.05885451682623022e-01,
    7.45877170715450677e-01,
    6.21261624902197052e-01,
    5.22994533400935269e-01,
    4.44469874827484512e-01,
    3.80923210837841919e-01,
    3.28912390446060132e-01,
    2.85949711597237033e-01,
    2.50280032040928802e-01,
    2.20224061101442048e-01,
    1.94686521675913549e-01,
    1.73082093051965591e-01,
    1.54536566013816490e-01,
    1.38539175663870029e-01,
    1.24652675945279762e-01,
    1.12585442799479921e-01,
    1.02045988802423507e-01,
    9.27433783998286437e-02
  };

/* BOSS spectrograph instrumental broadening */
/* R = 2000; */                                  /* resolving power          dimensionless */

/* width of instrument broadening in pixels */
/* pixel_spacing = 1e-4; */
/* pixel_sigma = 1 / (R * 2 * sqrt(2 * M_LN2) * (pow(10, pixel_spacing) - 1)); */

static const int width = 3;                      /* width of convolution     dimensionless */

/*
  total = 0;
  for (i = -width, j = 0; i <= width; i++, j++) {
    instrument_profile[j] = exp(-0.5 * i * i / (pixel_sigma * pixel_sigma));
    total += instrument_profile[j];
  }

  for (i = 0; i < 2 * width + 1; i++)
    instrument_profile[i] /= total;
*/

static const double instrument_profile[] =
  {
    2.17460992138080811e-03,
    4.11623059580451742e-02,
    2.40309364651846963e-01,
    4.32707438937454059e-01, /* center pixel */
    2.40309364651846963e-01,
    4.11623059580451742e-02,
    2.17460992138080811e-03
  };

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {

  double *lambdas, *profile, *multipliers, *raw_profile;
  double z, N, velocity, total;
  int num_lines, i, j, k;
  mwSize num_points;

  /* get input */
  lambdas = mxGetPr(LAMBDAS_ARG);                /* wavelengths             Å             */
  z       = mxGetScalar(Z_ARG);                  /* redshift                dimensionless */
  N       = mxGetScalar(N_ARG);                  /* column density          cm⁻²          */

  num_lines = (nrhs > 3) ? (int)(mxGetScalar(NUM_LINES_ARG)) : NUM_LINES;

  num_points = mxGetNumberOfElements(LAMBDAS_ARG);

  /* initialize output */
  PROFILE_ARG = mxCreateDoubleMatrix(num_points - 2 * width, 1, mxREAL);
  profile = mxGetPr(PROFILE_ARG);                /* absorption profile      dimensionless */

  /* to hold the profile before instrumental broadening */
  raw_profile = mxMalloc(num_points * sizeof(double));

  multipliers = mxMalloc(num_lines * sizeof(double));
  for (i = 0; i < num_lines; i++)
    multipliers[i] = c / (transition_wavelengths[i] * (1 + z)) / 1e8;

  /* compute raw Voigt profile */
  for (i = 0; i < num_points; i++) {
    /* apply each absorption line */
    total = 0;
    for (j = 0; j < num_lines; j++) {
      /* velocity relative to transition wavelength */
      velocity = lambdas[i] * multipliers[j] - c;
      total += -leading_constants[j] * voigt(velocity, sigma, gammas[j]);
    }

    raw_profile[i] = exp(N * total);
  }

  num_points = mxGetNumberOfElements(PROFILE_ARG);

  /* instrumental broadening */
  for (i = 0; i < num_points; i++)
    for (j = i, k = 0; j <= i + 2 * width; j++, k++)
      profile[i] += raw_profile[j] * instrument_profile[k];

  mxFree(raw_profile);
  mxFree(multipliers);

}
