catalog.mat
-----------

This file contains basic information related to each line of sight
extracted from previously published catalogs. There are a total of N =
297301 QSOs considered, representing the entire DR12Q catalog. The
variables are:

* `bal_visual_flags`: binary flags indicating whether the QSO was
  identified as having broad absorption line (BAL) features during
  visual inspection. 1 indicates BAL QSO, 0 indicates non-BAL
  QSO. (DR12Q catalog, variable `BAL_FLAG_VI`)
* `decs`: declination in degrees (J2000) (DR12Q catalog, variable
  `DEC`)
* `fiber_ids`: SDSS spectroscopic fiber identifier (DR12Q catalog,
  variable `FIBERID`)
* `filter_flags`: bitfields indicating whether a given line of sight
  was rejected during the DLA search and why. A value of 0 means the
  line of sight was searched for DLAs. For a nonzero value, the
  individual bits of the flags may be inspected:

  - bit 0 (least significant): 1 indicates z_QSO < 2.15
  - bit 1: 1 indicates BAL QSO (see `bal_visual_flags`)
  - bit 2: 1 indicates spectrum could not be normalized (no nonmasked
    pixels in [1310, 1325] Angstroms QSO rest)
  - bit 3: 1 indicates less than 200 pixels available
  - bit 4: 1 indicates the ZWARNING: Warning for SDSS spectra. Usually
    refer to problems of z estimations. We keep MANY_OUTLIERS in our
    catalog but exclude all the other ZWARNING spectra.

  More than one bit may be set.
* `in_dr10`: binary flags indicating whether the observation exists in
  DR10Q; 1 indicates observation was included in DR10Q
* `in_dr9`: binary flags indicating whether the observation exists in
  DR9Q; 1 indicates observation was included in DR9Q
* `mjds`: SDSS modified Julian date of observation (DR12Q catalog,
  variable `MJD`)
* `plates`: SDSS spectroscopic plate (DR12Q catalog, variable `PLATE`)
* `ras`: right ascension in degrees (J2000) (DR12Q catalog, variable
  `RA`)
* `sdss_names`: SDSS-DR12 designation (DR12Q catalog, variable
  `SDSS_NAME`)
* `snrs`: median signal-to-noise ratio (whole spectrum) (DR12Q
  catalog, variable `SNR_SPEC`)
* `thing_ids`: SDSS unique identifier (DR12Q catalog, variable
  `THING_ID`)
* `z_qsos`: visual inspection QSO redshift from DR12Q (DR12Q catalog,
  variable `Z_VI`)

The `catalog.mat` file also contains data relating to the results of
previous DLA searches, contained in a series of hash tables. Each
table has three keys corresponding to the search:

* `dr9q_concordance`, containing the concordance catalog from the
  [BOSS DR9 Lyman-alpha Forest Catalog and Sample](https://www.sdss3.org/dr9/algorithms/lyaf_sample.php)
  ([Lee, et al. 2013](http://adsabs.harvard.edu/abs/2013AJ....145...69L)),
* `dr12q_visual`, containing the results of a visual survey (extracted
   from
   ([Noterdaeme DR12 catalog v2, 2016](http://www2.iap.fr/users/noterdae/DLA/DLA_DR12_v2.tgz))),
* `dr12q_noterdaeme`, containing the results of the method from
   ([Noterdaeme et al., 2012c](http://adsabs.harvard.edu/abs/2012A%26A...547L...1N)),
   available at
   ([Noterdaeme DR12 catalog v2, 2016](http://www2.iap.fr/users/noterdae/DLA/DLA_DR12_v2.tgz)).

The hash tables are:

* `los_inds`: Each key maps to a binary vector of length N. A value
  of 1 indicates the line of sight was considered in the corresponding
  search.
* `dla_inds`: Each key maps to a binary vector of length N. A value
  of 1 indicates the line of sight was considered in the corresponding
  search and flagged as containing a DLA.
* `z_dlas`: Each key maps to a map of length N. The entries for which
  `dla_inds(key)(qso_ind) == 0`, that is, the lines of sight either
  not considered by the DLA search or considered but not flagged as
  containing a DLA, are empty. The others contain a list of absorber
  redshifts for any DLAs flagged along the line of sight by the
  search. For the key `dr12q_visual`, the flagged DLAs have absorber
  redshift coded as the QSO redshift, as the visual inspection did not
  record the absorber redshift. In the case of lines of sight
  containing multiple flagged DLAs, the order of the objects matches
  the order in `log_nhis` below.
* `log_nhis`: Each key maps to a map of length N. The entries for
  which `dla_inds(key)(qso_ind) == 0`, that is, the lines of sight
  either not considered by the DLA search or considered but not
  flagged as containing a DLA, are empty. The others contain a list of
  decimal logarithmic column densities (in units cm^-2) for any DLAs
  flagged along the line of sight by the search. For the key
  `dr12q_visual`, the flagged DLAs have absorber redshift coded as
  having column density log N_HI = 20.3, as the visual inspection did
  not record the column density. In the case of lines of sight
  containing multiple flagged DLAs, the order of the objects matches
  the order in `z_dlas`.

dla_samples_a03.mat
---------------

Information about the DLA parameter samples used to evaluate the model
evidence of the DLA model. N = 10000 such samples were used. The
variables are:

* `alpha`: the fraction of column density samples to draw from the
  data-driven model rather than the uniform model (value: 0.97)
  Note: Garnett (2017) used a=0.90
* `fit_max_log_nhi`: the maximum log_10 N_HI  column density
  considered for fitting the data-driven column-density prior (value:
  22)
* `fit_min_log_nhi`: the maximum log_10 N_HI column density
  considered for fitting the data-driven column-density prior (value:
  20)
* `log_nhi_samples`: a length-N vector containing the log_10 N_HI
  column density values corresponding to each parameter sample
* `nhi_samples`: a length-N vector containing the N_HI column density
  values corresponding to each parameter sample (equal to
  10^`log_nhi_samples`)
* `offset_samples`: a length-N vector containing the absorber redshift
  offsets corresponding to each parameter sample; the offset is a
  value in (0, 1) that determines a unique value in [z_min, z_max] via
  a simple affine map
* `uniform_max_log_nhi`: the maximum log_10 N_HI column density used
  for the uniform component of the prior (value: 23)
* `uniform_min_log_nhi`: the minimum log_10 N_HI column density used
  for the uniform component of the prior (value: 20)

learned_qso_model_lyseries_variance_kim_dr9q_minus_concordance.mat
--------------------------------------------

The GP model learned for the null model of quasar emission.  This
is the learned model includes lyman series forest derived from
Kim (2017) in the mean model and the noise variance (Omega).
The variables are:

* `M`: (size 1217x20) the learned "eigenspectra" for the null model
* `initial_M`: (size 1217x20) the initial guess for `M` fed to the
  optimization routine
* `initial_log_sigma`: (size 1217x1) the initial guess for `log_sigma`
  fed to the optimization routine
* `log_sigma`: (size 1217x1) the learned log "absorption noise"
  for the null model
* `log_likelihood`: the log likelihood of the training data given the
  learned parameters (value: 2.883 x 10^7)
* `max_noise_variance`: the maximum normalized pixel noise variance
  considered during the training routine (value: 1)
* `minFunc_options`: the options passed into the `minFunc`
  optimization routine; values:

  - `MaxFunEvals`: 4000
  - `MaxIter`: 2000

* `minFunc_output`: the output from the `minFunc` optimization
  routine, containing:

  - `algorithm`: the algorithm used (value: 5, L-BFGS)
  - `firstorderopt`: the first-order optimality at termination (value:
    0.0597)
  - `iterations`: the number of iterations (value: 1107)
  - `funcCount`: the number of function evaluations (value: 1150)
  - `message` the exit message (string, value: `Step Size below progTol`)
  - `trace`: trace of optimization process per iteration:

    * `fval`: the best function value seen by the termination of each
      iteration
	* `funcCount`: the cumulative total number of function evaluations
      by the termination of each iteration
	* `optCond`: the first-order optimality at the termination of each
      iteration

* `mu`: (size 1217 x 1) the learned mean vector
* `rest_wavelengths`: (size 1217x1) the rest wavelengths in Angstroms
  on which the model was learned
* `train_ind`: (size 297301x1) a binary vector indicating which
  observations from the DR12Q catalog were considered for building the
  training data; a 1 indicates that object was considered. The value
  here corresponds to all unfiltered objects in DR12Q that were also
  in DR9Q and considered for inclusion in the DR9Q concordance
  catalog, but not in the DR9Q concordance catalog.
* `training_release`: the name of the release used for training
  (string, value: `dr12q`, corresponds to objects in `catalog.mat` and
  `preloaded_qsos.mat`)

preloaded_qsos.mat
------------------

Preloaded, truncated, and normalized QSO spectra. Each spectrum not
filtered from the search is loaded and:

* pixels with inverse noise variance set to 0 by the SDSS pipeline or
  with flag `BRIGHTSKY` set in the "and" mask are masked
* the median flux among nonmasked pixels in the range [1270, 1290]
  Angstroms QSO rest is computed as a normalizer
* the flux and pipeline noise is normalized by dividing by the computed
  normalizer
* pixels outside the range [910, 1217] Angstroms QSO rest are discarded

The results of this initial preprocessing are stored here as arrays of
length N = 297301. All variables corresponding to filtered
observations are empty. The variables are:

* `all_flux`: normalized flux
* `all_noise_variance`: normalized noise variance
* `all_normalizers`: the normalizer computed during preprocessing in
  10^-17 erg s^-1 cm^-2 A^-1
* `all_pixel_mask`: the liberal pixel mask used (`IVAR == 0 ||
  BRIGHTSKY`)
* `all_wavelengths`: the observed wavelengths corresponding to each
  pixel in Angstroms

Additionally, several loading parameters are stored:

* `loading_max_lambda`: the maximum rest wavelength in Angstroms to
  load (value: 1217)
* `loading_min_lambda`: the minimum rest wavelength in Angstroms to
  load (value: 910)
* `min_num_pixels`: the minimum number of nonmasked pixels required to
  retain a spectrum (value: 200)
* `normalization_max_lambda`: the maximum rest wavelength in Angstroms
  considered for flux normalization (value: 1290)
* `normalization_min_lambda`: the minimum rest wavelength in Angstroms
  considered for flux normalization (value: 1270)

processed_qsos_multi_lyseries_a03_lyb_zwarn_occams_dr12q.mat
------------------------

The results of the DLA search on the 158825 nonfiltered lines
of sight (`filter_flags == 0` in `catalog.mat`). Note: we only
include the rest-frame wavelengths between Lyman beta to Lyman
alpha. The variables include:

* `all_exceptions`: all the spectra have the interpolation error, usually
  are the spectra which are empty (with no flux/lambda values)
* `base_sample_inds`: the inds used to Monte carlo resample the DLA samples
  based on the parameter priors. For a given QSO, the first row stored the
  inds we used to sample the DLA parameters for M_DLA2; the second row stored
  the  inds we used to sample the parameters for M_DLA3. For example, if we
  want to sample for M_DLA3, we would use `sample_z_dlas` and `nhi_samples` to
  sample one of the DLAs, and use the first and second rows in
  `base_sample_inds(quasar_ind,:,:)` to sample the other two DLAs.
  Here `quasar_ind` is the index of the QSO we sample.
* `dla_catalog_name`: the DLA catalog used to compute the DLA model
  prior Pr(M_DLA | z_QSO) (string, value: `dr9q_concordance`; see
  corresponding DLA catalog information in `catalog.mat`)
* `log_likelihoods_dla`: log p(D | M_DLA, z_QSO), it includes M_DLA1 to
  M_DLA4 in the form:
  [log p(D | M_DLA1, z_QSO), log p(D | M_DLA2, z_QSO), log p(D | M_DLA3, z_QSO)
  log p(D | M_DLA4, z_QSO)]
* `log_likelihoods_no_dla`: log p(D | M_!DLA, z_QSO)
* `log_likelihoods_lls`: log p(D | M_subDLA, z_QSO)
* `log_posteriors_dla`: log Pr(M_DLA | z_QSO) + log p(D | M_DLA, z_QSO)
* `log_posteriors_no_dla`: log Pr(M_!DLA | z_QSO) + log p(D | M_!DLA, z_QSO)
* `log_posteriors_lls`: log Pr(M_subDLA | z_QSO) + log p(D | M_subDLA, z_QSO)
* `log_priors_dla`: log Pr(M_DLA | z_QSO)
* `log_priors_no_dla`: log Pr(M_!DLA | z_QSO)
* `log_priors_lls`: log Pr(M_subDLA | z_QSO)
* `MAP_z_dlas`: the maximum a posteriori (MAP) predictions of the redshift
  of DLAs for a given M_DLA(k). For each QSO, we store MAP predictions from
  M_DLA1 to M_DLA4
* `MAP_log_nhis`: the maximum a posteriori (MAP) predictions of the log
  column densities of DLAs for a given M_DLA(k). For each QSO, we store
  MAP predictions from M_DLA1 to M_DLA4
* `max_z_cut`: the maximum absorber redshift considered is (z_QSO -
  max_z_cut) < (z_QSO + max_z_cut) (value: 3000 km s^1)
* `max_z_dlas`: the maximum absorber redshift considered
* `min_z_dlas`: the minimum absorber redshift considered
* `model_posteriors`: Pr(M | D, z_QSO); first column corresponds to
  M_!DLA, second to M_subDLA, third is M_DLA1, fourth is M_DLA2,
  fifth is M_DLA3, sixth is M_DLA4; M_!DLA is `p_no_dlas`, M_subDLA
  is `p_lls`, and the summation of M_DLA1 to M_DLA4 is `p_dla`
* `num_lines`: the number of members in the Lyman series to consider
  when computing Voigt profiles corresponding to DLA parameters
  (value: 3, consider Lyman-alpha, -beta, and -gamma absorption)
* `p_dlas`: Pr(M_DLA | D, z_QSO)
* `p_no_dlas`: Pr(M_!DLA | D, z_QSO)
* `p_lls`: Pr(M_subDLA | D, z_QSO)
* `prior_ind`: a binary vector indicating the quasars used to compute
  the model prior Pr(M | z_QSO) (value: all spectra in DR9Q considered
  for inclusion in the DR9Q concordance catalog)
* `prior_z_qso_increase`: the model prior is computed using all
  objects in the training catalog with z < (z_QSO +
  `prior_z_qso_increase`) (value: 30000 km s^-1)
* `release`: the data release considered (string, value: `dr12q`)
* `sample_log_likelihoods_dla`: the log likelihoods for each DLA
  parameter sample for each object; log p(D | M_DLA, z_QSO, theta_i);
  see DLA parameter sample values in `dla_samples_a03.mat`. Note that the
  absorber redshift for object i in parameter sample j is computed as
  `min_z_dlas`(i) + (`max_z_dlas`(i) - `min_z_dlas`(i)) * `offset_samples`(j);
  Note: we store sample likelihoods from M_DLA1 to M_DLA4
* `sample_log_likelihoods_lls`: the log likelihoods for each subDLA
  parameter sample for each object; log p(D | M_subDLA, z_QSO, theta_i);
  see DLA parameter sample values in `set_lls_parameters.m`. Note that the
  absorber redshift for object i in parameter sample j is computed as
  `min_z_dlas`(i) + (`max_z_dlas`(i) - `min_z_dlas`(i)) * `offset_samples`(j)
* `test_ind`: binary flags into the DR12Q catalog indicating which
  lines of sight were considered (equivalent to `filter_flags == 0`
  for the objects in `catalog.mat`)
* `test_set_name`: the name of the test set (string, value: `dr12q`)
* `training_release`: the name of the training release used to compute
  the null model parameters p(D | M_DLA, z_QSO) (string, value: `dr12q`)
* `training_set_name`: the name of the training set used used to
  compute the null model likelihoods p(D | M_DLA, z_QSO) (string,
  value: `dr9q_minus_concordance`)

processed_qsos_multi_lyseries_a03_lyb_zwarn_occams_trunc_dr12q.mat
------------------------

Same as processed_qsos_multi_lyseries_a03_lyb_zwarn_occams_dr12q.mat, but
without `sample_log_likelihoods_dla`, `sample_log_likelihoods_lls`, and
`base_sample_inds` these large matrices.

processed_qsos_multi_lyseries_a03_zwarn_occams_trunc_dr12q.mat
------------------------

Same as processed_qsos_multi_lyseries_a03_lyb_zwarn_occams_dr12q.mat, but
without `sample_log_likelihoods_dla`, `sample_log_likelihoods_lls`, and
`base_sample_inds` these large matrices. Also, it samples the DLA parameters
from Lyman limit to Lyman alpha.

================================================================================
(End)          Roman Garnett [Washington University in St. Louis]           2017
(Update)       Ming-Feng  Ho [University of California, Riverside]          2020
