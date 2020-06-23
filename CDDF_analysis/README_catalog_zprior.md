zqso_only_catalog.mat
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

The `zqso_only_catalog.mat` file also contains data relating to the results of
previous DLA searches, contained in a series of hash tables. Each
table has three keys corresponding to the search:
(Note: the code provided on GitHub: https://github.com/sbird/gp_qso_redshift
does not include DLA catalogue. The catalog file here runs on both zestiamtion
and DLA finding, that's the discrepancy.)

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

learned_zqso_only_model_outdata_full_dr9q_minus_concordance_norm_1176-1256.mat
--------------------------------------------

The GP model learned for the null model of quasar emission.

* `M`: (size 8361x20) the learned "eigenspectra" for the null model
* `initial_M`: (size 8361x20) the initial guess for `M` fed to the
  optimization routine
* `bluewards_mu`: the learned blueward mean value (assuming blueward flux is
    drawing from an iid Guassian). (value: 0.162)
* `bluewards_sigma`: the learned blueward standard deviation (assuming blueward
    flux is drawing from an iid Guassian). (value: 0.047)
* `log_likelihood`: the log likelihood of the training data given the
  learned parameters (value: -2.398 x 10^8)
* `max_noise_variance`: the maximum normalized pixel noise variance
  considered during the training routine (value: 16)
* `minFunc_options`: the options passed into the `minFunc`
  optimization routine; values:

  - `MaxFunEvals`: 8000
  - `MaxIter`: 4000

* `minFunc_output`: the output from the `minFunc` optimization
  routine, containing:

  - `algorithm`: the algorithm used (value: 5, L-BFGS)
  - `firstorderopt`: the first-order optimality at termination (value:
    1.17 x 10^5)
  - `iterations`: the number of iterations (value: 3)
  - `funcCount`: the number of function evaluations (value: 32)
  - `message` the exit message (string, value: `Step Size below progTol`)
  - `trace`: trace of optimization process per iteration:

    * `fval`: the best function value seen by the termination of each
      iteration
	* `funcCount`: the cumulative total number of function evaluations
      by the termination of each iteration
	* `optCond`: the first-order optimality at the termination of each
      iteration

* `mu`: (size 8361x1) the learned mean vector
* `rest_wavelengths`: (size 8361x1) the rest wavelengths in Angstroms
  on which the model was learned
* `redwards_mu`: the learned redward mean value (assuming redward flux is
    drawing from an iid Guassian). (value: 0.1281)
* `redwards_sigma`: the learned redward standard deviation (assuming redward
    flux is drawing from an iid Guassian). (value: 0.0028)
* `train_ind`: (size 297301x1) a binary vector indicating which
  observations from the DR12Q catalog were considered for building the
  training data; a 1 indicates that object was considered. The value
  here corresponds to all unfiltered objects in DR12Q that were also
  in DR9Q and considered for inclusion in the DR9Q concordance
  catalog, but not in the DR9Q concordance catalog.
  (Note: the code provided on GitHub: https://github.com/sbird/gp_qso_redshift
does not include DLA catalogue. So the code on Github train on all unfiltered
spectra with or without DLAs. The train_ind here runs the same on both
zestiamtion and DLA finding, that's the discrepancy. The zestimation result
ideally would not change too much)
* `training_release`: the name of the release used for training
  (string, value: `dr12q`, corresponds to objects in `zqso_only_catalog.mat` and
  `preloaded_zqso_only_qsos.mat`)

processed_zqso_only_qsos_dr12q-100.mat
------------------

We do not truncate or normalize QSO spectra during preloading. But we built the
pixel_mask based on:

* pixels with inverse noise variance set to 0 by the SDSS pipeline or
  with flag `BRIGHTSKY` set in the "and" mask are masked

The results of this initial preprocessing are stored here as arrays of
length N = 297301. All variables corresponding to filtered
observations are empty. The variables are:

* `all_flux`: normalized flux
* `all_noise_variance`: normalized noise variance
* `all_normalizers`: we don't normalize the flux during preload. So all values
    are zeros
* `all_pixel_mask`: the liberal pixel mask used (`IVAR == 0 ||
  BRIGHTSKY`)
* `all_wavelengths`: the observed wavelengths corresponding to each
  pixel in Angstroms

Additionally, several loading parameters are stored:

* `loading_max_lambda`: This is not used in the script anymore since we don't
    normalize during preloading.
* `loading_min_lambda`: This is not used in the script anymore since we don't
    normalize during preloading.
* `min_num_pixels`: the minimum number of nonmasked pixels required to
  retain a spectrum (value: 400)
* `normalization_max_lambda`: This is not used in the script anymore since we don't
    normalize during preloading.
* `normalization_min_lambda`: This is not used in the script anymore since we don't
    normalize during preloading.

processed_zqso_only_qsos_dr12q-100.mat
------------------------

The results of the redshift estimations on the 162861 (158821 after removing NaN
results) nonfiltered lines of sight (`filter_flags == 0` in
`zqso_only_catalog.mat`).
Note: the test_ind here is a little bit different from the uniform prior run,
due to data-driven prior run was an earlier version. The result in the paper
used the uniform prior run.

* `all_exceptions`: all the empty spectra in `preloaded_qsos.mat`.
* `all_posdeferrors`: all spectra encountered MATLAB's positive definite errors
    during log likelihood calculation of GP. Usually due to normalization range
    is outside the scope of modelling window or very noisy instrumental noise.
* `offset_samples_qso`: all the samples of zQSO used to calculate log
    likelihoods.
  * `sample_log_posteriors`: the log posteriors for each zQSO
  parameter sample for each object; log p(M | D, z_QSO_i);
* `test_ind`: binary flags into the DR12Q catalog indicating which
  lines of sight were considered (equivalent to `filter_flags == 0`
  for the objects in `zqso_only_catalog.mat`)
* `test_set_name`: the name of the test set (string, value: `dr12q`)
* `training_release`: the name of the training release used to compute
  the null model parameters p(D | M_DLA, z_QSO) (string, value: `dr12q`)
* `training_set_name`: the name of the training set used used to
  compute the null model likelihoods p(D | M_DLA, z_QSO) (string,
  value: `dr9q_minus_concordance`)
* `z_map`: the maximum a posteriori estimation of zQSO.
* `z_qsos`: catalogue values of zQSOs.
* `z_true`: same as `z_qsos`.

hist2d_z_map_vs_z_true_pure-z-log.pdf
----

A 2D histogram indicating the predictions versus ground truth.

pdf_dr12q_misfits
----

A folder of the plots of misfit spectra (|z_map - z_true| > 0.5).

================================================================================
(End)          Roman Garnett [Washington University in St. Louis]           2017
(Update)       Ming-Feng  Ho [University of California, Riverside]      Mar 2020
(Update)       Ming-Feng  Ho [University of California, Riverside]      Jun 2020
