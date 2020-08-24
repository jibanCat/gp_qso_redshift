% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
% load the DR14Q catalog
variables_to_load = {'z_qsos', 'plates', 'mjds', 'fiber_ids', 'filter_flags'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

num_quasars = numel(z_qsos);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);

for i = 1:num_quasars
  if (filter_flags(i) > 0)
    continue;
  end

  try
    [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask] ...
        = file_loader(plates(i), mjds(i), fiber_ids(i));
  catch
    warning('File non-existed on the SDSS website.');

    % bit 5: non-existed file
    filter_flags(i) = bitset(filter_flags(i), 6, true);
    continue;
  end

  % do not normalize flux: this is done in the learning and processing code.

  ind = (this_wavelengths >= (min_lambda * (z_qso_cut) + 1)) & ...
        (this_wavelengths <= (max_lambda * (z_qso_training_max_cut) + 1)) & ...
        (~this_pixel_mask);

  % bit 3: not enough pixels available
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_wavelengths{i}    =    this_wavelengths;%no longer (ind)
  all_flux{i}           =           this_flux;
  all_noise_variance{i} = this_noise_variance;
  all_pixel_mask{i}     =     this_pixel_mask;

  fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
          i, num_quasars, plates(i), mjds(i), fiber_ids(i));
end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers'};
save(sprintf('%s/preloaded_zqso_only_qsos', processed_directory(release)), ...
     variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
save(sprintf('%s/zqso_only_catalog', processed_directory(release)), ...
     'filter_flags', '-append');
