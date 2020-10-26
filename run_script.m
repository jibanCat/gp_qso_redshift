c4_set_parameters
c4_build_catalog
c4_preload_qsos
% load('data/dr7/processed/preloaded_qsos.mat');
cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
filter_flags= f.filter_flags;

half_ID = randsample(all_QSO_ID, int32(0.7*numel(all_QSO_ID)));
test_ind = ((~ismember(all_QSO_ID, half_ID)) & (filter_flags==0));

prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
               (ismember(all_QSO_ID, half_ID)));  
train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
                ismember(all_QSO_ID, half_ID) );
save('indeces-70%.mat', 'train_ind', 'prior_ind','test_ind')

learn_qso_model
generate_c4_samples
process_qsos
% generate_ascii_catalog