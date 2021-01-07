function sbl_rfx(electrode_num) 
% This function performs RFX BMS using previously computed free energy values. 
% As 'electrode_num' corresponds to broadcasted index of electrode [1:64] to 
% be analyzed, the function thus may be easily parallelized. It relies on 
% the SPM12 toolbox to perform the BMS inference.

% Performs the following analyses:
% A1) family RFX: DC model vs HMM vs null model 
% A2) family RFX: TP1 vs TP2 inference within DC model family
% A3) family RFX: SP vs AP vs TP inference within DC model family
% A4) model RFX : PS vs BS vs CS computation within DC TP1 model family
%
%   Inputs
%       electrode_number : number of electrode to be analyzed (1-64)
%       The function furthermore loads in a matrix containing free energy
%       values for all models, subjects, timebins, and electrodes.
%
%   Output
%       Four files containing the exdeence probabilities associated with
%       the aforementioned four analyses (A1-A4).
%
% -------------------------------------------------------------------------
% Written by: Sam Gijsen, Miro Grundei

% directory preparation
% -------------------------------------------------------------------------

addpath('')                                                                 ; % add path to SPM12
ffxdir      = fullfile('')                                                  ; % directory containing results from FFX analysis directory
resdir      = fullfile('')                                                  ; % results directory
verbose     = 0                                                             ; % verbosity

%% A1) family RFX: DC model vs HMM vs null model 
% -------------------------------------------------------------------------
elec    = str2num(electrode_num)                                            ; % electrode to be analyzed
n_f     = 3                                                                 ; % number of families (Null, DC, HMM)

% load previously approximated log-model evidence values for DC model
load([ffxdir, 'FFX_results.mat'], 'opt_F_elec_DC') 
opt_F_elec_DC = test;
F_DC               = squeeze(opt_F_elec_DC(:,:,:,elec,:))                   ; % select free energy values only for the relevant electrode

% load previously approximated log-model evidence values for HMM 
load([ffxdir, 'FFX_results.mat'], 'opt_F_elec_HMM') 
opt_F_elec_HMM = test;
F_HMM               = squeeze(opt_F_elec_HMM(:,:,:,elec,:))                 ; % select free energy values only for the relevant electrode

n_m                 = 12+12+1                                               ; % number of total models considered (12 for DC, 12 for HMM, 1 null)
n_e                 = 40                                                    ; % number of subjects
n_t                 = 359                                                   ; % number of peri-stimulus timebins
lme_fam             = nan(n_t,n_e,n_m)                                      ; % initialization of log-model evidence values by families
m_count             = 1                                                     ; % initialize count

% collect lme values in a single matrix
% DC model
for p = 1:4     % SP, AP, TP1, TP2
    for m = 2:4 % PS, BS, CS
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_fam(:,:,m_count)    = squeeze(F_DC(:,m,p,:))                    ; % select appropriate model regressor
    end
end

% HMM model
for p = 1:4
    for m = 2:4
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_fam(:,:,m_count)    = squeeze(F_HMM(:,m,p,:))                   ; % select appropriate model regressor
    end
end

% null model
lme_fam(:,:,1)      = squeeze(F_DC(:,1,1,:))                                ; % null-model comes first

% set RFX BMS parameters
family.infer        = 'RFX'                                                 ; 
family.partition	= [1, 2*ones(1,12), 3*ones(1,12)]                       ; 
family.names        = {'NULL', 'DC', 'HMM'}                                 ;
family.Nsamp     	= 5e5                                                   ; % amount of samples for RFX sampling procedure

fam_xp          	= NaN(n_t, n_f)                                         ; % initialize results matrix for family exceedence probabilities
mod_xp          	= NaN(n_t, n_m)                                         ; % initialize results matrix for model exceedence probabilities

for t = 1:n_t % loop over peri-stimulus timebins
   [f, m]         	= spm_compare_families(squeeze(lme_fam(t,:,:)), family) ; % perform family-level RFX
   fam_xp(t,:)    	= f.xp                                                  ; % family xp
   mod_xp(t,:)      = m.xp                                                  ; % model-specific xp
end

% set filename and save results
fn                  = sprintf('%sRFX_level1_elec-%02d.mat', resdir, str2num(electrode_num)); 
save(fn, 'fam_xp', 'mod_xp', '-v7.3')
clear fam_xp mod_xp F_HMM


%% A2) family RFX: TP1 vs TP2 inference within DC model family
% -------------------------------------------------------------------------
n_f                 = 2                                                   	; % number of families (Null, DC, HMM)
n_m                 = 3+3                                                   ; % number of total models considered (3 for TP1, 3 for TP2)
n_e                 = 40                                                    ; % number of subjects
n_t                 = 359                                                   ; % number of peri-stimulus timebins
lme_fam             = nan(n_t,n_e,n_m)                                      ; % initialization of log-model evidence values by families
m_count             = 0                                                     ; % initialize count

% collect lme values in a single matrix
for p = 3:4
    for m = 2:4
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_fam(:,:,m_count)    = squeeze(F_DC(:,m,p,:))                    ; % select appropriate model regressor
    end
end

% set RFX BMS parameters
family.infer        = 'RFX'                                                 ; 
family.partition	= [1*ones(1,3), 2*ones(1,3)]                            ; 
family.names        = {'TP1', 'TP2'}                                        ;
family.Nsamp     	= 5e5                                                   ; % amount of samples for RFX sampling procedure

fam_xp          	= NaN(n_t, n_f)                                         ; % initialize results matrix for family exceedence probabilities
mod_xp          	= NaN(n_t, n_m)                                         ; % initialize results matrix for model exceedence probabilities

for t = 1:n_t % loop over peri-stimulus timebins
   [f, m]         	= spm_compare_families(squeeze(lme_fam(t,:,:)), family) ; % perform family-level RFX
   fam_xp(t,:)    	= f.xp                                                  ; % family xp
   mod_xp(t,:)      = m.xp                                                  ; % model-specific xp
end

% set filename and save results
fn                  = sprintf('%sRFX_level2_elec-%02d.mat', resdir, str2num(electrode_num)); 
save(fn, 'fam_xp', 'mod_xp', '-v7.3')
clear fam_xp mod_xp


%% A3) family RFX: SP vs AP vs TP inference within DC model family
% -------------------------------------------------------------------------
n_f                 = 3                                                     ; % number of families (SP, AP, TP)
n_m                 = 3+3+3                                                 ; % number of total models considered (3 for SP, 3 for AP, 3 for TP)
n_e                 = 40                                                    ; % number of subjects
n_t                 = 359                                                   ; % number of peri-stimulus timebins
lme_fam             = nan(n_t,n_e,n_m)                                      ; % initialization of log-model evidence values by families
m_count             = 0                                                     ; % initialize count

% collect lme values in a single matrix
% DC model
for p = 1:3
    for m = 2:4
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_fam(:,:,m_count)    = squeeze(F_DC(:,m,p,:))                    ; % select appropriate model regressor
    end
end

% set RFX BMS parameters
family.infer        = 'RFX'                                                 ; 
family.partition	= [1*ones(1,3), 2*ones(1,3), 3*ones(1,3)]               ; 
family.names        = {'SP', 'AP', 'TP'}                                    ;
family.Nsamp     	= 5e5                                                   ; % amount of samples for RFX sampling procedure

fam_xp          	= NaN(n_t, n_f)                                         ; % initialize results matrix for family exceedence probabilities
mod_xp          	= NaN(n_t, n_m)                                         ; % initialize results matrix for model exceedence probabilities

for t = 1:n_t % loop over peri-stimulus timebins
   [f, m]         	= spm_compare_families(squeeze(lme_fam(t,:,:)), family) ; % perform family-level RFX
   fam_xp(t,:)    	= f.xp                                                  ; % family xp
   mod_xp(t,:)      = m.xp                                                  ; % model-specific xp
end

% set filename and save results
fn                  = sprintf('%sRFX_level3_elec-%02d.mat', resdir, str2num(electrode_num)); 
save(fn, 'fam_xp', 'mod_xp', '-v7.3')
clear fam_xp mod_xp


%% A4) model RFX : PS vs BS vs CS computation within DC TP1 model family
% -------------------------------------------------------------------------
n_f                 = 3                                                     ; % number of families (SP, AP, TP)
n_m                 = 3+3+3                                                 ; % number of total models considered (3 for SP, 3 for AP, 3 for TP)
n_e                 = 40                                                    ; % number of subjects
n_t                 = 359                                                   ; % number of peri-stimulus timebins
lme_fam             = nan(n_t,n_e,n_m)                                      ; % initialization of log-model evidence values by families
m_count             = 0                                                     ; % initialize count

% collect lme values in a single matrix
% DC model
for p = 3
    for m = 2:4
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_fam(:,:,m_count)    = squeeze(F_DC(:,m,p,:))                    ; % select appropriate model regressor
    end
end
% sampling - use sampling to compute exact alpha
% ecp      - 1 to compute exceedance probability
% set RFX BMS parameters
family.infer        = 'RFX'                                                 ; 
family.partition	= [1*ones(1,3), 2*ones(1,3), 3*ones(1,3)]               ; 
family.names        = {'PS', 'BS', 'CS'}                                    ;
family.Nsamp     	= 5e5                                                   ; % amount of samples for RFX sampling procedure

fam_xp          	= NaN(n_t, n_f)                                         ; % initialize results matrix for family exceedence probabilities
mod_xp          	= NaN(n_t, n_m)                                         ; % initialize results matrix for model exceedence probabilities

for t = 1:n_t % loop over peri-stimulus timebins
   [f, m]         	= spm_compare_families(squeeze(lme_fam(t,:,:)), family) ; % perform family-level RFX
   fam_xp(t,:)    	= f.xp                                                  ; % family xp
   mod_xp(t,:)      = m.xp                                                  ; % model-specific xp
end

% set RFX BMS parameters, this time no family-level inference is necessary
Nsamp               = 5e5                                                   ; % amount of samples
do_plot             = 0                                                     ; % plotting
sampling            = 0                                                     ; % use sampling to compute exact alpha
ecp                 = 1                                                     ; % 1 to compute exceedance probability

n_m                 = 3                                                     ; % number of total models considered (PS, BS, CS)
n_e                 = 40                                                    ; % number of subjects
n_t                 = 359                                                   ; % number of peri-stimulus timebins
lme_res             = nan(n_t,n_e,n_m)                                      ; % initialization of log-model evidence values
m_count             = 0                                                     ; % initialize count

% collect lme values in a single matrix
% DC TP model
for p = 3
    for m = 2:4
        m_count                 = m_count+1                                 ; % increment count by 1
        lme_res(:,:,m_count)    = squeeze(F(:,m,p,:))                       ; % select appropriate model regressor
    end
end

TP_mod_xp           = NaN(n_t, n_m)                                         ; % initialize results matrix for model exceedence probabilities

for t = 1:n_t % loop over peri-stimulus timebins
   [~,~,xp,~,~]     = spm_BMS(squeeze(lme_res(t,:,:)), Nsamp, do_plot, sampling, ecp); % perform RFX BMS  
   TP_mod_xp(t,:)  	= xp                                                    ; % store exceedence probabilities
end

% set filename and save results
fn = sprintf('%sRFX_level4_elec-%02d.mat', resdir, str2num(electrode_num))  ; 
save(fn, 'TP_mod_xp', '-v7.3')

end % function end