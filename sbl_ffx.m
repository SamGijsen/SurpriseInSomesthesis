function sbl_ffx(sub, model_class)
% This function performs FFX BMS on EEG electrode data.It relies on 
% spm_vi_glm.m script to approximate the log model evidence using a 
% variational inference algorithm.
%
%   Inputs
%       sub         : number of subject data file to be analyzed
%       model_class : 'DC' or 'HMM', specifying which model class the regressors 
%       belong to. In case of 'DC', it computes free energy values for 101
%       values of \tau.
%
%   Output
%       lme         : output matrix containing the free energy values for each
%       electrode, peri-stimulus timebin, model, and, if applicable,
%       value of the forgetting parameter \tau
% -------------------------------------------------------------------------
% Written by: Dirk Ostwald, Miro Grundei, Sam Gijsen


% directory prep
% -------------------------------------------------------------------------
addpath('')                                                                 ; % add path to spm_vi_glm.,
ddir        = fullfile('')                                                  ; % EEG data directory
resdir      = fullfile('')                                                  ; % results directory
btdir       = fullfile('')                                                  ; % EEG badtrials directory
regdir      = fullfile('')                                                  ; % model regressor directory
verbose     = 0                                                             ; % verbosity

% data analysis parameters
% -------------------------------------------------------------------------
n_t         = 359                                                           ; % number of peristimulus time bins
z_score     = 1                                                             ; % regressor zscoring flag
n_i         = 20                                                            ; % maximal number of VI algorithm iterations
delta       = 1e-5                                                          ; % variational free energy convergence criterion
alpha       = 1e3                                                           ; % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
beta_lambda = 1e1                                                           ; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
gama_lambda = 1e-1                                                          ; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
n_e         = 64                                                            ; % number of electrodes of interest

% check which model class is supplied
if strcmp(model_class, 'DC') 
    taus    = [0, logspace(-4,0,100)]                                       ; % values of forgetting parameter (\tau) used to previously generate regresors
    ntaus   = size(taus,2)                                                  ;
elseif strcmp(model_class, 'HMM')
    ntaus   = 1                                                             ; % \tau is not applicable for the HMM
else
    error('Incorrect "model_class" variable: DC or HMM')
end

% data analysis array initializations
% -------------------------------------------------------------------------
plab        = {'SP', 'AP', 'TP'}                                            ; % probability type / inference (SP, AP, TP)
n_p         = numel(plab)                                                   ;
mlab        = {'PS', 'BS', 'CS'}                                            ; % surprise model labels
mlab_full   = {'predictive_surprise','bayesian_surprise',                   ...
               'confidence_corrected_surprise'}                             ;
n_m         = 3                                                             ; % number of models of interest
X           = cell(1,n_m)                                                   ; % design matrix
n_trl       = NaN(1)                                                        ; % number of trials upon bad and catch trial removal
lme         = NaN(ntaus, n_t, n_m, n_p, n_e)                                ; % log model evidence array

% data preprocessing
% ---------------------------------------------------------------------
% participant-specific data filenames
dfile       = fullfile(ddir, sprintf('MAT_bletdhmar_%s_SBL.mat', sub))      ; % participant EEG file
bfile       = fullfile(btdir, sprintf('badtrials_%s_SBL.mat',sub))          ; % participant bad trials file

% electrode data loading
S           = load(dfile, 'Y', 'events', 'chanlabels', 'time')              ; % load EEG data file containing data (Y)
events      = S.events                                                      ; % EEG trigger values used to discern identity of stimuli

catch_trials = find(events == 33)                                           ; % catch-trials are coded with a '33'
load(bfile, 'badtrials')                                                    ; % bad trial (EEG) indices
Y           = S.Y                                                           ; % n_d x n_t x number of trials data array
rmv         = unique(sort([catch_trials, badtrials]))                       ; % indices of trials to be removed
n_trl(s)    = size(Y,3) - length(rmv)                                       ; % number of valid trials
Y(:,:,rmv)  = []                                                            ; % electrode data trial removal (both bad and catch trials)  

clear S

% electrode iterations
% ---------------------------------------------------------------------
for e = 1:64

    % user information
    if verbose
        fprintf('Processing electrode %e of %e\n', e, n_e) 
    end

    % probability type iterations
    % -----------------------------------------------------------------
    for p = 1:n_p

        % tau iterations
        % -------------------------------------------------------------
        for tau = 1:size(taus,2)

            % model preprocessing
            % -------------------------------------------------------------   
            x = zeros(4000,numel(mlab));            
            for m = 1:n_m  
                if strcmp(model_class, 'DC')                                                     % in case of DC model regressors, the tau value is needed for the correct filename
                    tau_label   = sprintf('%s_tau_%.10f_%s_CD.mat', sub,taus(tau),plab{p});      % construct part of regressor filename
                    mfile      	= fullfile(regdir, sprintf('sub-%02d', str2num(sub_number)), tau_label) ; % participant DC regressor file
                else
                    mfile      	= fullfile(regdir, sprintf('sub-%02d_%s_HMM.mat', sub, plab{p})); % participant HMM regressor file
                end

                reg             = load(mfile)                               ; % model regressor structure                                                                
                x(:,m)          = reg.(mlab_full{m})'                       ; % regressors of interest              
            end

            % remove bad/catch-trials
            x(rmv,:) = [];                      

            % regressor z scoring
            if z_score
                x = zscore(x);
            end                                        

            % model formulation & save DM
            n           = size(x,1)                                        ; % number of data points 
            DM          = cell(1,4)                                        ; % initialization
            DM{1}       = ones(n,1)                                        ; % null model design matrix
            DM{2}       = [ones(n,1) x(:,1)]                               ; % predictive suprise model design matrix
            DM{3}       = [ones(n,1) x(:,2)]                               ; % bayesian suprise model design matrix
            DM{4}       = [ones(n,1) x(:,3)]                               ; % confidence corrected suprise model design matrix       
            X{s,p}      = DM;             

            % peri-stimulus time bin iterations
            % -------------------------------------------------------------
            for t = 1:n_t 

                % analysis model iterations 
                % ---------------------------------------------------------

                num_mdl = size(X{1,p},2);
                for mdl = 1:num_mdl

                    % spm_vb_glm analysis
                    glm             = []                                  	; % structure initialization
                    glm.X           = X{1,p}{mdl}                         	; % design matrix
                    glm.y           = zscore(squeeze(Y(e,t,:)))           	; % data (constant over models)
                    glm.n           = size(glm.X,1)                       	; % number of data points
                    glm.p           = size(glm.X,2)                        	; % number of analysis model regression parameters
                    glm.mu_beta     = zeros(glm.p,1)                     	; % expectation parameter of p(\beta) = N(\beta;0_p,\alpha^{-1}I_p)
                    glm.alpha       = alpha                              	; % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
                    glm.beta_lambda = beta_lambda                         	; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                    glm.gama_lambda = gama_lambda                         	; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                    glm.p           = size(glm.X,2)                        	; % number of beta parameters
                    glm.n_i         = n_i                                 	; % maximum number of iterations
                    glm.delta       = delta                               	; % variational free energy convergence criterion
                    glm             = spm_vi_glm(glm)                     	; % estimation

                    % record converged variational free energy only
                    lme(tau,t,mdl,p,e)    = glm.F_max; 

                end
            end
        end
    end
end

% set filename and save
fn = sprintf('%sLME_%s_sub-%02d.mat', resdir, model_class, str2num(sub_number)); 
save(fn, 'lme', '-v7.3')
end
