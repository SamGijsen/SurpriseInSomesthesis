function SBL_modelrecovery_ffx(val)

% This function performs simulations to validate the log marginal 
% likelihood evaluation of sequential Bayesian learning regressors
% in light of EEG single trial data. In essence, it evaluates a free-
% form variational inference algorithm for general linear models with
% spherical covariance matrices as originally proposed in the context of 
% general linear models with autoregressive errors by Penny et al. 
% 'Variational Bayesian inference for fMRI time series' Neuroimage 19 
% (2003) 727-741. 
%
% The algorthm is used here to perform a model recovery study given fixed
% \beta and \alpha parameters, but varying values for \lambda (i.e., under
% various levels of SNR.) Set up to be parallelized; here 'val' codes for
% the iteration. In our paper, 100 iterations were ran each with 40 fake
% subjects.
%	
%   Authors - Miro Grundei, Sam Gijsen
%   VI algorithm Author - Dirk Ostwald
% -------------------------------------------------------------------------

% initialization
% -------------------------------------------------------------------------
resdir      = fullfile('/scratch/samgijsen/results/F/model_recovery/')      ; % results directory
seqdir      = fullfile('')                                                  ; % observation-sequence directory
regdir      = fullfile('')                                                  ; % readout regressor directory

% intialize random number generator then diversify for each job/run
rng('default')
rng(str2num(val))

% subject number so we can later load and  use the actual sequences
% administered to participants
SJs             = [1:17,20,23:44]                                           ; % subjects

              
% model formulation 
% -------------------------------------------------------------------------
n_m         = 25                                                            ; % number of models
alpha       = 1e3                                                           ; % % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
lambdas     = 1./[1 10 100 500 750 1000 1250 1500 2000 3500 5000 10000]     ; % true, but unknown, observation noise parameter
n_l         = numel(lambdas)                                                ; % number of noise levels
beta_tbu    = [1;1]                                                         ; % true, but unknown, beta parameter vector
beta_lambda = 1e1                                                           ; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
gama_lambda = 1e-1                                                          ; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda)


% simulation parameters
% -------------------------------------------------------------------------
n_s         = 40                                                            ; % number of simulations per parameter setting
n_i         = 20                                                            ; % maximum number of variational inference algorithm
delta       = 1e-5                                                          ; % variational free energy convergence criterion
S           = cell(n_l,n_s,n_m,n_m)                                         ; % simulation result array initialization
lme         = NaN(n_s, n_m, n_m, n_l)                                       ; % log-model evidence array initialization
beta        = NaN(n_s, n_m, n_m, n_l)                                       ; % \beta parameter array initialization
lambda      = NaN(n_s, n_m, n_m, n_l)                                       ; % \lambda parameter array initialization

% Model Compendium
% > ------------------------------------------------------------- <
% 1. Null model
% 2.  DC SP   PS tau=0.007 % 14. HMM SP  PS 
% 3.  DC SP   BS tau=0.007 % 15. HMM SP  BS 
% 4.  DC SP   CS tau=0.007 % 16. HMM SP  CS 
% 5.  DC AP   PS tau=0.007 % 17. HMM AP  PS 
% 6.  DC AP   BS tau=0.007 % 18. HMM AP  BS 
% 7.  DC AP   CS tau=0.007 % 19. HMM AP  CS
% 8.  DC TP1  PS tau=0.007 % 20. HMM TP1 PS
% 9.  DC TP1  BS tau=0.007 % 21. HMM TP1 BS
% 10. DC TP1  CS tau=0.007 % 22. HMM TP1 CS 
% 11. DC TP2  PS tau=0.007 % 23. HMM TP2 PS 
% 12. DC TP2  BS tau=0.007 % 24. HMM TP2 BS 
% 13. DC TP2  CS tau=0.007 % 25. HMM TP2 CS 

% simulation iterations
% -------------------------------------------------------------------------
for s = 1:n_s
    clear x X
    catch_i = [];

    % determine catch trials
    for run = 1:5
        load(sprintf('%s/sub-%02d_ses-1_ph_run-%d.mat', seqdir, SJs(s), run))

        % collect catch trials
        catch_i = [catch_i (find(C(:,2) == 2)' + 800*(run-1))];
    end
    inds = 1:4000;
    inds(catch_i) = [];

    % Models #2-4
    src_dir     = sprintf('%s/CD_logtaus/sub-%02d/sub-%02d_tau_0.0072208090_SP_CD.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,1:3)           = zscore([reg.predictive_surprise(inds)', ...            % regressors of interest
                      reg.bayesian_surprise(inds)', ...
                      reg.confidence_corrected_surprise(inds)'])        ; % z-scored regressors

    % Models #5-7
    src_dir     = sprintf('%s/CD_logtaus/sub-%02d/sub-%02d_tau_0.0072208090_AP_CD.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,4:6)    = zscore([reg.predictive_surprise(inds)', ...             % regressors of interest
                      reg.bayesian_surprise(inds)', ...
                      reg.confidence_corrected_surprise(inds)'])        ; % z-scored regressors

    % Models #8-10
    src_dir     = sprintf('%s/CD_logtaus/sub-%02d/sub-%02d_tau_0.0072208090_TP_CD.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,7:9)    = zscore([reg.predictive_surprise(inds)', ...             % regressors of interest
                      reg.bayesian_surprise(inds)', ...
                      reg.confidence_corrected_surprise(inds)'])        ; % z-scored regressors

    % Models #11-13
    src_dir     = sprintf('%s/CD_logtaus/sub-%02d/sub-%02d_tau_0.0072208090_TP_order2_CD.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,10:12)  = zscore([reg.predictive_surprise(inds)', ...             % regressors of interest
                      reg.bayesian_surprise(inds)', ...
                      reg.confidence_corrected_surprise(inds)'])        ; % z-scored regressors

    % Models #14-16
    src_dir     = sprintf('%s/HMM/paper/sub-%02d/sub-%02d_states-2_SP_nocatch_order1_paper_HMM.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,13:15)  = zscore([reg.predictive_surprise', ...                    % regressors of interest
                      reg.bayesian_surprise', ...
                      reg.confidence_corrected_surprise'])              ; % z-scored regressors

    % Models #17-19
    src_dir     = sprintf('%s/HMM/paper/sub-%02d/sub-%02d_states-2_AP_nocatch_order1_paper_HMM.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,16:18)  = zscore([reg.predictive_surprise', ...                     % regressors of interest
                      reg.bayesian_surprise', ...
                      reg.confidence_corrected_surprise'])              ; % z-scored regressors

    % Models #20-22
    src_dir     = sprintf('%s/HMM/paper/sub-%02d/sub-%02d_states-2_TP_nocatch_order1_paper_HMM.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,19:21)  = zscore([reg.predictive_surprise', ...                     % regressors of interest
                      reg.bayesian_surprise', ...
                      reg.confidence_corrected_surprise'])              ; % z-scored regressors

    % Models #23-25
    src_dir     = sprintf('%s/HMM/paper/sub-%02d/sub-%02d_states-2_TP_nocatch_order2_paper_HMM.mat', regdir, SJs(s), SJs(s)); % data source directory
    reg         = load(src_dir)                                         ; % model regressor structure                                                                
    x(:,22:24)  = zscore([reg.predictive_surprise', ...                     % regressors of interest
                      reg.bayesian_surprise', ...
                      reg.confidence_corrected_surprise'])              ; % z-scored regressors

    n           = size(x,1)                                             ; % number of data points
    X{1}        = ones(n,1)                                             ; % null model design matrix
    for k = 1:n_m-1
        X{k+1}    = [ones(n,1) x(:,k)]                                  ; % model-of-interest design matrix
    end

    % generative model iterations
    % ---------------------------------------------------------------------
    for m = 1:n_m
        for l = 1:n_l
            lambda_tbu = lambdas(l);

            % generative model sampling
            % -----------------------------------------------------------------
            y = mvnrnd(X{m}*beta_tbu(1:size(X{m},2)),1/(lambda_tbu)*eye(n))';

            % analysis model iterations
            % -----------------------------------------------------------------
            for i =  1:n_m

                % user information
                % fprintf('Generative model %d of %d, data sample %d of %d, analysis model %d of %d \n', m, n_m, s,n_s, i,n_m)

                % spm_vi_glm estimation
                % -------------------------------------------------------------
                glm             = []                                            ; % structure initialization
                glm.X           = X{i}                                          ; % design matrix
                glm.y           = zscore(y)                                     ; % data
                glm.n           = size(glm.X,1)                                 ; % number of data points
                glm.p           = size(glm.X,2)                                 ; % number of analysis model regression parameters
                glm.mu_beta     = zeros(glm.p,1)                                ; % expectation parameter of p(\beta) = N(\beta;0_p,\alpha^{-1}I_p)
                glm.alpha       = alpha                                         ; % precision parameter of p(\beta) = N(\beta; 0, \alpha^{-1}I_p)
                glm.beta_lambda = beta_lambda                                   ; % shape  parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                glm.gama_lambda = gama_lambda                                   ; % scalar parameter of p(\lambda)  = G(\lambda,\beta_\lambda, \gamma_\lambda) 
                glm.p           = size(glm.X,2)                                 ; % number of beta parameters
                glm.n_i         = n_i                                           ; % maximum number of iterations
                glm.delta       = delta                                         ; % variational free energy convergence criterion
                glm             = spm_vi_glm(glm)                               ; % estimation

                % recording
                lme(s,m,i,l)    = glm.F_max;
                lambda(s,m,i,l) = glm.bc_lambda_i(glm.i_c);
                if i > 1
                    beta(s,m,i,l) = glm.m_beta_i(2,glm.i_c);
                end

            end
        end
    end
end

fn = sprintf('%sSBL_modelrecovery_ffx_HPC_%03d.mat', resdir, str2num(val)); 
save(fn, 'lme', 'beta', 'lambda', '-v7.3')
end


% Auxillary functions, including spm_vi_glm and other support from SPM12.
% This may be useful for compiling this code and running it on a computing
% cluster.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function glm = spm_vi_glm(glm)

% This function implements a free-form variational inference approach for 
% general linear models of the form
%
%    y = X\beta + \varepsilon, \varepsilon ~ N(0_n,\lambda^{-1}I_n)
%
% where y,\varepsilon \in \mathbb{R}^n, X \in \mathbb{R}^{n \times p}, 
% \beta in \mathbb{R}^p, and \lambda > 0.
%   
%   Inputs
%       glm  : glm structure with required fields
%           .y              : n x 1 data vector
%           .X              : n x p design matrix
%           .n              : number of data points
%           .p              : number of regression parameters
%           .alpha          : precision parameter of p(\beta)
%           .beta_lambda    : shape parameter of p(\lambda)  
%           .gama_lambda    : scale parameter of p(\lambda)
%           .n_i            : maximum number of iterations
%           .delta          : variational free energy convergence criterion
%
%   Outputs
%       glm  : input structure with additional fields
%           .m_beta_i       : p x n_i     q(\beta) expectation iterands
%           .S_beta_i       : p x p x n_i q(\beta) covariance iterands
%           .b_lambda_i     : 1 x n_i     q(\lambda) shape iterands         
%           .c_lambda_i     : 1 x n_i     q(\lambda) scale iterands 
%           .bc_lambda_i    : 1 x n_i     \lambda expected value iterands
%           .F_i            : 1 x n_i     variational free energy iterands
%           .i_c            : iteration of convergence
%           .F_max          : maximized variational free energy
%
%   Author - Dirk Ostwald
% -------------------------------------------------------------------------

% precomputation of repeatedly used quantities
% -------------------------------------------------------------------------
glm.XTX                 = glm.X'*glm.X                                      ; % X'*X
glm.XTY                 = glm.X'*glm.y                                      ; % X'*y

% iterands array initialization
% -------------------------------------------------------------------------
glm.m_beta_i            = NaN(glm.p,glm.n_i+1)                              ; % variational expectation iterands array
glm.S_beta_i            = NaN(glm.p,glm.p,glm.n_i+1)                        ; % variational covariance iterands array
glm.b_lambda_i          = NaN(glm.p,glm.n_i+1)                              ; % variational shape parameter iterands array
glm.c_lambda_i          = NaN(glm.p,glm.p,glm.n_i+1)                        ; % variational scale parameter iterands array
glm.bc_lambda_i         = NaN(1,glm.n_i+1)                                  ; % lambda variational expectation iterands array
glm.F_i                 = NaN(1,glm.n_i+1)                                  ; % variational free energy iterands array

% algorithm initialization
% -------------------------------------------------------------------------
glm.m_beta              = glm.mu_beta                                       ; % m_\beta^{(0)}
glm.S_beta              = glm.alpha^(-1)*eye(glm.p)                         ; % S_\beta^{(0)}
glm.b_lambda            = glm.beta_lambda                                   ; % b_lambda^{(0)}
glm.c_lambda            = (glm.n)/2 + glm.gama_lambda                       ; % c_lambda^{(0)}
glm.bc_lambda           = glm.b_lambda.*glm.c_lambda                        ; % <\lambda>_{q^{(0)}(\lambda)}(\lambda)
glm.F                   = 0                                                 ; % F^{(0)}

% initialized iterands recording
glm.m_beta_i(:,1)       = glm.m_beta;
glm.S_beta_i(:,:,1)     = glm.S_beta;
glm.b_lambda_i(1)       = glm.b_lambda;
glm.c_lambda_i(1)       = glm.c_lambda;
glm.bc_lambda_i(1)      = glm.bc_lambda;
glm.F_i(1)              = glm.F;

% algorithm iterations
for i = 2:glm.n_i 
    
    % q^{(i)}(\beta) evaluation
    glm                 = spm_vi_beta(glm)                                  ; % regression parameter variational distribution update
    glm.m_beta_i(:,i)   = glm.m_beta                                        ; % variational expectation parameter recording
    glm.S_beta_i(:,:,i) = glm.S_beta                                        ; % variational covariance parameter recording
     
    % q^{(i)}(\lambda) evaluation
    glm                 = spm_vi_lambda(glm)                                ; % observation noise parameter variational distribution update
    glm.b_lambda_i(i)   = glm.b_lambda                                      ; % variational shape parameter recording
    glm.c_lambda_i(i)   = glm.c_lambda                                      ; % variational scale parameter recording
    glm.bc_lambda_i(i)  = glm.bc_lambda                                            ; % \lambda expectation recording
         
    % F^{(i)} evaluation 
    glm                 = spm_vi_F(glm)                                     ; % variational free energy update
    glm.F_i(i)          = glm.F                                             ; % variational free energy recording
    
    % convergence assessment
    delta_F             = glm.F_i(i) - glm.F_i(i-1);
      
    if i > 2
        if delta_F < 0
            fprintf('* Warning: decrease in F of %1.4f per cent* \n',100*(delta_F/glm.F));
            glm.i_c = i                                                     ; % iteration of convergence
            break;
        elseif abs(delta_F/glm.F) < glm.delta
            glm.i_c = i                                                     ; % iteration of convergence
            break;
        end
    end  
    glm.i_c = i                                                             ; % iteration of convergence

end

% converged variational free energy
glm.F_max = glm.F_i(glm.i_c);

% remove dynamic structure fields
% -------------------------------------------------------------------------
rf  = {'XTX','XTY','m_beta','S_beta','b_lambda','c_lambda','bc_lambda','F'} ; % fields to be removed
glm = rmfield(glm,rf)                                                       ; % field removal

end

% free-form variational inference algorithm subfunctions
% -------------------------------------------------------------------------
function glm = spm_vi_beta(glm)

% This function evaluates the variational parameter update for
%
%                q^(\beta) = N(\beta;m_beta,S_beta)
%           
%   Inputs
%           glm : glm structure with required fields
%               .p          : number of beta parameters
%               .XTX        : precomputed design matrix cross-product
%               .XTY        : percomputed design matrix-data product
%               .bc_lambda  : scalar expected value of \lambda
%               .alpha      : scalar precision parameter of p(\beta)
%   
%   Outputs
%           glm     : glm structure with updated fields
%               .m_beta     : p x 1 q(\beta) expectation parameter
%               .S_beta     : p x p q(\beta) covariance parameter
%
%   Author - Dirk Ostwald
% -------------------------------------------------------------------------
glm.S_beta = inv(glm.bc_lambda*glm.XTX + glm.alpha*eye(glm.p))              ; % S_\beta^{(i)}
glm.m_beta = glm.S_beta*glm.bc_lambda*glm.XTY                               ; % m_\beta^{(i)}
end

function glm = spm_vi_lambda(glm)

% This function evaluates the variational parameter update for
%
%                q^(\lambda) = G(\lambda;b_lambda,c_lambda)
%           
%   Inputs
%           glm : glm structure with required fields
%               .y          : n x 1 data vector
%               .X          : n x p design matrix
%               .m_beta     : p x 1 q(\beta) expectation parameter
%               .S_Beta     : p x p q(\beta) covariance parameter
%               .XTX        : precomputed design matrix cross-product

%   Outputs
%           glm     : glm structure with updated fields
%               .b_lambda   : scalar q(\lambda) parameter
%               .bc         : scalar q(\lambda) expected value
%
%   Author - Dirk Ostwald
% -------------------------------------------------------------------------
e                       = glm.y - glm.X*glm.m_beta                          ; % variational expectation-based data prediction error
G                       = trace(glm.S_beta*glm.XTX) + e'*e                  ; % b_lambda update term
glm.b_lambda            = 1./(G./2 + 1./glm.beta_lambda)                    ; % b_\lambda^{(i)}
glm.bc_lambda           = glm.b_lambda*glm.c_lambda                         ; % <\lambda>_{q^{(i)}(\lambda)}(\lambda)
end

function glm = spm_vi_F(glm)

% This function evaluates the variational free energy (evidence lower bound)
% for the free-form mean-field variational Bayes algorithm.
%
%   Inputs
%           glm     : glm structure with required fields
%               .y          : n x 1 data vector
%               .X          : n x p design matrix
%               .p          : number of beta parameters
%               .n          : number of data points
%               .XTX        : precomputed design matrix cross-product
%               .XTY        : percomputed design matrix-data product
%               .m_beta     : p x 1 q(\beta) expectation parameter
%               .S_Beta     : p x p q(\beta) covariance parameter
%               .b_lambda   : scalar q(\lambda) shape parameter          
%               .c_lambda   : scalar q(\lambda) scale parameter          
%               .bc         : scalar expected value of \lambda
%               .alpha      : scalar precision parameter of p(\beta)
%
%   Outputs
%           glm     : glm structure with additional field
%               .F          : scalar variational free energy
%
%   Author - Dirk Ostwald
% -------------------------------------------------------------------------

% average likelihood evaluation
% -------------------------------------------------------------------------
T1  = -(1/2)*glm.n*log(2*pi)                                                ; % average log likelihood term 1 
T2  = -(1/2)*glm.bc_lambda*(glm.y-glm.X*glm.m_beta)'*(glm.y-glm.X*glm.m_beta); % average log likelihood term 2 
T3  = -(1/2)*glm.bc_lambda*trace(glm.S_beta*glm.XTX)                        ; % average log likelihood term 3 
T4  = (1/2)*glm.n*(psi(glm.c_lambda) + log(glm.b_lambda))                   ; % average log likelihood term 4  
Lav = T1 + T2 + T3 + T4                                                     ; % average log likelihood 

% KL divergence terms
% -------------------------------------------------------------------------
KL_beta   = spm_kl_normal(glm.m_beta,glm.S_beta,glm.mu_beta,glm.alpha^(-1)*eye(glm.p)); % KL(q^{(i)}(\beta)||p(\beta))
KL_lambda = spm_kl_gamma(glm.b_lambda,glm.c_lambda,glm.beta_lambda,glm.gama_lambda); % KL(q^{(i)}(\lambda)||p(\lambda))

% variational free energy evaluation 
% -------------------------------------------------------------------------
glm.F               = Lav - (KL_beta + KL_lambda);
end

% KL divergence subfunctions
% -------------------------------------------------------------------------
function [d] = spm_kl_gamma (b_q,c_q,b_p,c_p)

% KL divergence between two Gamma densities
% FORMAT [d] = spm_kl_gamma (b_q,c_q,b_p,c_p)
%
% KL (Q||P) = <log Q/P> where avg is wrt Q
%
% b_q, c_q    Parameters of first Gamma density
% b_p, c_p    Parameters of second Gamma density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kl_gamma.m 2696 2009-02-05 20:29:48Z guillaume $
digamma_c_q = psi(c_q);
d           = (c_q-1)*digamma_c_q-log(b_q)-c_q-gammaln(c_q);
d           = d+gammaln(c_p)+c_p*log(b_p)-(c_p-1)*(digamma_c_q+log(b_q));
d           = d+b_q*c_q/b_p;

end

function [d] = spm_kl_normal(m_q,c_q,m_p,c_p)

% KL divergence between two multivariate normal densities
% FORMAT [d] = spm_kl_normal (m_q,c_q,m_p,c_p)
%
% KL (Q||P) = <log Q/P> where avg is wrt Q
%
% between two Normal densities Q and P
%
% m_q, c_q    Mean and covariance of first Normal density
% m_p, c_p    Mean and covariance of second Normal density
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kl_normal.m 2696 2009-02-05 20:29:48Z guillaume $

d           = length(m_q);
m_q         = m_q(:);
m_p         = m_p(:);
Term1       = 0.5*spm_logdet(c_p)-0.5*spm_logdet(c_q);
inv_c_p     = inv(c_p);
Term2       = 0.5*trace(inv_c_p*c_q)+0.5*(m_q-m_p)'*inv_c_p*(m_q-m_p);
d           = Term1 + Term2 - 0.5*d;
end

function H = spm_logdet(C)

% Compute the log of the determinant of positive (semi-)definite matrix C
% FORMAT H = spm_logdet(C)
% H = log(det(C))
%
% spm_logdet is a computationally efficient operator that can deal with
% full or sparse matrices. For non-positive definite cases, the determinant
% is considered to be the product of the positive singular values.
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston and Ged Ridgway
% $Id: spm_logdet.m 6321 2015-01-28 14:40:44Z karl $

% Note that whether sparse or full, rank deficient cases are handled in the
% same way as in spm_logdet revision 4068, using svd on a full version of C


% remove null variances
%--------------------------------------------------------------------------
i       = find(diag(C));
C       = C(i,i);
[i,j,s] = find(C);
if any(isnan(s)), H = nan; return; end

% TOL = max(size(C)) * eps(max(s)); % as in MATLAB's rank function
%--------------------------------------------------------------------------
TOL   = 1e-16;

if any(i ~= j)
    
    % assymetric matrix
    %------------------------------------------------------------------
    if norm(spm_vec(C - C'),inf) > TOL
        
        s = svd(full(C));
        
    else
        
        % non-diagonal sparse matrix
        %------------------------------------------------------------------
        if issparse(C)
            
            % Note p is unused but requesting it can make L sparser
            %--------------------------------------------------------------
            [L,nondef,p] = chol(C, 'lower', 'vector');
            if ~nondef
                
                % pos. def. with Cholesky decomp L, and det(C) = det(L)^2
                %----------------------------------------------------------
                H = 2*sum(log(full(diag(L))));
                return
                
            end
            s = svd(full(C));
            
        else
            
            % non-diagonal full matrix
            %--------------------------------------------------------------
            try
                R = chol(C);
                H = 2*sum(log(diag(R)));
                return
            catch
                s = svd(C);
            end
        end
    end
end

% if still here, singular values in s (diagonal values as a special case)
%--------------------------------------------------------------------------
H     = sum(log(s(s > TOL & s < 1/TOL)));
end

function [vX] = spm_vec(X,varargin)

% Vectorise a numeric, cell or structure array - a compiled routine
% FORMAT [vX] = spm_vec(X)
% X  - numeric, cell or stucture array[s]
% vX - vec(X)
%
% See spm_unvec
%__________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_vec.m 6110 2014-07-21 09:36:13Z karl $


%error('spm_vec.c not compiled - see Makefile')

% initialise X and vX
%--------------------------------------------------------------------------
if nargin > 1
    X = [{X},varargin];
end


% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    vX = X(:);

% vectorise logical arrays
%--------------------------------------------------------------------------
elseif islogical(X)
    vX = X(:);

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
elseif isstruct(X)
    vX = [];
    f   = fieldnames(X);
    X    = X(:);
    for i = 1:numel(f)
        vX = cat(1,vX,spm_vec({X.(f{i})}));
    end

% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
elseif iscell(X)
    vX   = [];
    for i = 1:numel(X)
        vX = cat(1,vX,spm_vec(X{i}));
    end
else
    vX = [];
end
end
