function SBL_modelrecovery_rfx(val)
% This function inputs the results from the FFX model recovery and 
% performs RFX model recovery for the four different levels of analyses:
% 1) Null Model vs DC vs HMM
% 2) DC TP1 vs TP2
% 3) DC SP vs AP vs TP1
% 4) DC TP1 PS vs BS vs CS
%
% Authors - Miro Grundei, Sam Gijsen
% --------------------------------------------------------------

% Load in FFX results for all iterations

n_m     = 25; % number of models
n_l     = 6;  % number of noise levels 

resdir      = fullfile('/scratch/samgijsen/results/rfx/model_recovery/')                       ; % results directory

% load in the FFX model recovery results
load(sprintf('/scratch/samgijsen/results/F/model_recovery/SBL_modelrecovery_ffx_HPC_%03d.mat', str2num(val)))

% penalize DC models
for l = 1:n_l
    lme_max(:,:,[2:13],l) = lme_max(:,:,[2:13],l) - mean_delta_F(l);
end
lme = lme_max;

clear lme_max mean_delta_F


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

% Perform RFX for Null vs DC vs HMM with following TBU models
% > ------------------------------------------------------------- <
% NM #1
% CD  SP  PS #2
% CD  AP  BS #6
% CD  TP1 CS #10
% CD  TP2 BS #12
% HMM SP  PS #14
% HMM AP  BS #18
% HMM TP1 CS #22
% HMM TP2 BS #24

% number of families and models to recover
n_f             = 3; 
tbu_data        = [1, 2, 6, 10, 12, 14, 18, 22, 24];

% set parameters
family.infer    = 'RFX';
family.partition = [1, repmat(2,1,12), repmat(3,1,12)];
family.names    = {'NULL', 'CD', 'HMM'};
family.Nsamp    = 5e5;

% initialize
l1_fam_xp       = NaN(numel(tbu_data), n_f, n_l);
l1_fam_expr     = NaN(numel(tbu_data), n_f, n_l);
l1_mod_xp       = NaN(numel(tbu_data), n_m, n_l);
l1_mod_expr     = NaN(numel(tbu_data), n_m, n_l);

% loop over TBU models and noise levels
for i = 1:numel(tbu_data)
    for l = 1:n_l
       [f, m]              = spm_compare_families(squeeze(lme(:,tbu_data(i),:,l)), family);
       l1_fam_xp(i,:,l)    = f.xp;  
       l1_fam_expr(i,:,l)  = f.exp_r;
       l1_mod_xp(i,:,l)    = m.xp;
       l1_mod_expr(i,:,l)  = m.exp_r;
    end
end

% save results
fn = sprintf('%sSBL_modelrecovery_RFX_HPC_l1_%03d.mat', resdir, str2num(val)); 
save(fn, 'l1_fam_xp', 'l1_fam_expr', 'l1_mod_xp', 'l1_mod_expr', '-v7.3')

clear l1_fam_xp l1_mod_xp l1_fam_expr l1_mod_expr

% TP1 vs TP2
% > ------------------------------------------------------------- <
% DC TP1 BS #9
% DC TP1 CS #10
% DC TP2 CS #13
% DC TP2 BS #12
% DC TP1 PS #8
% DC TP2 PS #11

% number of families and models to recover
n_f             = 2; 
tbu_data        = [9, 10, 13, 12, 8, 11];
relevant        = [8:13]; % only consider these models

% set parameters
family.infer    = 'RFX';
family.partition = [1 1 1 2 2 2];
family.names    = {'TP1', 'TP2'};
family.Nsamp    = 5e5;

% initialize
l2_fam_xp       = NaN(numel(tbu_data), n_f, n_l);
l2_fam_expr     = NaN(numel(tbu_data), n_f, n_l);
l2_mod_xp       = NaN(numel(tbu_data), numel(relevant), n_l);
l2_mod_expr     = NaN(numel(tbu_data), numel(relevant), n_l);

% loop over TBU models and noise levels
for i = 1:numel(tbu_data)
    for l = 1:n_l
       [f, m]              = spm_compare_families(squeeze(lme(:,tbu_data(i),relevant,l)), family);
       l2_fam_xp(i,:,l)    = f.xp;  
       l2_fam_expr(i,:,l)  = f.exp_r;
       l2_mod_xp(i,:,l)    = m.xp;
       l2_mod_expr(i,:,l)  = m.exp_r;
    end
end

% save results
fn = sprintf('%sSBL_modelrecovery_RFX_HPC_l2_%03d.mat', resdir, str2num(val)); 
save(fn, 'l2_fam_xp', 'l2_fam_expr', 'l2_mod_xp', 'l2_mod_expr', '-v7.3')

clear l2_fam_xp l2_mod_xp l2_fam_expr l2_mod_expr

% SP AP TP1
% > ------------------------------------------------------------- <
% DC SP  BS #3
% DC AP  BS #6
% DC TP1 BS #9
% DC SP  CS #4
% DC AP  CS #7
% DC TP1 CS #10
% DC SP  PS #2
% DC AP  PS #5
% DC TP1 PS #8

% number of families and models to recover
n_f             = 3; 
tbu_data        = [3, 6, 9, 4, 7, 10, 2, 5, 8];
relevant        = [2:10]; % only consider these models

% set parameters
family.infer    = 'RFX';
family.partition = [1 1 1 2 2 2 3 3 3];
family.names    = {'SP', 'AP', 'TP1'};
family.Nsamp    = 5e5;

% initialize
l3_fam_xp       = NaN(numel(tbu_data), n_f, n_l);
l3_fam_expr     = NaN(numel(tbu_data), n_f, n_l);
l3_mod_xp       = NaN(numel(tbu_data), numel(relevant), n_l);
l3_mod_expr     = NaN(numel(tbu_data), numel(relevant), n_l);

% loop over TBU models and noise levels
for i = 1:numel(tbu_data)
    for l = 1:n_l
       [f, m]              = spm_compare_families(squeeze(lme(:,tbu_data(i),relevant,l)), family);
       l3_fam_xp(i,:,l)    = f.xp;  
       l3_fam_expr(i,:,l)  = f.exp_r;
       l3_mod_xp(i,:,l)    = m.xp;
       l3_mod_expr(i,:,l)  = m.exp_r;
    end
end

% save results
fn = sprintf('%sSBL_modelrecovery_RFX_HPC_l3_%03d.mat', resdir, str2num(val)); 
save(fn, 'l3_fam_xp', 'l3_fam_expr', 'l3_mod_xp', 'l3_mod_expr', '-v7.3')

clear l3_fam_xp l3_mod_xp l3_fam_expr l3_mod_expr


% PS BS CS - within TP1
% > ------------------------------------------------------------- <
% DC TP1 PS #8
% DC TP1 BS #9
% DC TP1 CS #10

% models to recover 
tbu_data        = [8:10];
relevant        = [8:10]; % only consider these models

Nsamp    = 1e6;
do_plot  = 0;
sampling = 1;
ecp      = 1;

% only 3 models for non-family analysis
n_m     = 3;

% initialize
l4_mod_xp       = NaN(numel(tbu_data), numel(relevant), n_l);
l4_mod_expr     = NaN(numel(tbu_data), numel(relevant), n_l);
l4_mod_pxp      = NaN(numel(tbu_data), numel(relevant), n_l);

% loop over TBU models and noise levels
for i = 1:numel(tbu_data)
    for l = 1:n_l
       [~,exp_r,xp,pxp,~]   = spm_BMS(squeeze(lme(:,tbu_data(i),relevant,l)), Nsamp, do_plot, sampling, ecp);  
       l4_mod_xp(i,:,l)     = xp;
       l4_mod_expr(i,:,l)   = exp_r;
       l4_mod_pxp(i,:,l)   = pxp;
    end
end

% save results
fn = sprintf('%sSBL_modelrecovery_RFX_HPC_l4_%03d.mat', resdir, str2num(val)); 
save(fn, 'l4_mod_xp', 'l4_mod_expr', 'l4_mod_pxp', '-v7.3')
end

%% Additional SPM Support Functions - Compile Requirements
% This may be useful for compiling this code and running it on a computing
% cluster.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [family,model] = spm_compare_families (lme,family)
% Bayesian comparison of model families for group studies 
% FORMAT [family,model] = spm_compare_families (lme,family)
%
% INPUT:
%
% lme           - array of log model evidences 
%                   rows: subjects
%                   columns: models (1..N)
%
% family        - data structure containing family definition and inference parameters:
%                  .infer='RFX' or 'FFX' (default)
%                  .partition  [1 x N] vector such that partition(m)=k signifies that
%                              model m belongs to family k (out of K) eg. [1 1 2 2 2 3 3]
%                  .names      cell array of K family names eg, {'fam1','fam2','fam3'}
%                  .Nsamp      RFX only: Number of samples to get (default=1e4)
%                  .prior      RFX only: 'F-unity' alpha0=1 for each family (default)
%                              or 'M-unity' alpha0=1 for each model (not advised)
%
% OUTPUT:
%
% family        - RFX only:  
%                   .alpha0       prior counts 
%                   .exp_r        expected value of r
%                   .s_samp       samples from posterior
%                   .xp           exceedance probs
%                - FFX only: 
%                   .prior        family priors
%                   .post         family posteriors
%
% model          - RFX only: 
%                   .alpha0        prior counts
%                   .exp_r         expected value of r
%                   .r_samp        samples from posterior
%                - FFX only: 
%                   .subj_lme      log model ev without subject effects
%                   .prior         model priors
%                   .like          model likelihoods 
%                                  (likelihood scaled to unity for most
%                                  likely model)
%                   .posts         model posteriors
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_compare_families.m 6052 2014-06-17 09:38:13Z will $

try
    infer=family.infer;
catch
    disp('Error in spm_compare_families: inference method not specified');
    return
end

try
    partition=family.partition;
catch
    disp('Error in spm_compare_families: partition not specified');
    return
end

try
    names=family.names;
catch
    disp('Error in spm_compare_families: names not specified');
    return
end

if strcmp(infer,'RFX')
    try
        Nsamp=family.Nsamp;
    catch
        Nsamp=1e4;
        family.Nsamp=Nsamp;
    end
    
    try
        prior=family.prior;
    catch
        prior='F-unity';
        family.prior='F-unity';
    end
end

% Number of models
N=length(partition);

% Number of families in partition
K=length(unique(partition));

% Size of families 
for i=1:K,
    ind{i}=find(partition==i);
    fam_size(i)=length(ind{i});
end

if strcmp(infer,'FFX')
    
    family.prior = [];
    % Family priors
    for i=1:K,
        family.prior(i)=1/K;
    end
    
    % Model priors
    for i=1:N,
        model.prior(i)=1/fam_size(partition(i));
    end
    
    slme=sum(lme,1);
    slme=slme-(max(slme,[],2))*ones(1,N);
    
    % Model likelihoods
    model.subj_lme=slme;
    model.like=exp(slme);
    
    % Model posterior
    num=model.prior.*model.like;
    model.post=num/sum(num);
    
    % Family posterior
    for i=1:K,
        family.like(i)=sum(model.like(ind{i}));
        family.post(i)=sum(model.post(ind{i}));
    end
    return;
end
    
% Set model priors 
switch prior,
    case 'F-unity',
        for i=1:K,
            model.alpha0(ind{i})=1/fam_size(i);
        end
        family.alpha0=ones(1,K);
    case 'M-unity',
        model.alpha0=ones(1,N);
        for i=1:K,
            family.alpha0(i)=fam_size(i);
        end
    otherwise
        disp('Error in spm_compare_families:Unknown prior');
end

% Get model posterior
[exp_r,xp,r_samp,g_post]=spm_BMS_gibbs(lme,model.alpha0,Nsamp);
model.exp_r=exp_r;
model.xp=xp;
model.r_samp=r_samp;
model.g_post=g_post;

% Get stats from family posterior
for i=1:K,
    ri=r_samp(:,ind{i});
    family.s_samp(:,i)=sum(ri,2);
    family.exp_r(i)=mean(family.s_samp(:,i));
end

% Family exceedence probs
xp = zeros(1,K);
r=family.s_samp;
[y,j]=max(r,[],2);
tmp=histc(j,1:K)';
family.xp=tmp/Nsamp;
end


function [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp)
% Bayesian model selection for group studies using Gibbs sampling
% FORMAT [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp)
%
% INPUT:
% lme      - array of log model evidences 
%              rows: subjects
%              columns: models (1..Nk)
% alpha0   - [1 x Nk] vector of prior model counts
% Nsamp    - number of samples (default: 1e6)
% 
% OUTPUT:
% exp_r    - [1 x  Nk] expectation of the posterior p(r|y)
% xp       - exceedance probabilities
% r_samp   - [Nsamp x Nk] matrix of samples from posterior
% g_post   - [Ni x Nk] matrix of posterior probabilities with 
%            g_post(i,k) being post prob that subj i used model k
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_gibbs.m 6381 2015-03-17 17:55:09Z will $


if nargin < 3 || isempty(Nsamp)
    Nsamp = 1e4;
end

Ni = size(lme,1);  % number of subjects
Nk = size(lme,2);  % number of models

% prior observations
%--------------------------------------------------------------------------
if nargin < 2 || isempty(alpha0)
    alpha0 = ones(1,Nk);    
end
alpha0     = alpha0(:)';

% Initialise; sample r from prior
r  = zeros(1,Nk);
for k = 1:Nk
    r(:,k) = gamrnd(alpha0(k),1);
end
sr = sum(r,2);
for k = 1:Nk
    r(:,k) = r(:,k)./sr;
end
        
% Subtract max evidence for subject 
lme = lme - max(lme,[],2)*ones(1,Nk);

% Gibbs sampling 
r_samp = zeros(Nsamp,Nk);
g_post = zeros(Ni,Nk);

for samp = 1:2*Nsamp
    
    mod_vec = sparse(Ni,Nk);
    % Sample m's given y, r
    for i = 1:Ni
        % Pick a model for this subject
        u         = exp(lme(i,:) + log(r)) + eps;
        g         = u / sum(u);
        gmat(i,:) = g;
        modnum    = spm_multrnd(g,1);
        mod_vec(i,modnum) = 1;
    end
    
    % Sample r's given y, m
    beta          = sum(mod_vec,1);
    alpha         = alpha0+beta;
    for k = 1:Nk
        r(:,k)    = gamrnd(alpha(k),1);
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k)    = r(:,k) ./ sr;
    end

    % Only keep last Nsamp samples
    if samp > Nsamp
        r_samp(samp-Nsamp,:) = r;
        g_post    = g_post+gmat;
    end
    
    %if mod(samp,1e4)==0
        %fprintf('%d samples out of %d\n',samp,2*Nsamp);
    %end
    
end
g_post = g_post/Nsamp;

% Posterior mean
exp_r = mean(r_samp,1);

% Exceedence probs
xp    = zeros(1,Nk);
[y,j] = max(r_samp,[],2);
tmp   = histc(j,1:Nk)';
xp    = tmp / Nsamp;
end


function [m] = spm_multrnd(p,N)
% Sample from multinomial distribution
% FORMAT [m] = spm_multrnd(p,N)
%
% p    - [M x 1] vector of probabilities
% N    - Number of samples to generate
% 
% m    - [N x 1] vector of samples, where each sample is number from 1 to M
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_multrnd.m 3190 2009-06-08 17:13:36Z guillaume $

cp = [0; cumsum(p(:))];
m  = zeros(N,1);
for n=1:N
    m(n) = find(rand > cp, 1, 'last');
end
end

function r = spm_gamrnd(a,b,varargin)
% Random arrays from gamma distribution - a compiled routine
% FORMAT r = spm_gamrnd(a,b,m,n,...)
%
% a        - shape parameter
% b        - scale parameter
% m,n,...  - dimensions of the output array [optional]
%
% r        - array of random numbers chosen from the gamma distribution
%__________________________________________________________________________
%
% Reference
% 
% George Marsaglia and Wai Wan Tsang, "A Simple Method for Generating Gamma
% Variables": ACM Transactions on Mathematical Software, Vol. 26, No. 3,
% September 2000, Pages 363-372
% http://portal.acm.org/citation.cfm?id=358414
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_gamrnd.m 3251 2009-07-06 17:29:44Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_gamrnd.c not compiled - see Makefile');
end


function [alpha,exp_r,xp,pxp,bor] = spm_BMS (lme, Nsamp, do_plot, sampling, ecp, alpha0)
% Bayesian model selection for group studies
% FORMAT [alpha,exp_r,xp,pxp,bor] = spm_BMS (lme, Nsamp, do_plot, sampling, ecp, alpha0)
% 
% INPUT:
% lme      - array of log model evidences 
%              rows: subjects
%              columns: models (1..Nk)
% Nsamp    - number of samples used to compute exceedance probabilities
%            (default: 1e6)
% do_plot  - 1 to plot p(r|y)
% sampling - use sampling to compute exact alpha
% ecp      - 1 to compute exceedance probability
% alpha0   - [1 x Nk] vector of prior model counts
% 
% OUTPUT:
% alpha   - vector of model probabilities
% exp_r   - expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% pxp     - protected exceedance probabilities
% bor     - Bayes Omnibus Risk (probability that model frequencies 
%           are equal)
% 
% REFERENCES:
%
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009)
% Bayesian Model Selection for Group Studies. NeuroImage 46:1004-1017
%
% Rigoux, L, Stephan, KE, Friston, KJ and Daunizeau, J. (2014)
% Bayesian model selection for group studies�Revisited. 
% NeuroImage 84:971-85. doi: 10.1016/j.neuroimage.2013.08.065
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Klaas Enno Stephan, Will Penny
% $Id: spm_BMS.m 6442 2015-05-21 09:13:44Z will $

if nargin < 2 || isempty(Nsamp)
    Nsamp = 1e6;
end
if nargin < 3 || isempty(do_plot)
    do_plot = 0;
end
if nargin < 4 || isempty(sampling)
    sampling = 0;
end
if nargin < 5 || isempty(ecp)
    ecp = 1;
end

Ni      = size(lme,1);  % number of subjects
Nk      = size(lme,2);  % number of models
c       = 1;
cc      = 10e-4;

% prior observations
%--------------------------------------------------------------------------
if nargin < 6 || isempty(alpha0)
    alpha0  = ones(1,Nk);    
end
alpha   = alpha0;

% iterative VB estimation
%--------------------------------------------------------------------------
while c > cc,

    % compute posterior belief g(i,k)=q(m_i=k|y_i) that model k generated
    % the data for the i-th subject
    for i = 1:Ni,
        for k = 1:Nk,
            % integrate out prior probabilities of models (in log space)
            log_u(i,k) = lme(i,k) + psi(alpha(k))- psi(sum(alpha));
        end
        
        % exponentiate (to get back to non-log representation)
        u(i,:)  = exp(log_u(i,:)-max(log_u(i,:)));
        
        % normalisation: sum across all models for i-th subject
        u_i     = sum(u(i,:));
        g(i,:)  = u(i,:)/u_i;
    end
            
    % expected number of subjects whose data we believe to have been 
    % generated by model k
    for k = 1:Nk,
        beta(k) = sum(g(:,k));
    end

    % update alpha
    prev  = alpha;
    for k = 1:Nk,
        alpha(k) = alpha0(k) + beta(k);
    end
    
    % convergence?
    c = norm(alpha - prev);

end


% Compute expectation of the posterior p(r|y)
%--------------------------------------------------------------------------
exp_r = alpha./sum(alpha);


% Compute exceedance probabilities p(r_i>r_j)
%--------------------------------------------------------------------------
if ecp
    if Nk == 2
        % comparison of 2 models
        xp(1) = spm_Bcdf(0.5,alpha(2),alpha(1));
        xp(2) = spm_Bcdf(0.5,alpha(1),alpha(2));
    else
        % comparison of >2 models: use sampling approach
        xp = spm_dirichlet_exceedance(alpha,Nsamp);
    end
else
        xp = [];
end

posterior.a=alpha;
posterior.r=g';
priors.a=alpha0;
bor = spm_BMS_bor (lme',posterior,priors);

% Compute protected exceedance probs - Eq 7 in Rigoux et al.
pxp=(1-bor)*xp+bor/Nk;

% Graphics output (currently for 2 models only)
%--------------------------------------------------------------------------
if do_plot && Nk == 2
    % plot Dirichlet pdf
    %----------------------------------------------------------------------
    if alpha(1)<=alpha(2)
       alpha_now =sort(alpha,1,'descend');
       winner_inx=2;
    else
        alpha_now =alpha;
       winner_inx=1;
    end
    
    x1  = [0:0.0001:1];
    for i = 1:length(x1),
        p(i)   = spm_Dpdf([x1(i) 1-x1(i)],alpha_now);
    end
    fig1 = figure;
    axes1 = axes('Parent',fig1,'FontSize',14);
    plot(x1,p,'k','LineWidth',1);
    % cumulative probability: p(r1>r2)
    i  = find(x1 >= 0.5);
    hold on
    fill([x1(i) fliplr(x1(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
    v = axis;
    plot([0.5 0.5],[v(3) v(4)],'k--','LineWidth',1.5);
    xlim([0 1.05]);
    xlabel(sprintf('r_%d',winner_inx),'FontSize',18);
    ylabel(sprintf('p(r_%d|y)',winner_inx),'FontSize',18);
    title(sprintf('p(r_%d>%1.1f | y) = %1.3f',winner_inx,0.5,xp(winner_inx)),'FontSize',18);
    legend off
end


% Sampling approach ((currently implemented for 2 models only):
% plot F as a function of alpha_1
%--------------------------------------------------------------------------
if sampling
    if Nk == 2
        % Compute lower bound on F by sampling
        %------------------------------------------------------------------
        alpha_max = size(lme,1) + Nk*alpha0(1);
        dx        = 0.1;
        a         = [1:dx:alpha_max];
        Na        = length(a);
        for i=1:Na,
            alpha_s                = [a(i),alpha_max-a(i)];
            [F_samp(i),F_bound(i)] = spm_BMS_F(alpha_s,lme,alpha0);
        end
        if do_plot
        % graphical display
        %------------------------------------------------------------------
        fig2 = figure;
        axes2 = axes('Parent',fig2,'FontSize',14);
        plot(a,F_samp,'Parent',axes2,'LineStyle','-','DisplayName','Sampling Approach',...
            'Color',[0 0 0]);
        hold on;
        yy = ylim;
        plot([alpha(1),alpha(1)],[yy(1),yy(2)],'Parent',axes2,'LineStyle','--',...
            'DisplayName','Variational Bayes','Color',[0 0 0]);
        legend2 = legend(axes2,'show');
        set(legend2,'Position',[0.15 0.8 0.2 0.1],'FontSize',14);
        xlabel('\alpha_1','FontSize',18);
        ylabel('F','FontSize',18);
        end
    else
        fprintf('\n%s\n','Verification of alpha estimates by sampling not available.')
        fprintf('%s\n','This approach is currently only implemented for comparison of 2 models.');
    end
end
end

function F = spm_Bcdf(x,v,w)
% Inverse Cumulative Distribution Function (CDF) of Beta distribution
% FORMAT F = spm_Bcdf(x,v,w)
%
% x   - Beta variates (Beta has range [0,1])
% v   - Shape parameter (v>0)
% w   - Shape parameter (w>0)
% F   - CDF of Beta distribution with shape parameters [v,w] at points x
%__________________________________________________________________________
%
% spm_Bcdf implements the Cumulative Distribution Function for Beta
% distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Beta distribution has two shape parameters, v and w, and is
% defined for v>0 & w>0 and for x in [0,1] (See Evans et al., Ch5).
% The Cumulative Distribution Function (CDF) F(x) is the probability
% that a realisation of a Beta random variable X has value less than
% x. F(x)=Pr{X<x}: This function is usually known as the incomplete Beta
% function. See Abramowitz & Stegun, 26.5; Press et al., Sec6.4 for
% definitions of the incomplete beta function.
%
% Variate relationships:
%--------------------------------------------------------------------------
% Many: See Evans et al., Ch5
%
% Algorithm:
%--------------------------------------------------------------------------
% Using MATLAB's implementation of the incomplete beta finction (betainc).
%
% References:
%--------------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Bcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for x in [0,1] & strictly positive v & w.
% Return NaN if undefined.
md = ( x>=0  &  x<=1  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: F=1 when x=1
F(md & x==1) = 1;

%-Non-zero where defined & x>0, avoid special cases
Q  = find( md  &  x>0  &  x<1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
F(Q) = betainc(x(Qx),v(Qv),w(Qw));
end


function xp = spm_dirichlet_exceedance(alpha,Nsamp)
% Compute exceedance probabilities for a Dirichlet distribution
% FORMAT xp = spm_dirichlet_exceedance(alpha,Nsamp)
% 
% Input:
% alpha     - Dirichlet parameters
% Nsamp     - number of samples used to compute xp [default = 1e6]
% 
% Output:
% xp        - exceedance probability
%__________________________________________________________________________
%
% This function computes exceedance probabilities, i.e. for any given model
% k1, the probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% 
% Refs:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (in press)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $

if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);

% Perform sampling in blocks
%--------------------------------------------------------------------------
blk = ceil(Nsamp*Nk*8 / 2^28);
blk = floor(Nsamp/blk * ones(1,blk));
blk(end) = Nsamp - sum(blk(1:end-1));

xp = zeros(1,Nk);
for i=1:length(blk)
    
    % Sample from univariate gamma densities then normalise
    % (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
    % 209-230)
    %----------------------------------------------------------------------
    r = zeros(blk(i),Nk);
    for k = 1:Nk
        r(:,k) = gamrnd(alpha(k),1,blk(i),1);
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end
    
    % Exceedance probabilities:
    % For any given model k1, compute the probability that it is more
    % likely than any other model k2~=k1
    %----------------------------------------------------------------------
    [y, j] = max(r,[],2);
    xp = xp + histc(j, 1:Nk)';
    
end
xp = xp / Nsamp;
end


function [F_samp,F_bound] = spm_BMS_F (alpha,lme,alpha0)
% Compute two lower bounds on model evidence p(y|r) for group BMS
% 
% FORMAT [F_samp,F_bound] = spm_BMS_smpl_me (alpha,lme,alpha0)
% 
% INPUT:
% alpha     parameters of p(r|y)
% lme       array of log model evidences 
%              rows:    subjects
%              columns: models (1..Nk)
% alpha0    priors of p(r)
% 
% OUTPUT:
% F_samp  -  sampling estimate of <ln p(y_n|r>
% F_bound -  lower bound on lower bound of <ln p(y_n|r>
% 
% REFERENCE: See appendix in
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (under review)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_F.m 2507 2008-11-30 14:45:22Z klaas $


alpha0 = sort(alpha0);
if alpha0(1) ~= alpha0(end)
    error('Error in function spm_BMS_F: alpha0 should have identical values.')
end
alpha0 = alpha0(1);

a_sum    = sum(alpha);
psi_sum  = psi(a_sum);
psi_diff = psi(alpha) - psi_sum;
gm       = gammaln(alpha);

[s_samp,s_bound] = spm_BMS_F_smpl(alpha,lme,alpha0);

K = length(alpha);
F = 0;
for k = 1:K,
    F = F - (alpha(k) - alpha0)*psi_diff(k) + gm(k);
end
F = F - gammaln(a_sum);

F_bound = F + s_bound;
F_samp  = F + s_samp;

return
end

function [s_samp,s_bound] = spm_BMS_F_smpl (alpha,lme,alpha0)
% Get sample and lower bound approx. for model evidence p(y|r)
% in group BMS; see spm_BMS_F.
% 
% FORMAT [s_samp,s_bound] = spm_BMS_F_smpl (alpha,lme,alpha0)
% 
% REFERENCE: See appendix in
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (under review)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_F_smpl.m 2626 2009-01-20 16:30:08Z maria $


% prevent numerical problems 
max_val = log(realmax('double'));
for i=1:size(lme,1),
        lme(i,:) = lme(i,:) - mean(lme(i,:));
        for k = 1:size(lme,2),
            lme(i,k) = sign(lme(i,k)) * min(max_val,abs(lme(i,k)));
        end
end

% Number of samples per alpha bin (0.1)
Nsamp = 1e3;

% Sample from univariate gamma densities then normalise
% (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
% 209-230)
Nk = length(alpha);
for k = 1:Nk,
    alpha_samp(:,k) = gamrnd(alpha(k),1,Nsamp,1);
end

Ni = size(lme,1);
for i = 1:Ni,
    s_approx(i) = sum((alpha./sum(alpha)).*lme(i,:));
    
    s(i) = 0;
    for n = 1:Nsamp,
        s(i) = s(i) + si_fun(alpha_samp(n,:),lme(i,:));
    end
    s(i) = s(i)/Nsamp;
end

s_bound = sum(s_approx);

s_samp = sum(s);

return


%=========================================================================
function [si] = si_fun (alpha,lme)
% Check a lower bound
% FORMAT [si] = si_fun (alpha,lme)

esi = sum((exp(lme).*alpha)/sum(alpha));
si  = log(esi);

return
end
end

function f = spm_Dpdf(x,a)
% Probability Density Function (PDF) of Dirichlet distribution
% FORMAT f = spm_Dpdf(x,a)
% 
% x - Dirichlet variate
% a - Dirichlet parameters (a>0)
% f - PDF of Dirichlet-distribution at point x
%__________________________________________________________________________
%
% spm_Dpdf implements the Probability Density Function for Dirichlet 
% distribution.
%
% Definition:
%--------------------------------------------------------------------------
% See http://en.wikipedia.org/wiki/Dirichlet_distribution
%
% Algorithm:
%--------------------------------------------------------------------------
% Direct computation using logs and MATLAB's implementation of the log of 
% the gamma function (gammaln).
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_Dpdf.m 4182 2011-02-01 12:29:09Z guillaume $

%-Check enough arguments
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

%-Computation
%--------------------------------------------------------------------------
a = a(:);
x = x(:);

f = exp( gammaln(sum(a)) + sum((a-1).*log(x+eps)) - sum(gammaln(a)) );
end

function [bor,F0,F1] = spm_BMS_bor(L,posterior,priors,C)
% Compute Bayes Omnibus Risk
% FORMAT [bor,F0,F1] = spm_BMS_bor(L,posterior,priors,C)
%
% L         Log model evidence table (models x  subjects)
% posterior .a model counts, .r model-subject probs
% priors    .a model counts
% C         if this field is specified then BOR under family prior 
%           is computed, otherwise BOR under model prior is computed.
%           C(k,f) = 1 if model k belongs to family f (0 otherwise)
%
% REFERENCES:
%
% Rigoux, L, Stephan, KE, Friston, KJ and Daunizeau, J. (2014)
% Bayesian model selection for group studies - Revisited. 
% NeuroImage 84:971-85. doi: 10.1016/j.neuroimage.2013.08.065
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_bor.m 6444 2015-05-21 11:15:48Z guillaume $


if nargin < 4
    options.families = 0;
    % Evidence of null (equal model freqs)
    F0 = FE_null(L,options); 
else
    options.families = 1;
    options.C = C;
    % Evidence of null (equal model freqs) under family prior
    [tmp,F0] = FE_null(L,options); 
end

% Evidence of alternative
F1 = FE(L,posterior,priors); 

% Implied by Eq 5 (see also p39) in Rigoux et al.
% See also, last equation in Appendix 2
bor = 1/(1+exp(F1-F0)); 
end

function [F,ELJ,Sqf,Sqm] = FE(L,posterior,priors)
% derives the free energy for the current approximate posterior
% This routine has been copied from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% and was written by Lionel Rigoux and J. Daunizeau
%
% See equation A.20 in Rigoux et al. (should be F1 on LHS)

[K,n] = size(L);
a0 = sum(posterior.a);
Elogr = psi(posterior.a) - psi(sum(posterior.a));
Sqf = sum(gammaln(posterior.a)) - gammaln(a0) - sum((posterior.a-1).*Elogr);
Sqm = 0;
for i=1:n
    Sqm = Sqm - sum(posterior.r(:,i).*log(posterior.r(:,i)+eps));
end
ELJ = gammaln(sum(priors.a)) - sum(gammaln(priors.a)) + sum((priors.a-1).*Elogr);
for i=1:n
    for k=1:K
        ELJ = ELJ + posterior.r(k,i).*(Elogr(k)+L(k,i));
    end
end
F = ELJ + Sqf + Sqm;
end

function [F0m,F0f] = FE_null (L,options)
% Free energy of the 'null' (H0: equal frequencies)
%
% F0m       Evidence for null (ie. equal probs) over models 
% F0f       Evidence for null (ie. equal probs) over families
%
% This routine derives from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% written by Lionel Rigoux and J. Daunizeau
%
% See Equation A.17 in Rigoux et al.

[K,n] = size(L);
if options.families
    f0 = options.C*sum(options.C,1)'.^-1/size(options.C,2);
    F0f = 0;
else
    F0f = [];
end
F0m = 0;
for i=1:n
    tmp = L(:,i) - max(L(:,i));
    g = exp(tmp)./sum(exp(tmp));
    for k=1:K
        F0m = F0m + g(k).*(L(k,i)-log(K)-log(g(k)+eps));
        if options.families
            F0f = F0f + g(k).*(L(k,i)-log(g(k)+eps)+log(f0(k)));
        end
    end
end
end

