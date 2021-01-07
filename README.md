# Neural Surprise In Somesthesis
Code accompanying the paper 'Neural surprise in somatosensory Bayesian learning'

## Authors: Sam Gijsen*, Miro Grundei*, Robert Tjarko Lange, Dirk Ostwald, Felix Blankenburg

## * Equal contribution


/models/ includes:
- Dirichlet-Categorical model
- Hidden Markov model 
in python3 with dependencies:


/analysis/ includes 
- fixed-effects (FFX) Bayesian model selection for EEG data
- random-effects (RFX) Bayesian model selection for EEG data
- FFX model recovery using simulations
- RFX model recovery using simulations

in MATLAB with dependencies:
+ SPM12

## Repository Structure
```
NeuralSurpriseInSomesthesis
├── figures: Store plots
├── report: Write up and presentations
├── results: Surprise regressors and modeling results
├── models: The compared learning models
    +- seq_gen.py: Samples sequence according to graphical model
    +- seq_analysis.py: Analyze the sampled sequence - get empirical stats
├── sbl_agents: Different sequential bayesian learning agents
    +- sbl_cat_dir: Categorical-Dirichlet SBL agent
    +- sbl_hmm: Hidden Markov Model Agent
    ├── utils: Helper files
        +- helpers.py: Define helper functions for loading files, visualization and post-processing
├── .gitignore: jadajada
├── README.md: Documentation
├── requirements.txt: Dependencies
├── run_in_parallel.txt: Run trial-by-trial analysis in parallel from command line
├── workspace.ipynb: Main workspace - Check it out!
