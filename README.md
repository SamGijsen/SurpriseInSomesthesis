# Neural Surprise In Somesthesis

## Code accompanying the paper 'Neural surprise in somatosensory Bayesian learning'

## Authors: Sam Gijsen*, Miro Grundei*, Robert Tjarko Lange, Dirk Ostwald, Felix Blankenburg

### * Equal contribution


/models/ in python3 with dependencies:


/analysis/ in MATLAB with dependencies:
+ SPM12

## Repository Structure
```
NeuralSurpriseInSomesthesis
├── models: The compared learning models
    +- cat_dir.py: Dirichlet-Categorical model
    +- hmm.py: Hidden Markov model
├── sbl_agents: Different sequential bayesian learning agents
    +- sbl_cat_dir: Categorical-Dirichlet SBL agent
    +- sbl_hmm: Hidden Markov Model Agent
    ├── utils: Helper files
        +- helpers.py: Define helper functions for loading files, visualization and post-processing
├── analysis: Various scripts used for EEG analysis
    +- sbl_ffx.m: fixed-effects (FFX) Bayesian model selection
    +- sbl_rfx.m: random-effects (RFX) Bayesian model selection
    +- sbl_modelrecovery_ffx.m: FFX model recovery using simulations
    +- sbl_modelrecovery_rfx.m: RFX model recovery using simulations
├── README.md: Documentation
├── requirements.txt: Dependencies

