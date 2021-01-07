import os
import argparse
import numpy as np
import math
import random
from hmmlearn import hmm
from utils.helpers import *


class SBL_HMM():
    """
    DESCRIPTION: Hidden Markov Model Bayesian Sequential Learner
        * Learner parses a categorical sequence
        * Updates an HMM posterior with new evidence
        * Calculates different surprise measures as the events come in
    INPUT: Sequence of observations, catch: inclusion of catch trial,
    type: estimated statistic (SP,AP,TP), n_states: number of latent states of HMM,
    fix_tm: whether transition matrix is estimated
    OUTPUT: Predictive surprisal, Bayesian surprisal, Confidence-corrected surprisal
    [t, o_t, s_t, Prediction_Surprise, Bayesian_Surprise, Confidence_Corrected_Surprise]
    """
    def __init__(self, seq, hidden, n_states, model_type, fix_tm, catch, order, verbose):
        # Initialize SBL-learned sequence and exponential forgetting parameter
        self.sequence = seq.astype(int)
        self.catch = catch
        self.order = order
        # remove catch trials from observation sequence
        if self.catch == False:
            catch_ind = np.argwhere(self.sequence==2)
            self.sequence = np.delete(self.sequence, catch_ind)

        self.hidden = hidden
        self.T = len(self.sequence)
        self.type = model_type
        self.n_states = n_states
        self.fix_tm = fix_tm
        self.verbose = verbose
        self.no_obs = np.unique(self.sequence).shape[0]
        self.stim_ind = np.zeros((self.T, self.no_obs))

        # SP: Construct matrix where col represents binary index of specific stim at t
        for t in range(self.T):
            self.stim_ind[t, self.sequence[t]] = 1

        # AP: Generate T-dim vector indicating no-alternation from t-1 to t
        self.repetition = np.zeros(self.T)
        for t in range(1, self.T):
            if self.sequence[t] == self.sequence[t-1]:
                self.repetition[t] = 1

        # TP: Generate T-dim vectors indicating transition from state i
        self.all_transitions = np.zeros(self.T)
        for t in range(1, self.T):
            # All transitions, with catch trials
            if self.catch == True:
                if self.sequence[t-1] == 0 and self.sequence[t] == 0:
                    self.all_transitions[t] = 0
                if self.sequence[t-1] == 0 and self.sequence[t] == 1:
                    self.all_transitions[t] = 1
                if self.sequence[t-1] == 0 and self.sequence[t] == 2:
                    self.all_transitions[t] = 2
                if self.sequence[t-1] == 1 and self.sequence[t] == 0:
                    self.all_transitions[t] = 3
                if self.sequence[t-1] == 1 and self.sequence[t] == 1:
                    self.all_transitions[t] = 4
                if self.sequence[t-1] == 1 and self.sequence[t] == 2:
                    self.all_transitions[t] = 5
                if self.sequence[t-1] == 2 and self.sequence[t] == 0:
                    self.all_transitions[t] = 6
                if self.sequence[t-1] == 2 and self.sequence[t] == 1:
                    self.all_transitions[t] = 7
                if self.sequence[t-1] == 2 and self.sequence[t] == 2:
                    self.all_transitions[t] = 8
                self.number_of_transitions = self.no_obs**2

            # Same procedure for a binary sequence (without catch trials)
            elif self.catch == False:
                # First Order transitions
                if self.order == 1:
                    if self.sequence[t-1] == 0 and self.sequence[t] == 0:
                        self.all_transitions[t] = 0
                    if self.sequence[t-1] == 0 and self.sequence[t] == 1:
                        self.all_transitions[t] = 1
                    if self.sequence[t-1] == 1 and self.sequence[t] == 0:
                        self.all_transitions[t] = 2
                    if self.sequence[t-1] == 1 and self.sequence[t] == 1:
                        self.all_transitions[t] = 3
                    self.number_of_transitions = self.no_obs**2

                # Second Order transitions
                elif self.order == 2:
                    if self.sequence[t-2] == 0 and self.sequence[t-1] == 0 and self.sequence[t] == 0:
                        self.all_transitions[t] = 0
                    if self.sequence[t-2] == 0 and self.sequence[t-1] == 0 and self.sequence[t] == 1:
                        self.all_transitions[t] = 1
                    if self.sequence[t-2] == 0 and self.sequence[t-1] == 1 and self.sequence[t] == 0:
                        self.all_transitions[t] = 2
                    if self.sequence[t-2] == 1 and self.sequence[t-1] == 0 and self.sequence[t] == 0:
                        self.all_transitions[t] = 3
                    if self.sequence[t-2] == 0 and self.sequence[t-1] == 1 and self.sequence[t] == 1:
                        self.all_transitions[t] = 4
                    if self.sequence[t-2] == 1 and self.sequence[t-1] == 0 and self.sequence[t] == 1:
                        self.all_transitions[t] = 5
                    if self.sequence[t-2] == 1 and self.sequence[t-1] == 1 and self.sequence[t] == 0:
                        self.all_transitions[t] = 6
                    if self.sequence[t-2] == 1 and self.sequence[t-1] == 1 and self.sequence[t] == 1:
                        self.all_transitions[t] = 7
                    self.number_of_transitions = self.no_obs**3
            else:
                raise "Order Not Supported (1st or 2nd)"


    def init_hmm(self):
        """
        Input: number of desired hidden states for HMM model.
        Output: Initialized matrices for HMM training.
        """
        startprob = np.repeat(1./self.n_states, self.n_states)

        # No perfect uniformity is used as this may lead to instabilities
        # in the hidden state inference. As such, we initialize close-to-uniform.
        if self.type == "SP":
            temp = np.repeat(1./self.no_obs, self.no_obs)
            if self.catch == False:
                emissionprob = np.array([[0.495, 0.505], [0.505, 0.495]])
            elif self.catch == True:
                emissionprob = np.array([[0.333, 0.334, 0.333], [0.334, 0.333, 0.333]])

        elif self.type == "AP":
            temp = np.repeat(1./2, 2)
            emissionprob = np.array([[0.49, 0.51], [0.51, 0.49]])

        elif self.type == "TP":
            if self.catch == False:
                if self.order == 1:
                    emissionprob = np.array([[0.49, 0.51, 0.51, 0.49], [0.51, 0.49, 0.49, 0.51]])
                elif self.order == 2:
                    emissionprob = np.array([[0.1245, 0.1245, 0.1255, 0.125, 0.1255, 0.1255, 0.125, 0.1245],
                    [0.1255, 0.1245, 0.1255, 0.125, 0.1245, 0.1245, 0.125, 0.1255]])

        # In case the transition matrix is fixed we provide the values,
        # otherwise it'll be estimated.
        if self.fix_tm == False:
            transmat = np.repeat([startprob], self.n_states, axis=0)
        elif self.fix_tm == True:
            # fix transition matrix close to true values used for sequence generation
            transmat = np.array([[0.99, 0.01], [0.01, 0.99]])

        return startprob, transmat, emissionprob

    def calc_all_posteriors(self, t, random_state):
        """
        Input: Unique state id transformed data, length per ep trace and HMM inits
        Output: The trained HMM model with all attributes
        - Option: If desired get AIC and BIC for model as well
        """
        startprob, transmat, emissionprob = self.init_hmm()

        # iteration parameters: number of iterations, epsilon (minimum change in fit), verbosity
        n_iter = 20000
        eps = 1e-10
        v = False

        # Initialization of HMM
        if self.fix_tm == False: # Fit transition matrix
                model = hmm.MultinomialHMM(n_components=self.n_states,
                random_state=random_state, n_iter=n_iter, tol=eps, init_params="")

        elif self.fix_tm == True: # Do not fit transition matrix (thus fixed to prior)
                model = hmm.MultinomialHMM(n_components=self.n_states,
                random_state=random_state, n_iter=n_iter, tol=eps, params="e", init_params="")

        model.startprob_ = startprob
        model.transmat_ = transmat
        model.emissionprob_ = emissionprob

        # The HMM requires each unique event to occur in the supplied sequence
        # and thus we concatenate a vector containing each unique event once with the actual sequence.
        if self.type == "SP":
            valid_sample = np.unique(self.sequence)
            temp = self.sequence[:t+1].reshape(1, -1).T
            valid_sample = valid_sample.reshape(1,-1).T

        elif self.type == "AP":
            valid_sample = np.array([0, 1])
            temp = self.repetition[:t+1].reshape(1, -1).T.astype(int)
            valid_sample = valid_sample.reshape(1,-1).T

        elif self.type == "TP":
            valid_sample = np.arange(self.number_of_transitions)
            temp = self.all_transitions[:t+1].reshape(1, -1).T.astype(int)
            valid_sample = valid_sample.reshape(1,-1).T

        valid_seq = np.vstack((valid_sample, temp))

        # Afterwards we can fit the HMM on the sequence observed so far
        self.model = model.fit(valid_seq)
        logprob, state_sequence = model.decode(valid_seq)
        logprob, posteriors = model.score_samples(valid_seq)

        return posteriors, state_sequence

    # Functions for computation of three surprise readouts: PS, BS, CS

    # The posterior predictive distribution is used for the computation of predictive surprise.
    def posterior_predictive(self, posterior, emissionprob):
        return np.matmul(emissionprob.T,
                         np.matmul(self.model.transmat_.T,
                                   posterior.T))

    def predictive_surprisal(self, posterior, emissionprob, ind):
        return -np.log(self.posterior_predictive(posterior, emissionprob)[ind])

    def bayesian_surprisal(self, posterior_old, posterior):
        return kl_general(posterior_old, posterior)

    def compute_surprisal(self, verbose_surprisal, max_T):
        print("{}: Computing different surprisal measures for {} timesteps.".format(self.type, max_T))
        results = []

        hmm_init_posterior = np.repeat(1./self.n_states, self.n_states)

        for t in range(self.T):
            # Loop over the full sequence and compute surprise iteratively using the previously defined functions
            if t == 0:
                posterior_old = hmm_init_posterior

                if self.type == "SP":
                    ep_old = np.ones([self.n_states, self.no_obs]) / self.no_obs
                elif self.type == "AP":
                    ep_old = np.ones([self.n_states, 2]) / 2
                elif self.type == "TP":
                    if self.order == 1:
                        ep_old = np.ones([self.n_states, self.no_obs**2]) / self.no_obs**2
                    elif self.order == 2:
                        ep_old = np.ones([self.n_states, self.no_obs**3]) / self.no_obs**3

            posteriors, state_sequence = self.calc_all_posteriors(t, 1)
            posterior = posteriors[-1]

            # Fetch the event that occurred at time t (inference-type dependent)
            if self.type == "SP":
                ind = int(self.sequence[t])
            elif self.type == "AP":
                ind = int(self.repetition[t])
            elif self.type == "TP":
                # from and to stimulus transition
                ind = int(self.all_transitions[t])
            else:
                raise "Provide right model type (SP, AP, TP)"

            PS_temp = self.predictive_surprisal(posterior_old, ep_old, ind) # Equation 12
            BS_temp = self.bayesian_surprisal(posterior_old, posterior)     # Equation 13

            # CS computation; Equation 14
            O_term = 0
            C_term = 0
            for s in range(self.n_states):
                O_term += ep_old[s, ind]
                C_term += posterior_old[s]*np.log(posterior_old[s])

            CS_temp = PS_temp + BS_temp + np.log(O_term) + C_term

            # extended output
            state1_temp = self.model.emissionprob_[0]
            state2_temp = self.model.emissionprob_[1]
            states_temp = posteriors[-1]
            transmat_temp = self.model.transmat_

            # State posterior and emission probabilities are saved for use a t+1
            posterior_old = posterior[:]
            ep_old = self.model.emissionprob_

            # Print computed quantities (debugging)
            if verbose_surprisal:
                print("time:{} --- posterior: {} --- state: {}".format(t, posteriors[-1], state_sequence[-1]))
                print("{} - t={}: PS={}, BS={}, CS={}".format(self.type, t+1, round(PS_temp, 4),  round(BS_temp, 4), round(CS_temp, 4)))

            temp = [t, self.sequence[t], self.hidden[t], PS_temp, BS_temp, CS_temp, state1_temp, state2_temp, states_temp, transmat_temp, state_sequence[-1]]
            results.append(temp)

        print("{}: Done computing surprisal measures for all {} timesteps.".format(self.type, self.T))
        return np.asarray(results)


def main(seq, hidden, n_states, model_type, fix_tm, catch, order,
         save_results, verbose, title="temp"):

    # Compute Surprise for all time steps
    HMM_SBL_temp = SBL_HMM(seq, hidden, n_states, model_type, fix_tm, catch, order, verbose)
    results = HMM_SBL_temp.compute_surprisal(verbose, max_T=HMM_SBL_temp.T, )

    # Save results in formatted manner
    time = results[:, 0]
    sequence = results[:, 1]
    hidden = results[:, 2]
    PS = results[:, 3]
    BS = results[:, 4]
    CS = results[:, 5]
    state1 = results[:, 6]
    state2 = results[:, 7]
    states = results[:, 8]
    transmat = results[:, 9]

    results_formatted = {"time": time,
                         "sequence": sequence,
                         "hidden": hidden,
                         "predictive_surprise": PS,
                         "bayesian_surprise": BS,
                         "confidence_corrected_surprise": CS,
                         "state1": state1,
                         "state2": state2,
                         "states": states,
                         "transmat": transmat}

    if save_results:
        save_obj(results_formatted, results_dir + title)
        print("Saved in File: {}".format(results_dir + title))
    else:
        return PS, BS, CS


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', '--sample_file', action="store",
                        default="temporary_sample_title", type=str,
                        help='Title of file in which sequence in stored')
    parser.add_argument('-states', '--n_states', action="store",
                        default=2, type=int,
                        help='Number of Hidden States in HMM')
    parser.add_argument('-model', '--model', action="store", default="SP",
                        type=str,
                        help='Categorical Dirichlet Probability Model (SP, AP, TP)')
    parser.add_argument('-fix_tm', '--fixed_tm', action="store", default=1,
                        type=int, help='Freely estimated (0) or Fixed (1) transition matrix')
    parser.add_argument('-catch', '--catch_trials', action="store_true", default=False,
                        help='Exclusion (0) or Inclusion (1) of catch trials')
    parser.add_argument('-order', '--order_TP', action="store", default=1, type=int,
                        help='Order for Transition Probability model (1 or 2)')
    parser.add_argument('-pkl_in', '--pickle', action="store_true", help='Load matlab sequence file.')
    parser.add_argument('-V', '--verbose', action="store_true", default=False,
                        help='Get status printed out')
    parser.add_argument('-S', '--save', action="store_true", default=False,
                        help='Save results to array.')
    args = parser.parse_args()

    if args.pickle:
        sample = load_obj(results_dir + args.sample_file + ".pkl")
    else:
        sample = load_obj(results_dir + args.sample_file + ".mat")

    # Fetch sequence (observations) and regime (hidden state) information
    seq = sample[:, 2]
    hidden = sample[:, 1]

    n_states = args.n_states
    model = args.model
    catch = args.catch_trials
    order = args.order_TP
    fix_tm = args.fixed_tm
    save_results = args.save
    v = args.verbose
    s = args.save

    main(seq, hidden, n_states, model, fix_tm, catch, order,
         save_results=s, verbose=v, title="HMM_" + model + "_" + args.sample_file)
