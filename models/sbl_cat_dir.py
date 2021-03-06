import os
import argparse
import numpy as np
from utils.helpers import *


class SBL_Cat_Dir():
    """
    DESCRIPTION: Categorical Dirichlet Bayesian Sequential Learner
        * Learner parses a categorical sequence
        * Updates a conjugate Cat-Dirichlet posterior with new evidence
        * Calculates different surprise measures as the events come in
    INPUT:  Sequence, type: estimated statistics (SP, AP, TP), order: order of TP inference,
    catch: in- or exclusion of catch trials, tau: exponential forgetting parameter
    OUTPUT: Predictive surprisal, Bayesian surprisal, Confidence-corrected surprisal
    [t, o_t, s_t, Prediction_Surprise, Bayesian_Surprise, Confidence_Corrected_Surprise]
    """
    def __init__(self, seq, hidden, tau, model_type, order, catch, verbose=False):
        self.catch = True
        # Initialize SBL-learned sequence and exponential forgetting parameter
        self.sequence = seq.astype(int)
        # remove catch trials from observation sequence
        if self.catch == False:
            catch_ind = np.argwhere(self.sequence==2)
            self.sequence = np.delete(self.sequence, catch_ind)

        self.hidden = hidden
        self.T = len(self.sequence)
        self.order = order
        self.type = model_type
        self.tau = tau
        self.verbose = verbose
        self.no_obs = np.unique(self.sequence).shape[0]
        self.stim_ind = np.zeros((self.T, self.no_obs))

        # Construct matrix where col represents binary ind of specific stim at t
        for t in range(self.T):
            self.stim_ind[t, self.sequence[t]] = 1

        # AP: Generate T-dim vector indicating no-alternation from t-1 to t
        self.repetition = np.zeros(self.T)
        for t in range(1, self.T):
            if self.sequence[t] == self.sequence[t-1]:
                self.repetition[t] = 1

        # TP: Generate T-dim vectors indicating transition from state i
        if self.order == 1:
            self.transitions = np.zeros((self.T, self.no_obs))
            for t in range(1, self.T):
                self.transitions[t, 0] = (self.sequence[t-1] == 0)
                self.transitions[t, 1] = (self.sequence[t-1] == 1)
                if self.catch == True:
                    self.transitions[t, 2] = (self.sequence[t-1] == 2)

        # TP: For second order also include transitions from t-2
        elif self.order == 2:
            self.transitions = np.zeros((self.T, self.no_obs, self.order))
            for t in range(2, self.T):
                self.transitions[t, 0, 0] = (self.sequence[t-1] == 0)
                self.transitions[t, 0, 1] = (self.sequence[t-2] == 0)
                self.transitions[t, 1, 0] = (self.sequence[t-1] == 1)
                self.transitions[t, 1, 1] = (self.sequence[t-2] == 1)
                if self.catch == True:
                    self.transitions[t, 2, 0] = (self.sequence[t-1] == 2)
                    self.transitions[t, 2, 1] = (self.sequence[t-2] == 2)
        else:
            raise "Provide supported order (1st or 2nd)"

        # Generate one T matrix with all discounting values
        self.exp_forgetting = np.exp(-self.tau*np.arange(self.T)[::-1])

        # Initialize prior alpha matrices; Equation 3
        if self.type == "SP":
            self.alphas = np.ones(self.no_obs)
        elif self.type == "AP":
            self.alphas = np.ones(2)
        elif self.type == "TP":
            if self.order == 1:
                self.alphas = np.ones((self.no_obs, self.no_obs))
            elif self.order == 2:
                self.alphas = np.ones((self.no_obs, self.no_obs, self.no_obs))
        else:
            raise "Provide right model type (SP, AP, TP)"

    def update_posterior(self):
        """
        Updates the matrix of alphas coding  the posterior distrution.
        Input: Inference Type (SP, AP, TP)
               Order (1, 2)
               Time point (1:T)
               Prior vector of alphas
        Output:
        """
        # We first compute the vector describing the exponential weighting of events.
        # See equation (8).
        exp_weighting = self.exp_forgetting[-(self.t+1):]

        # Next, this vector is combined with the observed sequence of events.
        if self.type == "SP":
            for i in range(self.no_obs):
                self.alphas[i] = 1 + np.dot(exp_weighting,
                                            self.stim_ind[:self.t+1, i])

        elif self.type == "AP":
            # Be careful 0th entry is set to zero - no repetition possible if only 1 obs!
            if self.t == 0: # Equation 3
                self.alphas[0] = 1
                self.alphas[1] = 1
            else:
                self.alphas[0] = 1 + np.dot(exp_weighting[:self.t+1],
                                            1-self.repetition[:self.t+1])
                self.alphas[1] = 1 + np.dot(exp_weighting[:self.t+1],
                                            self.repetition[:self.t+1])

        elif self.type == "TP":
            if self.t < self.order: # Equation 3
                if self.order == 1:
                    self.alphas = np.ones((self.no_obs, self.no_obs))
                elif self.order == 2:
                    self.alphas = np.ones((self.no_obs, self.no_obs, self.no_obs))
            else:
                for i in range(self.no_obs):
                    for j in range(self.no_obs):
                        # First order
                        if self.order == 1:
                             # from-to alphas
                             self.alphas[i, j] = 1 + np.dot(exp_weighting,
                             self.stim_ind[:self.t+1, j]*self.transitions[:self.t+1, i])
                        # Second order
                        elif self.order == 2:
                            for k in range(self.no_obs):
                                self.alphas[k, i, j] = 1 + np.dot(exp_weighting,
                                self.stim_ind[:self.t+1, j]*self.transitions[:self.t+1, i, 0]*self.transitions[:self.t+1, k, 1])

    def posterior_predictive(self, alphas):
        # Used for PS computation; Equation 6
        return np.array([alpha/alphas.sum(axis=0) for alpha in alphas])

    def naive_posterior(self, ind):
        # Performs single step update on naive/flat prior by adding 1 to alpha
        # corresponding to observed event. Used for CS computation.
        naive_update = self.alpha_naive_init.copy()
        naive_update[ind] += 1
        return naive_update

    # Surprise readout computation
    def predictive_surprisal(self, alphas, ind):
        # Equation 9
        return -np.log(self.posterior_predictive(alphas)[ind])

    def bayesian_surprisal(self, alphas_old, alphas):
        # Equation 10
        return kl_dir(alphas_old, alphas)

    def corrected_surprisal(self, alphas_old, ind):
        # Equation 11
        return kl_dir(alphas_old, self.naive_posterior(ind))

    def compute_surprisal(self, max_T, verbose_surprisal=False):
        if verbose_surprisal:
            print("{}: Computing different surprisal measures for {} timesteps.".format(self.type, max_T))
        results = []

        self.alpha_naive_init = self.alphas.copy()
        for t in range(max_T):
            # Loop over the full sequence and compute surprisal iteratively
            alphas_old = self.alphas.copy()
            self.t = t
            self.update_posterior()

            if self.type == "SP":
                ind = int(self.sequence[self.t])
            elif self.type == "AP":
                ind = int(self.repetition[self.t])
            elif self.type == "TP":
                # from and to stimulus transition
                if self.order == 1:
                    ind = (np.argmax(self.transitions[self.t, :]), np.argmax(self.stim_ind[self.t, :]))
                elif self.order == 2:
                    ind = (np.argmax(self.transitions[self.t, :, 1]), np.argmax(self.transitions[self.t, :, 0]), np.argmax(self.stim_ind[self.t, :]))
            else:
                raise "Provide right model type (SP, AP, TP)"

            # Compute surprise for y_t
            PS_temp = self.predictive_surprisal(alphas_old, ind)
            BS_temp = self.bayesian_surprisal(alphas_old, self.alphas)
            CS_temp = self.corrected_surprisal(alphas_old, ind)

            if verbose_surprisal:
                print("{} - t={}: PS={}, BS={}, CS={}".format(self.type, t+1, round(PS_temp, 4),  round(BS_temp, 4), round(CS_temp, 4)))

            # Save t, observation, regime, and surprise values
            temp = [t, self.sequence[t], self.hidden[t], PS_temp, BS_temp, CS_temp]
            distr_params = list(self.alphas.reshape(1, -1)[0])
            results.append(temp + distr_params)

        print("{}: Done computing surprisal measures for all {} timesteps.".format(self.type, self.T))
        return np.asarray(results)


def main(seq, hidden, tau, model_type, order, catch,
         save_results=False, title="temp", verbose=False):
    # Compute Surprisal for all time steps for DirCat Model
    CD_SBL_temp = SBL_Cat_Dir(seq, hidden, tau, model_type, order, catch, verbose)
    results = CD_SBL_temp.compute_surprisal(max_T=CD_SBL_temp.T, verbose_surprisal=verbose)

    # Prepare to save formatted results
    time = results[:, 0]
    sequence = results[:, 1]
    hidden = results[:, 2]
    PS = results[:, 3]
    BS = results[:, 4]
    CS = results[:, 5]

    if save_results:
        results_formatted = {"time": time,
                             "sequence": sequence,
                             "hidden": hidden,
                             "predictive_surprise": PS,
                             "bayesian_surprise": BS,
                             "confidence_corrected_surprise": CS}

        save_obj(results_formatted, results_dir + title)
        #print("Saved in File: {}".format(results_dir + title))
    else:
        return PS, BS, CS

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', '--sample_file', action="store",
                        default="temporary_sample_title", type=str,
                        help='Title of file in which sequence in stored')
    parser.add_argument('-tau', '--forget_param', action="store",
                        default=0., type=float,
                        help='Exponentially weighting parameter for memory/posterior updating')
    parser.add_argument('-model', '--model', action="store", default="SP",
                        type=str,
                        help='Categorical Dirichlet Probability Model (SP, AP, TP)')
    parser.add_argument('-order', '--order', action="store", default=1, type=int,
                        help='Order for Transition Probability model (1 or 2)')
    parser.add_argument('catch', '--catch_trials', action="store_true", default=False,
                        help='Exclusion (0) or Inclusion (1) of catch trials')
    parser.add_argument('-pkl_in', '--pickle', action="store_true", help='Load matlab sequence file.')
    parser.add_argument('-S', '--save', action="store_true", default=False, help='Save results to array.')
    parser.add_argument('-V', '--verbose',
                        action="store_true",
                        default=False,
						help='Get status printed out')
    args = parser.parse_args()

    if args.pickle:
        sample, meta = load_obj(results_dir + args.sample_file + ".pkl")
    else:
        sample = load_obj(results_dir + args.sample_file + ".mat")

    # Fetch sequence (observations) and regime (hidden state) information
    seq = sample[:, 2]
    hidden = sample[:, 1]

    tau = args.forget_param
    model = args.model
    save_results = args.save
    order = args.order
    catch = args.catch
    v = args.verbose

    main(seq, hidden, tau, model, order, catch,
         save_results=save_results, title="CD_" + model + "_" + args.sample_file,
         verbose=v)
