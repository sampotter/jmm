import logging
import matplotlib.pyplot as plt
import numpy as np

from utd_wedge_problem import UtdWedgeProblem

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    log = logging.getLogger('make_errors_plots.py')

    n = 2/3
    sp = 5/np.sqrt(2)
    phip = 2*np.pi/3
    w = 10
    h = 2
    R = 1
    r = 1.5
    omega = 10000

    target_edge_lengths = [1, 0.75, 0.5, 0.25, 0.2, 0.15, 0.1]

    probs = dict()

    max_edge_lengths = []
    avg_edge_lengths = []

    MAE_pt_src_tau = []
    MAE_near_refl_tau = []
    MAE_far_refl_tau = []
    MAE_diff_tau = []

    RMSE_pt_src_tau = []
    RMSE_near_refl_tau = []
    RMSE_far_refl_tau = []
    RMSE_diff_tau = []

    for target_edge_length in target_edge_lengths:
        log.info('target_edge_length: %g', target_edge_length)

        maxvol = target_edge_length**3
        log.info('maxvol: %g', maxvol)

        prob = UtdWedgeProblem(maxvol, n, sp, phip, w, h, R, r, omega)
        probs[target_edge_length] = prob # just in case

        max_edge_lengths.append(prob.max_edge_length)
        avg_edge_lengths.append(prob.avg_edge_length)

        MAE_pt_src_tau.append(prob.compute_MAE(prob.error_pt_src_tau))
        MAE_near_refl_tau.append(prob.compute_MAE(prob.error_near_refl_tau))
        MAE_far_refl_tau.append(prob.compute_MAE(prob.error_far_refl_tau))
        MAE_diff_tau.append(prob.compute_MAE(prob.error_diff_tau))

        RMSE_pt_src_tau.append(prob.compute_RMSE(prob.error_pt_src_tau))
        RMSE_near_refl_tau.append(prob.compute_RMSE(prob.error_near_refl_tau))
        RMSE_far_refl_tau.append(prob.compute_RMSE(prob.error_far_refl_tau))
        RMSE_diff_tau.append(prob.compute_RMSE(prob.error_diff_tau))
