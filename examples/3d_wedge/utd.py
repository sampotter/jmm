import numpy as np
import scipy.special

def F(x):
    sqrtx = np.sqrt(x)
    return 2j*sqrtx*np.exp(1j*x)*scipy.special.modfresnelm(sqrtx)[0]

def D(sign, n_o, t_o, t_e, k, n, t_in, t_out, D2T, s):
    t_in_perp = t_in - np.outer(t_in@t_e, t_e)
    t_out_perp = t_out - np.outer(t_out@t_e, t_e)

    beta0 = np.arccos(t_out@t_e)
    phi_in = np.arctan2(t_in_perp@n_o, t_in_perp@t_o)
    phi_out = np.arctan2(t_out_perp@n_o, t_out_perp@t_o)

    t_in_perp = t_in - np.outer(t_in@t_e, t_e)
    t_out_perp = t_out - np.outer(t_out@t_e, t_e)

    beta0 = np.arccos(t_out@t_e)
    phi_in = np.arctan2(t_in_perp@n_o, t_in_perp@t_o)
    phi_out = np.arctan2(t_out_perp@n_o, t_out_perp@t_o)

    L = s*(rho_e + s)*rho_1*rho_2/(rho_e*(rho_1 + s)*(rho_2 + s))
    L *= np.sin(beta0)**2

    def N(sign, beta):
        I = np.array([-1, 0, 1], dtype=int)
        return I[np.argmin(abs(beta + sign*np.pi - 2*np.pi*n*I))]

    def a(sign, beta):
        return 2*np.cos(2*np.pi*n*N(sign, beta) - beta)**2

    def Di(sign1, sign2):
        _ = -np.exp(1j*np.pi/4)/(2*n*np.sqrt(2*np.pi*k)*np.sin(beta0))
        _ /= np.tan((np.pi + sign1*(phi_out + sign2*phi_in))/(2*n))
        _ *= F(k*L*a(sign1, phi_out + sign2*phi_in))
        return _

    D1 = Di(+1, -1)
    D2 = Di(-1, -1)
    D3 = Di(+1, +1)
    D4 = Di(-1, +1)

    return D1 + D2 + sign*(D3 + D4)
