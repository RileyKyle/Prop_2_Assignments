import numpy as np
from scipy.optimize import fsolve

# Constants
R = 287
Cp = 1004.5
y = 1.4
h = 4.45E7
T_STP = 288

pi_b = 1
pi_d = 0.989163
A_1 = 0.00318
mdotf_ondesign = 0.008726071272118822
M_4 = 1
M_8 = 1
A4 = 0.0005141364620689669
A8 = 0.0012745633749949147

# pi_t/tau_t are constant because A8/A4 and n_t are constant
A8_A4 = A8 / A4
n_t = 0.90
pi_t = 0.35368326100748115
tau_t = 0.7687676014665914

# ECS Limits
picmax = 18
picmin = 7.5
RPMmax = 66000
Tt4max = 1800

# Flight Characteristics On-Design
M_0 = np.array([0.85])
mdotf = np.array([mdotf_ondesign])
T_0 = 236.23
P_0 = 35651

# Initializing Systems of Equations Answers
Tt_4 = np.zeros((len(mdotf), len(M_0)))
tau_comp = np.zeros((len(mdotf), len(M_0)))
pi_c = np.zeros((len(mdotf), len(M_0)))
n_comp = np.zeros((len(mdotf), len(M_0)))
m_dot = np.zeros((len(mdotf), len(M_0)))
spill = np.zeros((len(mdotf), len(M_0)))
RPM = np.zeros((len(mdotf), len(M_0)))
rpm_answer =  np.zeros((len(mdotf), len(M_0)))

def equations(vars):
    Tt4, pic, mdot, rpm = vars
    rho_0 = P_0 / (R * T_0)

    u_0 = M_0 * (y * R * T_0) ** 0.5
    mdot_NS = rho_0 * u_0 * A_1
    Tt_0 = T_0 * (1 + ((y - 1) / 2) * M_0 ** 2)
    Pt_0 = P_0 * (1 + ((y - 1) / 2) * M_0 ** 2) ** (y / (y - 1))

    Tt_2 = Tt_0
    Pt_2 = pi_d * Pt_0

    tau_c = 1 + (Tt4 / Tt_0) * (1 - tau_t)
    aa = 22352.47
    bb = 3718.445
    cc = -52.5576
    dd = -1.0899
    RPM = (aa + bb * (Tt4 / Tt_2) + cc * (Tt4 / Tt_2) ** 2 + dd * (Tt4 / Tt_2) ** 3)
    pi_r = (1 + (y - 1) / 2 * M_0 ** 2) ** (y / (y - 1))

    n_c = -2.3219 + 1.7094 * (Tt4 / Tt_2) + -0.2897 * (Tt4 / Tt_2) ** 2 + 0.0155 * (Tt4 / Tt_2) ** 3
    n_c2 = (pic ** ((y - 1) / y) - 1) / (tau_c - 1)
    mdotf_eq = (mdot * Cp * (Tt4 - Tt_0 * tau_c)) / h
    G = y ** (1 / 2) * (2 / (y + 1)) ** ((y + 1) / 2 / (y - 1))
    mdot_eq = G / np.sqrt(R) * A4 * P_0 * pi_r * pi_d * pi_b * pic / np.sqrt(Tt4)
    rpm_eq = np.sqrt(Tt_0 / T_STP) * RPM

    return [(n_c - n_c2)[0],(mdotf_eq - mdotf)[0],(mdot_eq - mdot)[0],(rpm_eq - rpm)[0] ]


for i in range(1):
    for j in range(1):
        initial_guess = [1000, 10, .5, 1000]
        x = fsolve(equations, initial_guess)
        Tt4, pic, mdot, rpm = x
        Tt_4[j, i] = Tt4
        pi_c[j, i] = pic
        m_dot[j, i] = mdot
        rpm_answer[j, i] = rpm
        # Station 3
        Pt_3 = pi_c[j, i] * Pt_2

        # Station 4
        T_4 = Tt_4[j, i] / (1 + ((y - 1) / 2) * M_4 ** 2)
        Pt_4 = Pt_3
        P_4 = Pt_4 / (1 + ((y - 1) / 2) * M_4 ** 2) ** (y / (y - 1))
        rho_4 = P_4 / (R * T_4)
        a_4 = (y * R * T_4) ** 0.5
        u_4 = M_4 * a_4

        # Station 5
        Tt_5 = tau_t * Tt_4[j, i]
        Pt_5 = pi_t * Pt_4
        # Station 8
        Tt_8 = Tt_5
        T_8 = Tt_8 / (1 + ((y - 1) / 2) * M_8 ** 2)
        u_8 = M_8 * ((y * R * T_8) ** 0.5)
        rho_8 = (rho_4 * u_4) / (u_8 * A8_A4)
        P_8 = rho_8 * R * T_8

        # Plotting Analysis
        Thrust[j, i] = m_dot[j, i] * (u_8 - u_0) + (P_8 - P_0) * A8
        spill[j, i] = mdot_NS - m_dot[j, i]