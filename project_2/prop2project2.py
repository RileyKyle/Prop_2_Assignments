# %%
import numpy as np
import scipy
import matplotlib.pyplot as plt

# %%
gamma = 1.4
Cp = 1004.5
R = 287
h = 4.45e+7

class r:
    M0 = 0.85
    T0 = 236.23 # K
    P0 = 35651 # Pa
    M2 = 0.5 
    A1 = 0.00318 # m^2
    A0 = A1
    
    pi_c = 14.43
    pi_b = 1.0
    pi_d = 0.989163
    eta_t = 0.90
    
    Tt4 = 1500 # K


def eta_c_fit(Tt2,Tt4):
    x = Tt4/Tt2
    a = -2.3219 
    b = 1.7094 
    c = -0.2897 
    d = 0.0155
    return a + b*x + c*x**2 + d*x**3

def corr_rpm_fit(pi_c):
    x = pi_c
    aa = 22352.47 
    bb = 3718.445 
    cc = -52.5576 
    dd = -1.0899
    return aa + bb*x + cc*x**2 + dd*x**3


# %% [markdown]
# # On Design

# %%
r.a0 = np.sqrt(gamma*R*r.T0)
r.u0 = r.a0 * r.M0

r.rho0 = r.P0/(R*r.T0)
r.m_dot_air = r.rho0 * r.u0 * r.A0

# %%
r.tau_r = (1+(gamma-1)/2 * r.M0**2)
r.pi_r = r.tau_r**(gamma/(gamma-1))
r.u1 = r.u0

# %% [markdown]
# ## inlet

# %%
r.tau_d = 1
r.u2 = r.u1

# %% [markdown]
# ## compressor

# %%
r.tau_c = r.pi_c **(gamma-1)/gamma
r.u3 = r.u2

# %% [markdown]
# ## burner

# %%
r.tau_lambda = r.Tt4 / r.T0
r.tau_b = r.tau_lambda / (r.tau_r * r.tau_d * r.tau_c)
r.u4 = r.u3

# %% [markdown]
# ## turbine

# %%
r.tau_t = 1-(r.tau_c - 1)*(r.tau_r/r.tau_lambda)
r.pi_t = (1-1/r.eta_t*(1-r.tau_t))**(gamma/(gamma-1))
r.u5 = 0.8 * r.u4

# %% [markdown]
# ## nozzle

# %%
r.tau_n = 1
r.pi_n = 1

# %%
r.Tt0 = r.T0 * r.tau_r
r.Tt1 = r.Tt0
r.Tt2 = r.Tt1 * r.tau_d
r.Tt3 = r.Tt2 * r.tau_c
# r.Tt4 is given
r.Tt5 = r.Tt4 * r.tau_t
r.Tt8 = r.Tt5 * r.tau_n

# %%
r.Pt0 = r.P0 * r.pi_r
r.Pt1 = r.Pt0
r.Pt2 = r.Pt1 * r.pi_d
r.Pt3 = r.Pt2 * r.pi_c
r.Pt4 = r.Pt3 * r.pi_b
r.Pt5 = r.Pt4 * r.pi_t
r.Pt8 = r.Pt5 * r.pi_n

# %%
r.T8 = r.Tt8 / (1+(gamma-1)/2*r.M2**2)
r.P8 = r.Pt8 / (1+(gamma-1)/2*r.M2**2)**(gamma/(gamma-1))


r.u8 = r.M2 * np.sqrt(gamma*R*r.T8)
r.rho8 = r.P8 / (r.T8*R)
r.A8 = r.m_dot_air / (r.rho8 * r.u8)
r.F = (r.m_dot_air) * (r.u8 - r.u0) + (r.P8 - r.P0) * r.A8

# %% [markdown]
# # Off Design

# %%


# %%



