# %%
import numpy as np
import scipy
import matplotlib.pyplot as plt


# %%
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

# %%
gamma = 1.4
Cp = 1004.5
R = 287
h = 4.45e+7

class given:
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

# %% [markdown]
# # On Design

# %%
u = np.zeros(9)
P = np.zeros(9)
T = np.zeros(9)
M = np.zeros(9)
rho = np.zeros(9)
A = np.zeros(9)
Pt = np.zeros(9)
Tt = np.zeros(9)

# %% [markdown]
# ## Farfield

# %%
M[0] = given.M0
T[0] = given.T0
P[0] = given.P0
rho[0] = P[0]/(R*T[0])
A[0] = given.A0

u[0] = M[0]*np.sqrt(gamma*R*T[0])

tau_r = 1+(gamma-1)/2 * M[0]**2
pi_r = tau_r**(gamma/(gamma-1))

Tt[0] = T[0] * tau_r
Pt[0] = P[0] * pi_r

m_dot = rho[0] * u[0] * A[0]
m_dot

pi_r

# %% [markdown]
# ## inlet

# %%
M[1] = given.M0
T[1] = given.T0
P[1] = given.P0
rho[1] = P[1]/(R*T[1])
A[1] = given.A1

u[1] = M[1]*np.sqrt(gamma*R*T[1])


Tt[1] = Tt[0]
Pt[1] = Pt[0]

# %% [markdown]
# ## Compressor Face

# %%
M[2] = given.M2

tau_d = 1
pi_d = given.pi_d

Tt[2] = Tt[1] * tau_d
Pt[2] = Pt[1] * pi_d

T[2] = Tt[2]/(1+(gamma-1)/2 * M[2]**2)
P[2] = Pt[2]/(1+(gamma-1)/2 * M[2]**2)**(gamma/(gamma-1))

rho[2] = P[2]/(R*T[2])

u[2] = M[2]*np.sqrt(gamma*R*T[2])

A[2] = m_dot / (rho[2] * u[2])

# %% [markdown]
# ## compressor

# %%
u[3] = u[2]

pi_c = given.pi_c
eta_c = eta_c_fit(Tt[2],given.Tt4)
tau_c =  ((pi_c**((gamma-1)/gamma)-1)/eta_c)+1

Tt[3] = Tt[2] * tau_c
Pt[3] = Pt[2] * pi_c

T[3] = (Cp*Tt[3] - u[3]**2/2)/Cp
M[3] = u[3]/np.sqrt(gamma*R*T[3])
P[3] = Pt[3]/(1+(gamma-1)/2 * M[3]**2)**(gamma/(gamma-1))

rho[3] = P[3]/(R*T[3])

A[3] = m_dot / (rho[3] * u[3])

# %% [markdown]
# ## turbine face

# %%
M[4] = 1
Tt[4] = given.Tt4

pi_b = given.pi_b

Pt[4] = Pt[3] * pi_b

T[4] = Tt[4]/(1+(gamma-1)/2 * M[4]**2)
P[4] = Pt[4]/(1+(gamma-1)/2 * M[4]**2)**(gamma/(gamma-1))

rho[4] = P[4]/(R*T[4])

u[4] = M[4]*np.sqrt(gamma*R*T[4])

A[4] = m_dot / (rho[4] * u[4])

tau_lambda = Tt[4] / T[0]

m_dot_f = m_dot * Cp * (Tt[4]-Tt[3]) / h


# %% [markdown]
# ## turbine exit

# %%
u[5] = 0.8 * u[4]

eta_t = given.eta_t

tau_t = 1-(tau_c - 1)*(tau_r/tau_lambda)
pi_t = (1-1/eta_t*(1-tau_t))**(gamma/(gamma-1))

Tt[5] = Tt[4] * tau_t
Pt[5] = Pt[4] * pi_t

T[5] = (Cp*Tt[5] - u[5]**2/2)/Cp
M[5] = u[5]/np.sqrt(gamma*R*T[5])
P[5] = Pt[5]/(1+(gamma-1)/2 * M[5]**2)**(gamma/(gamma-1))

rho[5] = P[5]/(R*T[5])

A[5] = m_dot / (rho[5] * u[5])

# %% [markdown]
# ## nozzle

# %%
M[8] = 1

tau_n = 1
pi_n = 1

Tt[8] = Tt[5] * tau_n
Pt[8] = Pt[5] * pi_n

T[8] = Tt[8]/(1+(gamma-1)/2 * M[8]**2)
P[8] = Pt[8]/(1+(gamma-1)/2 * M[8]**2)**(gamma/(gamma-1))

rho[8] = P[8]/(R*T[8])

u[8] = M[8]*np.sqrt(gamma*R*T[8])

A[8] = m_dot / (rho[8] * u[8])



# %%
# r.T8 = r.Tt8 / (1+(gamma-1)/2*r.M2**2)
# r.P8 = r.Pt8 / (1+(gamma-1)/2*r.M2**2)**(gamma/(gamma-1))


# r.u8 = r.M2 * np.sqrt(gamma*R*r.T8)
# r.rho8 = r.P8 / (r.T8*R)
# r.A8 = r.m_dot_air / (r.rho8 * r.u8)
# r.F = (r.m_dot_air) * (r.u8 - r.u0) + (r.P8 - r.P0) * r.A8

# %% [markdown]
# # Table

# %%
table = np.array([
    u,
    P,
    T,
    M,
    rho,
    A,
    Pt,
    Tt 
])

table = np.delete(table,(6,7),1)
print(table)

# %%


# %%


