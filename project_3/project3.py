import numpy as np
import matplotlib.pyplot as plt
import ion_drive_solver as ion_drive

Isp = 3500 # sec
prop_charge = 2500 # C/kg
max_V_gradient = 1.0e+5
drive_power = 0.35 # percent/100

m_0 = 80_000 # Kg
m_L = 5000 # Kg

alpha_prime = 600 # W/kg

# # Trajectory

# %%
g0 = 9.81
C = Isp * g0 # effective velocity

dV_accel = 35_000
dV_deccel = 35_000
dV_tot = dV_accel + dV_deccel



m_p_func = lambda dV, m_0: m_0 * (1-np.exp(-dV/C))

m_p_accel = m_p_func(dV_accel,m_0)
m_p_deccel = m_p_func(dV_deccel,m_0-m_p_accel)
m_p_tot = m_p_accel + m_p_deccel

engine_mass = m_0 - m_L - m_p_tot

# # Ion Drive
epsilon0 = 8.85e-12

W_engine = engine_mass * alpha_prime

W_ion =  W_engine / (1 + drive_power)
W_ionizer= W_ion * drive_power

voltage = (Isp * g0)**2 / (2 * prop_charge)
current = W_ion / voltage
acc_distance_cm = voltage / max_V_gradient
d = acc_distance_cm/100

j_max = 0.4444 * epsilon0 * np.sqrt(2*prop_charge) * voltage**(3/2) / d**2

A_beam = current/j_max
D_beam = np.sqrt(A_beam * 4 / np.pi)
F = 0.8888 * epsilon0 * (voltage/d)**2 * A_beam
m_dot_ion = F/C

# # Brayton cycle
Ru = 8314.4598 # J/mol/K
# nitrogen
molar_mass = 28 # kg/kmol
gamma = 1.4
R = Ru/molar_mass
Cp = gamma*R/(gamma-1)

T3 = 3000 # K
Nsize = 500
pi_c_list = np.linspace(1,9,Nsize)
Ta_list = np.linspace(800,1600, Nsize)

pi_c_list, Ta_list = np.meshgrid(pi_c_list,Ta_list)
pi_c_list = pi_c_list.flatten()
Ta_list = Ta_list.flatten()
eta_c = 0.90
eta_t = 0.90

epsilon = 1.0
sigma = 5.670374419e-8

def brayton_cycle(pi_c, Ta):
    pi_t = 1/pi_c
    T4 = T3 * (1 - eta_t * (1 - pi_t**((gamma-1)/gamma)))
    T2 = Ta * (1 + 1/eta_c * (pi_c**((gamma-1)/gamma) - 1))
    eta = 1 - (T4 - Ta)/(T3 - T2)
    Arad = (1/eta - 1)/(3*(T4 - Ta)*sigma*epsilon) * (1/Ta**3 - 1/T4**3)
    return T2,T4,eta,Arad

T2,T4,eta,Arad = brayton_cycle(pi_c_list,Ta_list)
for i in range(pi_c_list.size):
    if Arad[i] <= 0:
        Arad[i] = np.nan 
Arad2 = Arad[~np.isnan(Arad)]
idx = np.where(Arad == np.min(Arad2))

Ta = Ta_list[idx][0]
pi_c = pi_c_list[idx][0]
T2 = T2[idx][0]
T4 = T4[idx][0]
eta = eta[idx][0]
Arad = Arad[idx][0]


min_rad_area = Arad*W_engine
Q_added = W_engine/eta
m_dot = (1/eta-1)/(T4-Ta)/Cp* W_engine

W_c = m_dot * Cp * (T2-Ta)
W_t = m_dot * Cp * (T3-T4)
net_power = W_t - W_c

# # Rocket Trajectory
TOF1 = m_p_accel / m_dot_ion#/3600/24
TOF2 = m_p_deccel / m_dot_ion#/3600/24

g0 = 9.81
C = Isp * g0 # effective velocity

dV_accel = 35_000
dV_deccel = 35_000
dV_tot = dV_accel + dV_deccel

N_itr = 100_000_00

m_accel_list = np.linspace(m_0,m_0-m_p_accel,N_itr+1)
t_accel_list = np.linspace(0,TOF1,N_itr+1) 


def distance(V0,TOF,m0,mp,N_itr = 10_000_000,accel = True):
    m_list = np.linspace(m0,m0-mp,N_itr+1)
    t_list = np.linspace(0,TOF,N_itr+1) 

    mult = 1
    if accel == False:
        mult = -1
    Vi = V0
    Dist = 0
    for i in range(N_itr):
        dt = t_list[i+1] - t_list[i]
        
        mi = m_list[i]
        dm = m_list[i+1] - mi
        dv = -C * (dm/mi) * mult
        Vi += dv
        Dist += Vi * dt
    print(Dist)


print(f'{TOF1 = }')
print(f'{TOF2 = }')

distance(0,TOF1,m_0,m_p_accel,N_itr=20_000_000)
TOF_cruise = 4000*24*3600
distance(35_000,TOF_cruise,m_0-m_p_accel,0,N_itr=20_000_000)
distance(35_000,TOF2,m_0-m_p_accel,m_p_deccel,N_itr=20_000_000,accel=False)