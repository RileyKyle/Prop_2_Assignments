# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:32:26 2024

@author: Riley Kyle
"""

import compressible_aero_calc as aero
import numpy as np
import matplotlib.pyplot as plt
#%% vehicle info

width = 1 # m
mass = 10000.0 # kg

class inlet:
    # will see an oblique shock going from the tip to the cowl
    delta = np.deg2rad(6.5)
    xl = 5 # m
    Ai = 0.91506138 # m^2
    A = xl * width / np.cos(delta)


class comb:
    # constant area heat addition
    xl = 0.5 # m
    class i:
        h = 0.1506297 # m
        A = h * width # m^2
        
    class e:
        h = 0.1506297 # m
        A = h * width # m^2


class cowl:
    # oblique shock
    xl = 1.2599182 # m
    A = xl * width

class noz:
    # isentropic
    xl = 3 # m
    class i:
        h = comb.e.h # m
        A = h * width # m^2
    
    class e:
        h = 1.14383 # m
        A = h * width # m^2
        

class vehicle:
    class top:
        delta = np.deg2rad(2.8524466)
        xl = inlet.xl + comb.xl + noz.xl # m
        A = xl * width / np.cos(delta)
    
    class bottom:
        delta = np.deg2rad(0)
        xl = cowl.xl + comb.xl + noz.xl # m
        A = xl * width
    
class wing:
   # delta = np.deg2rad(7)
   # A = 12 # m^2
    class top:
        # expansion fan
        delta = np.deg2rad(7)
        A = 12 # m^2
        
    class bottom:
        # oblique shock
        delta = np.deg2rad(7)
        A = 12 # m^2
        
#%% freestream conditions
gamma = 1.4 
R = 287 # J/kg-K
Cp = 1004.5 # J/kg-K


alpha = np.deg2rad(3)
V_inf = 3048.8 # m/s
T_inf = 231.24 # K
P_inf = 1185.5 # Pa
rho_inf = P_inf/(R*T_inf)
g = 9.81 # m/s^2

a_inf = np.sqrt(gamma*R*T_inf)
M_inf = V_inf / a_inf

m_dot = inlet.Ai * V_inf * rho_inf

[P_P0_inf, T_T0_inf, rho_rho0_inf] = aero.isentropic_relations(M_inf)

Tt_inf = T_inf/T_T0_inf 
Pt_inf = P_inf/P_P0_inf 

#%% region downstream of incident inlet shock
inlet.Tt = Tt_inf
inlet.theta = inlet.delta + alpha

[M2, P2_P1, T2_T1, rho2_rho1, P02_P01, beta] = aero.oblique_shock_relations(M_inf, inlet.theta) 

inlet.M = M2
inlet.P = P2_P1 * P_inf
inlet.T = T2_T1 * T_inf
inlet.rho = rho2_rho1 * rho_inf
inlet.Pt = P02_P01 * Pt_inf
inlet.beta = beta

inlet.a = np.sqrt(gamma*R*inlet.T)
inlet.V = inlet.M * inlet.a

# find force magnitude and direction
inlet.n_hat = np.array([np.sin(inlet.theta),np.cos(inlet.theta)])
inlet.F_mag = inlet.P * inlet.A
inlet.F = inlet.F_mag * (inlet.n_hat)
#%% reflected shock
cowl.Tt = Tt_inf
cowl.theta = inlet.delta

[M2, P2_P1, T2_T1, rho2_rho1, P02_P01, beta] = aero.oblique_shock_relations(inlet.M, cowl.theta) 

cowl.M = M2
cowl.P = P2_P1 * inlet.P
cowl.T = T2_T1 * inlet.T
cowl.rho = rho2_rho1 * inlet.rho
cowl.Pt = P02_P01 * inlet.Pt
cowl.beta = beta

cowl.a = np.sqrt(gamma*R*cowl.T)
cowl.V = cowl.M * cowl.a

# find force magnitude and direction
cowl.n_hat = np.array([-np.sin(alpha),-np.cos(alpha)])
cowl.F_mag = cowl.P * cowl.A
cowl.F = cowl.F_mag * (cowl.n_hat)
#%% top surface of vehicle

vehicle.top.Tt = Tt_inf
vehicle.top.Pt = Pt_inf
vehicle.top.theta = alpha - vehicle.top.delta 

vehicle.top.M = aero.PM_mach(M_inf,vehicle.top.theta)

[P_P0, T_T0, rho_rho0] = aero.isentropic_relations(vehicle.top.M)

vehicle.top.P = Pt_inf * P_P0  
vehicle.top.T = Tt_inf * T_T0
vehicle.top.rho = vehicle.top.P/(R*vehicle.top.T)

vehicle.top.a = np.sqrt(gamma*R*vehicle.top.T)
vehicle.top.V = vehicle.top.M * vehicle.top.a

# find force magnitude and direction
vehicle.top.n_hat = np.array([-np.sin(vehicle.top.theta),-np.cos(vehicle.top.theta)])
vehicle.top.F_mag = vehicle.top.P * vehicle.top.A
vehicle.top.F = vehicle.top.F_mag * (vehicle.top.n_hat)
#%% bottom of vehicle

vehicle.bottom.Tt = Tt_inf
vehicle.bottom.Pt = Pt_inf
vehicle.bottom.theta = alpha

[M2, P2_P1, T2_T1, rho2_rho1, P02_P01, beta] = aero.oblique_shock_relations(M_inf, vehicle.bottom.theta) 

vehicle.bottom.M = M2
vehicle.bottom.P = P2_P1 * P_inf
vehicle.bottom.T = T2_T1 * T_inf
vehicle.bottom.rho = rho2_rho1 * rho_inf
vehicle.bottom.Pt = P02_P01 * Pt_inf
vehicle.bottom.beta = beta


vehicle.bottom.a = np.sqrt(gamma*R*vehicle.bottom.T)
vehicle.bottom.V = vehicle.bottom.M * vehicle.bottom.a

# find force magnitude and direction
vehicle.bottom.n_hat = np.array([np.sin(alpha),np.cos(alpha)])
vehicle.bottom.F_mag = vehicle.bottom.P * vehicle.bottom.A
vehicle.bottom.F = vehicle.bottom.F_mag * (vehicle.bottom.n_hat)


#%% top of wing

wing.top.Tt = Tt_inf
wing.top.Pt = Pt_inf
wing.top.theta = alpha + wing.bottom.delta 

wing.top.M = aero.PM_mach(M_inf,wing.top.theta)

[P_P0, T_T0, rho_rho0] = aero.isentropic_relations(wing.top.M)

wing.top.P = Pt_inf * P_P0  
wing.top.T = Tt_inf * T_T0
wing.top.rho = wing.top.P/(R*wing.top.T)

wing.top.a = np.sqrt(gamma*R*wing.top.T)
wing.top.V = wing.top.M * wing.top.a

# find force magnitude and direction
wing.top.n_hat = np.array([-np.sin(wing.top.theta),-np.cos(wing.top.theta)])
wing.top.F_mag = wing.top.P * wing.top.A
wing.top.F = wing.top.F_mag * (wing.top.n_hat)

#%% bottom of wing

wing.bottom.Tt = Tt_inf
wing.bottom.theta = alpha + wing.bottom.delta 

[M2, P2_P1, T2_T1, rho2_rho1, P02_P01, beta] = aero.oblique_shock_relations(M_inf, wing.bottom.theta) 

wing.bottom.M = M2
wing.bottom.P = P2_P1 * P_inf
wing.bottom.T = T2_T1 * T_inf
wing.bottom.rho = rho2_rho1 * rho_inf
wing.bottom.Pt = P02_P01 * Pt_inf
wing.bottom.beta = beta


wing.bottom.a = np.sqrt(gamma*R*wing.bottom.T)
wing.bottom.V = wing.bottom.M * wing.bottom.a

# find force magnitude and direction
wing.bottom.n_hat = np.array([np.sin(wing.bottom.theta),np.cos(wing.bottom.theta)])
wing.bottom.F_mag = wing.bottom.P * wing.bottom.A
wing.bottom.F = wing.bottom.F_mag * (wing.bottom.n_hat)

#%% combustor/nozzle
Fwing = wing.bottom.F + wing.top.F
Finlet = inlet.F + cowl.F
Ftop = vehicle.top.F
Fbottom = vehicle.bottom.F
F_grav = np.array([0,-mass*g])
F_aero = Fwing + Finlet + Ftop + Fbottom 

# the inlet to the combustor has the same properties as that behind the second shock
comb.i.M = cowl.M
comb.i.P = cowl.P
comb.i.T = cowl.T
comb.i.rho = cowl.rho
comb.i.Pt = cowl.Pt
comb.i.Tt = cowl.Tt
comb.i.a = cowl.a
comb.i.V = cowl.V

# direction vectors for the inlet and outlet of the engine
comb.i.n_hat = np.array([-np.cos(alpha), np.sin(alpha)])
noz.e.n_hat = np.array([  np.cos(alpha),-np.sin(alpha)])

# rayleigh flow equations
f_ray = lambda M: ((1+(gamma-1)/2*M**2)/(1+gamma*M**2)**2)*M**2
M_ray_super = lambda f: np.sqrt(2*f/(1-2*gamma*f-((1-2*(gamma+1)*f)**(1/2)))) 

# nozzle area relations
G_noz = lambda M: 1/M*(2/(gamma+1)*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))


# start = 102000000
# stop = 102500000

# Q_dot_list = np.arange(start, stop, 1)

Q_required = 102457683 * .5
Q_dot_list  = np.array([Q_required])

F_axial = np.zeros(Q_dot_list.size)

# I originally used this to solve for the Q_dot value iteritively
# After finding Q_dot I just removed that iteritive part and plugged in the solution
for i in range(Q_dot_list.size):
    Q_dot = Q_dot_list[i]
    q = Q_dot/m_dot
    f1 = f_ray(comb.i.M)
    f2 = (1+q/(Cp*comb.i.Tt))*f1
    comb.e.M = M_ray_super(f2)
    
    [P_P0, T_T0, rho_rho0] = aero.isentropic_relations(comb.e.M)
    
    comb.e.P = comb.i.P * (1+gamma*comb.i.M**2)/(1+gamma*comb.e.M**2)
    comb.e.Pt = comb.e.P / P_P0
    comb.e.Tt = comb.i.Tt * f2/f1
    comb.e.T = comb.e.Tt * T_T0
    comb.e.a = np.sqrt(gamma*R*comb.e.T)
    comb.e.V = comb.e.M * comb.e.a
    
    G2G1 = noz.e.A / noz.i.A
    G1 = G_noz(comb.e.M)
    G2 = G2G1 * G1
    
    # used secant method to solve for the Mach number at the exit
    noz.e.M = aero._secant_method(G_noz,G2, 11, 6)
    
    [P_P0, T_T0, rho_rho0] = aero.isentropic_relations(noz.e.M)
    
    noz.e.Pt = comb.e.Pt
    noz.e.P = noz.e.Pt * P_P0
    noz.e.Tt = comb.e.Tt
    noz.e.T = noz.e.Tt * T_T0
    
    noz.e.a = np.sqrt(gamma*R*noz.e.T)
    noz.e.V = noz.e.M * noz.e.a
    
    # quasi 1D flow from the combustor enterance to the nozzle exit
    Fx_prime = -((m_dot*noz.e.V+noz.e.P*noz.e.A)-(m_dot*comb.i.V+comb.i.P*comb.i.A)) 
    
    # Multiply the qauasi 1D force by the direction vector
    F_thrust = Fx_prime * (noz.e.n_hat)
    
    F_total = F_aero + F_thrust + F_grav
    
    # if F_total[0] < 1e-4:
    #     Q_required = Q_dot_list[i]
    #     break
    
    # F_axial[i] = F_total[0] 

#%% forces



print(f'inlet (N):          {Finlet}')
print(f'wing top (N):       {wing.top.F}')
print(f'wing bottom (N):    {wing.bottom.F}')
print(f'wing total (N):     {Fwing}')
print(f'bottom surface (N): {Fbottom}')
print(f'top surface (N):    {Ftop}\n')
print(f'aero surfaces (N):  {F_aero}')
print(f'thrust (N):         {F_thrust}')
print(f'total (N):          {F_total}\n\n')

print(f'Q_dot (MW):         {Q_required*10**-6}')
print(f'm_dot_f (kg/s):     {Q_required/(1.2*10**8)}')


#%% print outputs
print('free-stream conditions')
print('----------------------')
print(f'Mach number:                {M_inf}')
print(f'static pressure (Pa):       {P_inf}')
print(f'static temperature (K):     {T_inf}')
print(f'velocity magnitude (m/s):   {V_inf}')
print(f'u (m/s):                    {V_inf}')
print(f'v (m/s):                    {0}')
print(f'total temperature (K):      {Tt_inf}')
print(f'total pressure (kPa):       {Pt_inf/1000}')
print('\n\n')

print('inlet conditions')
print('----------------------')
print(f'beta (deg):                 {np.rad2deg(inlet.beta)}')
print(f'theta (deg):                {np.rad2deg(inlet.theta)}')
print(f'Mach number:                {inlet.M}')
print(f'static pressure (Pa):       {inlet.P}')
print(f'static temperature (K):     {inlet.T}')
print(f'velocity magnitude (m/s):   {inlet.V}')
print(f'u (m/s):                    {inlet.V*np.cos(inlet.theta)}')
print(f'v (m/s):                    {-inlet.V*np.sin(inlet.theta)}')
print(f'total temperature (K):      {inlet.Tt}')
print(f'total pressure (kPa):       {inlet.Pt/1000}')
print('\n\n')

print('cowl conditions')
print('----------------------')
print(f'beta (deg):                 {np.rad2deg(cowl.beta)}')
print(f'theta (deg):                {np.rad2deg(cowl.theta)}')
print(f'Mach number:                {cowl.M}')
print(f'static pressure (Pa):       {cowl.P}')
print(f'static temperature (K):     {cowl.T}')
print(f'velocity magnitude (m/s):   {cowl.V}')
print(f'u (m/s):                    {cowl.V*np.cos(alpha)}')
print(f'v (m/s):                    {-cowl.V*np.sin(alpha)}')
print(f'total temperature (K):      {cowl.Tt}')
print(f'total pressure (kPa):       {cowl.Pt/1000}')
print('\n\n')

print('combustor exit')
print('----------------------')
print(f'Mach number:                {comb.e.M}')
print(f'static pressure (Pa):       {comb.e.P}')
print(f'static temperature (K):     {comb.e.T}')
print(f'velocity magnitude (m/s):   {comb.e.V}')
print(f'u (m/s):                    {comb.e.V*np.cos(alpha)}')
print(f'v (m/s):                    {comb.e.V*np.sin(-alpha)}')
print(f'total temperature (K):      {comb.e.Tt}')
print(f'total pressure (kPa):       {comb.e.Pt/1000}')
print('\n\n')

print('nozzle exit')
print('----------------------')
print(f'Mach number:                {noz.e.M}')
print(f'static pressure (Pa):       {noz.e.P}')
print(f'static temperature (K):     {noz.e.T}')
print(f'velocity magnitude (m/s):   {noz.e.V}')
print(f'u (m/s):                    {noz.e.V*np.cos(alpha)}')
print(f'v (m/s):                    {noz.e.V*np.sin(-alpha)}')
print(f'total temperature (K):      {noz.e.Tt}')
print(f'total pressure (kPa):       {noz.e.Pt/1000}')
print('\n\n')

print('top surface conditions')
print('----------------------')
print(f'theta (deg):                {np.rad2deg(vehicle.top.theta)}')
print(f'Mach number:                {vehicle.top.M}')
print(f'static pressure (Pa):       {vehicle.top.P}')
print(f'static temperature (K):     {vehicle.top.T}')
print(f'velocity magnitude (m/s):   {vehicle.top.V}')
print(f'u (m/s):                    {vehicle.top.V*np.cos(vehicle.top.theta)}')
print(f'v (m/s):                    {-vehicle.top.V*np.sin(vehicle.top.theta)}')
print(f'total temperature (K):      {vehicle.top.Tt}')
print(f'total pressure (kPa):       {vehicle.top.Pt/1000}')
print('\n\n')

print('bottom surface conditions')
print('----------------------')
print(f'beta (deg):                 {np.rad2deg(vehicle.bottom.beta)}')
print(f'theta (deg):                {np.rad2deg(vehicle.bottom.theta)}')
print(f'Mach number:                {vehicle.bottom.M}')
print(f'static pressure (Pa):       {vehicle.bottom.P}')
print(f'static temperature (K):     {vehicle.bottom.T}')
print(f'velocity magnitude (m/s):   {vehicle.bottom.V}')
print(f'u (m/s):                    {vehicle.bottom.V*np.cos(vehicle.bottom.theta)}')
print(f'v (m/s):                    {-vehicle.bottom.V*np.sin(vehicle.bottom.theta)}')
print(f'total temperature (K):      {vehicle.bottom.Tt}')
print(f'total pressure (kPa):       {vehicle.bottom.Pt/1000}')
print('\n\n')

print('wing top conditions')
print('----------------------')
print(f'theta (deg):                {np.rad2deg(wing.top.theta)}')
print(f'Mach number:                {wing.top.M}')
print(f'static pressure (Pa):       {wing.top.P}')
print(f'static temperature (K):     {wing.top.T}')
print(f'velocity magnitude (m/s):   {wing.top.V}')
print(f'u (m/s):                    {wing.top.V*np.cos(wing.top.theta)}')
print(f'v (m/s):                    {-vehicle.top.V*np.sin(-wing.top.theta)}')
print(f'total temperature (K):      {wing.top.Tt}')
print(f'total pressure (kPa):       {wing.top.Pt/1000}')
print('\n\n')

print('wing bottom conditions')
print('----------------------')
print(f'beta (deg):                 {np.rad2deg(wing.bottom.beta)}')
print(f'theta (deg):                {np.rad2deg(wing.bottom.theta)}')
print(f'Mach number:                {wing.bottom.M}')
print(f'static pressure (Pa):       {wing.bottom.P}')
print(f'static temperature (K):     {wing.bottom.T}')
print(f'velocity magnitude (m/s):   {wing.bottom.V}')
print(f'u (m/s):                    {wing.bottom.V*np.cos(wing.bottom.theta)}')
print(f'v (m/s):                    {-wing.bottom.V*np.sin(wing.bottom.theta)}')
print(f'total temperature (K):      {wing.bottom.Tt}')
print(f'total pressure (kPa):       {wing.bottom.Pt/1000}')
print('\n\n')


