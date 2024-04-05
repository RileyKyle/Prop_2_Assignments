import numpy as  np

class air:
    R = 287
    Cp = 1004.5
    gamma = 1.4


def Friction_comp_circle(rho, u, A, dx, Cf):
    if Cf == None or Cf == 0:
        return 0
    else:
        r = np.sqrt(A/np.pi)
        c = 2*np.pi*r
        dS = c * dx
        tau_w = 1/2*rho*u**2*Cf
        return tau_w*dS/A/rho
    
def dq_conv_circle(Tw,Tt,dx,A,Cp,Cf):
    if Cf == None or Cf == 0  or Tw == None:
        return 0
    else:
        D = 2 * np.sqrt(A/np.pi)
        return 2*Cp*Cf*(Tw-Tt)*(dx/D)
    

def quasi_1D_flow(
        inlet: list[float], 
        Af = lambda x: 1 , 
        L: float = 1.0 , 
        gas: str = 'air',
        N: int   = 100_000, 
        w: float = 0, 
        q: float = 0,
        Cf: float = None,
        eta: float = 1,
        Tw: float = None
        ):
    
    # load gas
    gas = eval(f'{gas}')
    Cp, R, gamma = gas.Cp, gas.R, gas.gamma

    # initialize arrays
    u, T, P, M = np.empty(N+1), np.empty(N+1), np.empty(N+1), np.empty(N+1)

    # add inlet conditions to the arrays
    u[0], T[0], P[0] = inlet[0], inlet[1], inlet[2] 
    M[0] = u[0]/np.sqrt(gamma*R*T[0])
    rho =  P[0]/(R*T[0])

    # this assumes that shaft work and radiative heat transfer are spread evenly through the tube
    dw = w/N
    dq_rad = q/N
    
    # break flow tube into finite chunks
    xl = np.linspace(0,L,N+1)
    dx = float(L / N)
    A = Af(0)
    
    for j in range(N):
        
        x = xl[j+1] # I found that this ran the same speed as if I just added dx
                    # this is good because just adding dx can introduce some floating point error
        Ap1 = Af(x)
        dA = Ap1 - A

        # skin friction component
        SF = Friction_comp_circle(rho, u[j], A, dx, Cf)
        
        # convective and radiative heating
        Tt = T[j]*(1+(gamma-1)/2*M[j]**2)
        dq = dq_rad + dq_conv_circle(Tw,Tt,dx,A,Cp,Cf)
            
        # change in velocity
        du = -(
            (rho/P[j]*(eta*dw -SF) - (dq + dw)/(Cp*T[j]) + dA/A) /
            (1/u[j] - u[j]*rho/(P[j]) + u[j]/(Cp*T[j]))
        )
        
        # change in pressure
        dP = rho*(eta*dw - SF - u[j]*du)

        # change in temperature
        dT = (dq + dw - u[j]*du)/Cp

        # change in density
        drho = rho*(dP/P[j] - dT/T[j])

        # add the change
        u[j+1] = u[j] + du
        T[j+1] = T[j] + dT
        P[j+1] = P[j] + dP
        rho    = rho  + drho

        M[j+1] = u[j+1]/np.sqrt(gamma*R*T[j+1])
        
        A = Ap1
        
        # error checking statements
        if u[j+1] < 0 or T[j+1] < 0 or P[j+1] < 0:
            print('error 1')
            break

        if M[0] < 1 and M[j+1] > 1:
            print('error 2') 
            break

        if M[0] > 1 and M[j+1] < 1:
            print('error 3')
            break 
        
    return [u,T,P], M, x, A