import numpy as  np

class air:
    R = 287
    Cp = 1004.5
    gamma = 1.4

def quasi_1D_flow(
        inlet: list[float], 
        Af = lambda x: 1 , 
        L: float = 1 , 
        gas: str = 'air',
        N: int   = 100_000, 
        w: float = 0, 
        q: float = 0,
        Cf: float = 0,
        eta: float = 1):
    
    gas = eval(f'{gas}')
    Cp, R, gamma = gas.Cp, gas.R, gas.gamma

    u, T, P, M = np.empty(N+1), np.empty(N+1), np.empty(N+1), np.empty(N+1)

    u[0] =     inlet[0]
    T[0] =     inlet[1] 
    P[0] =     inlet[2]
    M[0] = u[0]/np.sqrt(gamma*R*T[0])
    print(M[0])
    rho =  P[0]/(R*T[0])


    x = 0
    dx = L / N
    A = Af(x)

    dw = w/N
    dq = q/N
    m_dot = rho*u[0]*A
    

    dS = 1
    dW_dot = dw/m_dot

    for j in range(N+1):
        
        x += dx

        tau_w = 1/2*rho*u[j]**2*Cf

        if tau_w > 0:
            r = np.sqrt(A/np.pi)
            c = 2*np.pi*r
            dS = c * dx
            
        Ap1 = Af(x)
        dA = Ap1 - A
        
        du = -(
            ((1/P[j] * (rho*eta*dW_dot - tau_w*dS/A) - 1/(Cp*T[j]) * (dq + dw)) + dA/A)/
            (1/u[j] - u[j]/(R*T[j]) + u[j]/(Cp*T[j]))
        )

        dP = rho*eta*dW_dot - tau_w*dS/A - rho*u[j]*du

        dT = 1/(Cp*T[j]) * (dq + dw - u[j]*du)

        drho = rho*(dP/P[j] - dT/T[j])

        A = Ap1
        
        u[j+1]    = u[j]   + du
        T[j+1]    = T[j]   + dT
        P[j+1]    = P[j]   + dP
        rho  = rho + drho

        M[j+1] = u[j+1]/np.sqrt(gamma*R*T[j+1])

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