import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# analytical solution functions
def sod_func(P):
    return (P - pr)*(( ((1 - mu**2)**2)*((rr*(P + mu*mu*pr))**(-1)) )**0.5) - 2*(np.sqrt(gam)/(gam - 1))*(1 - P**((gam - 1)/(2*gam)))
def analytic_sod(t):
    c_l = np.sqrt(gam*pl/rl)
    c_r = np.sqrt(gam*pr/rr)
    P_post = fsolve(sod_func,0.31)
    v_post = 2*(np.sqrt(gam)/(gam - 1))*(1 - P_post**((gam - 1)/(2*gam)))
    rho_post = rr*(( (P_post/pr) + mu**2 )/(1 + mu*mu*(P_post/pr)))
    v_shock = v_post*((rho_post/rr)/( (rho_post/rr) - 1))
    rho_middle = (rl)*(P_post/pl)**(1/gam)

    #Key Positions
    x1 = x0 - c_l*t
    x3 = x0 + v_post*t
    x4 = x0 + v_shock*t
    c_2 = c_l - ((gam - 1)/2)*v_post
    x2 = x0 + (v_post - c_2)*t

    #start setting values
    #boundaries (can be set)
    x_min = 0.00
    x_max = L

    rho = np.zeros(N)
    P = np.zeros(N)
    u = np.zeros(N)
    e = np.zeros(N)

    for index in range(0,N):
        if x[index] < x1:
            rho[index] = rl
            P[index] = pl
            u[index] = ul
        elif (x1 <= x[index] and x[index] <= x2):
            c = mu*mu*((x0 - x[index])/t) + (1 - mu*mu)*c_l
            rho[index] = rl*(c/c_l)**(2/(gam - 1))
            P[index] = pl*(rho[index]/rl)**gam
            u[index] = (1 - mu*mu)*( (-(x0-x[index])/t) + c_l)
        elif (x2 <= x[index] and x[index] <= x3):
            rho[index] = rho_middle
            P[index] = P_post
            u[index] = v_post
        elif (x3 <= x[index] and x[index] <= x4):
            rho[index] = rho_post
            P[index] = P_post
            u[index] = v_post
        elif x4 < x[index]:
            rho[index] = rr
            P[index] = pr
            u[index] = ur
        e[index] = P[index]/((gam - 1)*rho[index])
        #e[index] = ener(P[index],rho[index],u[index])
    return rho,P,u,e

# load parameters
T,CFL,x0,rl,pl,u0,rr,pr,delx,gam,L = np.loadtxt('calibrationShock.dat')
ul = u0
ur = u0
mu = np.sqrt( (gam-1)/(gam+1) )
N=int(L/delx)

# load numeric solution
data = np.loadtxt('dataShock.dat')
x = data[:,0]
r = data[:,1]
v = data[:,2]
p = data[:,3]
e = data[:,4]

# get the analytical solution
ra,pa,va,ea = analytic_sod(T)

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, r,'ro-',linewidth=1.5,label='Numerical Sol')
plt.plot(x, ra,'b-',linewidth=2.0,label='Analytical Sol')
plt.ylim(np.min(r)-0.1, np.max(r)+0.1)
plt.xlabel(r"Position ($x$)")
plt.ylabel(r"$\rho$")
plt.title(r"Density - Shock Problem")
plt.axvline(0.5, color='k', linestyle='solid')
plt.legend()
plt.savefig('DensidadShock.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, p,'ro-',linewidth=1.5,label='Numerical Sol')
plt.plot(x, pa,'b-',linewidth=2.0,label='Analytical Sol')
plt.ylim((np.min(p)-0.1, np.max(p)+0.1))
plt.xlabel(r"Position ($x$)")
plt.ylabel(r"$P$")
plt.title(r"Pressure - Shock Problem")
plt.axvline(0.5, color='k', linestyle='solid')
plt.legend()
plt.savefig('PresionShock.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, v,'ro-',linewidth=1.5,label='Numerical Sol')
plt.plot(x, va,'b-',linewidth=2.0,label='Analytical Sol')
plt.ylim((np.min(v)-0.1, np.max(v)+0.1))
plt.xlabel(r"Position ($x$)")
plt.ylabel(r"$v$")
plt.title(r"Velocity - Shock Problem")
plt.axvline(0.5, color='k', linestyle='solid')
plt.legend()
plt.savefig('VelocidadShock.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, e,'ro-',linewidth=1.5,label='Numerical Sol')
plt.plot(x, ea,'b-',linewidth=2.0,label='Analytical Sol')
plt.ylim((np.min(e)-0.1, np.max(e)+0.1))
plt.xlabel(r"Position ($x$)")
plt.ylabel(r"$e$")
plt.title(r"Internal Energy - Shock Problem")
plt.axvline(0.5, color='k', linestyle='solid')
plt.legend()
plt.savefig('EnergiaShock.pdf', format='pdf')
plt.close()
