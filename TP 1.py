
# -*- coding: utf-8 -*-
"""
CONTROL OF AIRCRAFT - PRACTICAL WORK
Author:
    LAUPIES Raphaêl
    CHARDON DU RANQUET Quentin
    Choosed number : 82
"""
from __future__ import unicode_literals

import numpy as np
from control import StateSpace
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy import interpolate
from control import matlab
import control
from pylab import * 
from atm_std import *
from sisopy31 import *
#%% 1)
# Constants

g = 9.81
m = 8400        # kg
S = 34          # m^2
alt_ft = 21230  # ft (Subject 82)
alt_m = alt_ft * 0.3048
M = 0.95        # Mach (Subject 82)
R = 287.05      # J/(kg.K)
l_ref = 5.24
rg = 2.65
c = 0.52 

# Graphical data

Cx0 = 0.018       # [M=0.95]
k = 0.28
Cz_alpha = 2.85   # 1/rad
Cz_delta_m = 0.75 # 1/rad
delta_m0 = 0.002  # rad
alpha0 = 0.013    # rad
f = 0.545
f_delta = 0.82
Cmq = -0.55       # s/rad


# Computed parameters

# Values will be computed from US Standard Atmosphere 1976 model
# For alt = 21230 ft = 6471 m, using interpolation from atm_std:
hgeo, rho, a = get_cte_atm(alt_m)
T = 249.2 - (alt_m - 6000) * (249.2 - 236.2) / 2000
P = rho * R * T
V_eq = M * a                 # m/s
Q = 0.5 * rho * V_eq**2      # Pa
l_t = (3/2) * l_ref
XF = -f * l_t
XG = -c * l_t
X = XF - XG
XF_delta = -f_delta * l_t
Y = XF_delta - XG


# Equilibrium point algorithm
eps=1e-6
alpha_eq = 0.0
Fpx_eq = 0.0

diff = 100000000  # initial value to enter in the while
iter_count = 0

while diff > eps:

    Czeq = (m*g - Fpx_eq * np.sin(alpha_eq)) / (Q * S)

    Cxeq = Cx0 + k * Czeq**2

    Cx_delta_m = 2 * k * Czeq * Cz_delta_m

    num = Cxeq * np.sin(alpha_eq) + Czeq * np.cos(alpha_eq)
    den = Cx_delta_m * np.sin(alpha_eq) + Cz_delta_m * np.cos(alpha_eq)
    delta_meq = delta_m0 - (num / den) * (X / (Y - X))

    alpha_next = alpha0 + (Czeq / Cz_alpha) - (Cz_delta_m / Cz_alpha) * delta_meq

    Fpx_next = Q * S * Cxeq / np.cos(alpha_next)

    diff = abs(alpha_next - alpha_eq)

    alpha_eq = alpha_next
    Fpx_eq = Fpx_next
    iter_count += 1
    
print("number iteration: ", iter_count)
print("alpha_eq =", alpha_eq * 180/np.pi)
print("delta_meq =", delta_meq)
print("Fpx_eq =", Fpx_eq)


#%% 2) State space representation
# Definitions

gamma_eq = 0.0
F_tau = 0.0

F_eq = Q * S * Cxeq / np.cos(alpha_eq)


# Aerodynamic derivatives
Cx_alpha = 2 * k * Czeq * Cz_alpha
CN_alpha = Cx_alpha * np.sin(alpha_eq) + Cz_alpha * np.cos(alpha_eq)
CN_delta_m = Cx_delta_m * np.sin(alpha_eq) + Cz_delta_m * np.cos(alpha_eq)

Cm_alpha = (X / l_ref) * CN_alpha
Cm_delta_m = (Y / l_ref) * CN_delta_m

Iyy = m * rg**2 


# X coefficients
X_V = (2 * Q * S * Cxeq) / (m * V_eq)
X_alpha = (F_eq * np.sin(alpha_eq)) / (m * V_eq) + (Q * S * Cx_alpha) / (m * V_eq)
X_gamma = (g * np.cos(gamma_eq)) / V_eq
X_delta_m = (Q * S * Cx_delta_m) / (m * V_eq)
X_tau = -(F_tau * np.cos(alpha_eq)) / (m * V_eq)

# m coefficients
m_V = 0
m_alpha = (Q * S * l_ref * Cm_alpha) / Iyy
m_q = (Q * S * l_ref**2 * Cmq) / (V_eq * Iyy)
m_delta_m = (Q * S * l_ref * Cm_delta_m) / Iyy


# Z coefficients
Z_V = (2 * Q * S * Czeq) / (m * V_eq)
Z_alpha = (F_eq * np.cos(alpha_eq) / (m * V_eq)) + (Q * S * Cz_alpha) / (m * V_eq)
Z_gamma = (g * np.sin(gamma_eq)) / V_eq
Z_delta_m = (Q * S * Cz_delta_m) / (m * V_eq)
Z_tau = (F_tau * np.sin(alpha_eq)) / (m * V_eq)

# State matrices
A = np.array([
    [-X_V,      -X_gamma,   -X_alpha,   0, 0, 0],
    [ Z_V,       0,          Z_alpha,    0, 0, 0],
    [-Z_V,       0,         -Z_alpha,    1, 0, 0],
    [ 0,         0,          m_alpha,    m_q, 0, 0],
    [ 0,         0,          0,          1,   0, 0],
    [ 0,         V_eq,       0,          0,   0, 0]
])

B = np.array([
    [0],
    [Z_delta_m],
    [-Z_delta_m],
    [m_delta_m],
    [0],
    [0]
])

C = np.eye(6)
D = np.zeros((6,1))


sys = StateSpace(A, B, C, D)

#%% 3)

import control.matlab as cmat
wn, zeta, pole = cmat.damp(sys)


#%% 4)
#%% reduced model 
import math
def mdamp(A):
    roots = np.linalg.eigvals(A)
    
    ri = []
    a = []
    b = []
    w = []
    xi = []
    st = []
    
    for i in range(roots.size):
        ri.append(roots[i])
        
        a.append(roots[i].real)
        b.append(roots[i].imag)
        
        w.append(math.sqrt(a[i]**2 + b[i]**2))
        xi.append(-a[i] / w[i])
        
        signb = '+' if b[i] > 0 else '-'
        
        st.append(
            f"{a[i]:.5f}{signb}j{abs(b[i]):.5f}  "
            f"xi = {xi[i]:.5f}  w = {w[i]:.5f} rad/s"
        )
    
    print(st)
    
    
Ar = np.array([
    [-X_V,     -X_gamma,   -X_alpha,   0],
    [ Z_V,      0,          Z_alpha,    0],
    [-Z_V,      0,         -Z_alpha,    1],
    [ 0,        0,          m_alpha,    m_q]
])
# mdamp(Ar)


Br = np.array([
    [0],
    [Z_delta_m],
    [-Z_delta_m],
    [m_delta_m]
])
# mdamp(Br)

Cr = np.eye(4)
Dr = np.zeros((4,1))

sys_r = StateSpace(Ar, Br, Cr, Dr)
wn_r, zeta_r, pole_r = cmat.damp(sys_r)

#%% 
 

Ai = Ar[2:4, 2:4]
Bi = Br[2:4, 0:1]

mdamp(Ai)

Cia = np.matrix([[1, 0]])
Ciq = np.matrix([[0, 1]])
Di = np.matrix([[0]])

TaDm_ss = control.ss(Ai, Bi, Cia, Di)
print("Transfer function alpha / delta m =")
TaDm_tf = control.tf(TaDm_ss)
print(TaDm_tf)

print("Static gain of alpha / delta m = %f" % (control.dcgain(TaDm_tf)))

TqDm_ss = control.ss(Ai, Bi, Ciq, Di)
print("Transfer function q / delta m =")
TqDm_tf = control.ss2tf(TqDm_ss)
print(TqDm_tf)

print("Static gain of q / delta m = %f" % (dcgain(TqDm_tf)))

figure(1)

Ya, Ta = control.matlab.step(TaDm_tf, arange(0, 10, 0.01))
Yq, Tq = control.matlab.step(TqDm_tf, arange(0, 10, 0.01))

plot(Ta, Ya, 'b', Tq, Yq, 'r', lw=2)

plot([0, Ta[-1]], [Ya[-1], Ya[-1]], 'k--', lw=1)
plot([0, Ta[-1]], [1.05 * Ya[-1], 1.05 * Ya[-1]], 'k--', lw=1)
plot([0, Ta[-1]], [0.95 * Ya[-1], 0.95 * Ya[-1]], 'k--', lw=1)

plot([0, Tq[-1]], [Yq[-1], Yq[-1]], 'k--', lw=1)
plot([0, Tq[-1]], [1.05 * Yq[-1], 1.05 * Yq[-1]], 'k--', lw=1)
plot([0, Tq[-1]], [0.95 * Yq[-1], 0.95 * Yq[-1]], 'k--', lw=1)

minorticks_on()
grid(visible=True, which='both')


title(r"Step response $\alpha/\delta m$ and $q/\delta m$")
legend((r'$\alpha/\delta m$', r'$q/\delta m$'))
xlabel("Time (s)")
ylabel(r'$\alpha$ (rad) & $q$ (rad/s)')

Osa, Tra, Tsa = step_info(Ta, Ya)
Osq, Trq, Tsq = step_info(Tq, Yq)

yya = interp1d(Ta, Ya)
plot(Tsa, yya(Tsa), 'bs')
text(Tsa, yya(Tsa) - 0.2, Tsa)

yyq = interp1d(Tq, Yq)
plot(Tsq, yyq(Tsq), 'r s')
text(Tsq, yyq(Tsq) - 0.2, Tsq)

print("Alpha Settling time 5%% = %f s" % Tsa)
print("q Settling time 5%% = %f s" % Tsq)

savefig("stepalphaq.pdf")

#%%
# question 1 page 26

# sisotool(TqDm_ss)  # Commented: requires slycot, Kr already calculated

#%% 5) CLOSED LOOP WITH GAIN Kr page 26

# q output selection
Cq = np.matrix([[0, 1]])
xi_desired = 0.65
TqDm_tf = control.ss2tf(control.ss(Ai, Bi, Cq, 0))
# Kr determined via sisotool for xi=0.65
Kr = 0.10

# Closed Loop State Space
# Control law: delta_m = Kr * (q_c - q)
# Note: m_delta_m < 0, signs adjusted for feedback consistency

# A_k = A + Kr * B * C_q (adjusted for negative m_delta_m)
Ak = Ai + Kr * Bi @ Cq
# B_k = Kr*B (input gain for reference command)
Bk = Kr * Bi
# C_k = C_out (output selection - here we want all states as output)
Ck_2x2 = np.eye(2)  # All states as output
# D_k = 0 (no direct feedthrough)
Dk_2x2 = np.zeros((2, 1))

# Create closed-loop state space system (2x2)
sys_cl = control.ss(np.array(Ak), np.array(Bk), Ck_2x2, Dk_2x2)

print("\nClosed-loop system characteristics (2x2):")
wn_cl, zeta_cl, pole_cl = control.matlab.damp(sys_cl)
print(f"Natural frequencies: {wn_cl}")
print(f"Damping ratios: {zeta_cl}")

# Analyze closed-loop poles
print("\nClosed-loop pole analysis:")
print("Poles of the closed-loop system:")
mdamp(np.array(Ak))

# Test the closed-loop system
print("\nTesting closed-loop step response")

# Closed-loop transfer function for q output
Cq_array = np.array([[0, 1]])  # Select second state (q or theta)
Dk_q = np.zeros((1, 1))
sys_cl_q = control.ss(np.array(Ak), np.array(Bk), Cq_array, Dk_q)

# Step response
t_cl = np.arange(0, 10, 0.01)
Y_cl, T_cl = control.matlab.step(sys_cl_q, t_cl)

# Calculate performance metrics
OS_cl, Tr_cl, Ts_cl = step_info(T_cl, Y_cl)

print(f"Closed-loop q response (with Kr = {Kr:.6f}):")
print(f"  Overshoot (OS):        {OS_cl:.3f} %")
print(f"  Rise time (Tr):        {Tr_cl:.3f} s")
print(f"  Settling time (Ts 5%): {Ts_cl:.3f} s")
 
# %% 6) Washout filter and comparison
# Washout filter (tau > 0.5s)
tau = 1
W = control.tf([tau, 0], [tau, 1])

# Comparison: Open-loop vs Closed-loop (washout zeros DC feedback)

# Using state-space for proper closed-loop (avoids transfer function sign issues)
# Build closed-loop from state-space representation
Ak_cl = Ai + Kr * Bi @ Ciq  # Closed-loop A matrix with q feedback (+ for neg m_delta_m)
Bk_cl = Kr * Bi
Ck_alpha = np.array([[1, 0]])  # Output = alpha

sys_cl_alpha = control.ss(np.array(Ak_cl), np.array(Bk_cl), Ck_alpha, 0)

# Open-loop (just the plant from delta_m to alpha, scaled by Kr)
sys_open_alpha = control.ss(Ai, Kr * Bi, Cia, 0)

# Simulate step responses
t = np.arange(0, 10, 0.01)
_, Y_open = control.step_response(sys_open_alpha, t)
_, Y_cl = control.step_response(sys_cl_alpha, t)

# For washout, we simulate with a filtered q feedback
# This is more complex and would require augmented state-space
# For simplicity, we show only open-loop vs closed-loop comparison

# Plot comparison
plt.figure(figsize=(8,5))
plt.plot(t, Y_open, 'k--', lw=2, label='Open-loop (Kr * alpha/delta_m)')
plt.plot(t, Y_cl, 'b-', lw=2, label='Closed-loop (with q feedback)')
plt.axhline(y=Y_cl[-1], color='g', linestyle=':', label=f'Steady-state = {Y_cl[-1]:.3f}')
plt.xlabel('Time (s)')
plt.ylabel(r'$\alpha$ (rad) per unit $q_c$')
plt.title('Comparison: Open-loop vs Closed-loop alpha response')
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()

# %% 7) Gamma feedback loop

# 4x4 Model [V, gamma, alpha, q]
C_gamma4 = np.array([[0., 1., 0., 0.]])
C_q4     = np.array([[0., 0., 0., 1.]])

A_q = Ar + Br @ (Kr * C_q4)
B_q = Br*Kr

# Plant: delta_m_qc -> gamma
sys_gamma_plant = control.ss(A_q, B_q, C_gamma4, 0)

# K selection (damping ~= 0.5)
Kgamma = 5.0

# Closed-loop (gamma_c -> gamma)
A_gamma = A_q + B_q @ (Kgamma * C_gamma4)
B_gamma = B_q * Kgamma
C_gamma = C_gamma4
D_gamma = np.array([[0.]])

# 1) State space
sys_gamma_cl = control.ss(A_gamma, B_gamma, C_gamma, D_gamma)

# 2) Transfer function
T_gamma_cl = control.ss2tf(sys_gamma_cl)

# 3) Poles/Damping
wn_g, zeta_g, poles_g = control.damp(sys_gamma_cl)

print("Closed-loop poles (gamma loop):")
print(poles_g)
print("Damping ratios:")
print(zeta_g)
print("Natural frequencies (rad/s):")
print(wn_g)

# 4) Step response gamma_c -> gamma
t = np.linspace(0, 15, 1000)
_, y_gamma = control.step_response(sys_gamma_cl, t)

plt.figure(figsize=(8,4))
plt.plot(t, y_gamma, lw=2)
plt.grid(True)
plt.xlabel("Time (s)")
plt.ylabel(r"$\gamma$ (rad)")
plt.title(r"Step response of closed-loop $\gamma$ control")
plt.tight_layout()
plt.show()

# (optional check for specs)
OS_g, Tr_g, Ts_g = step_info(t, y_gamma)
print(f"Overshoot OS = {OS_g:.2f} %")
print(f"Settling time Ts(5%) = {Ts_g:.2f} s")
print(f"Rise time Tr = {Tr_g:.2f} s")

# %% 8) Z feedback loop

# 6x6 Model [V, gamma, alpha, q, theta, z]
C_q6      = np.array([[0., 0., 0., 1., 0., 0.]])
C_gamma6  = np.array([[0., 1., 0., 0., 0., 0.]])
C_z6      = np.array([[0., 0., 0., 0., 0., 1.]])

# 1) Plant seen from gamma_c
A_q6 = A + B @ (Kr * C_q6)
B_q6 = B * Kr

A_gamma6 = A_q6 + B_q6 @ (Kgamma * C_gamma6)
B_gamma6 = B_q6 * Kgamma

# %% Z feedback loop — CORRECTED VERSION

# 1) Plant: gamma_c -> z (aircraft + q loop + gamma loop)
Gz = control.ss(A_gamma6, B_gamma6, C_z6, 0)

# 2) Altimeter model (1st order sensor)
tau_alt = 0.65
H_alt = control.tf([1.0], [tau_alt, 1.0])

# 3) Plant including sensor (z -> z_meas)
Gz_measured = control.series(Gz, H_alt)

# 4) Altitude controller gain
Kz = 0.001  # Subject 82

# 5) CLOSED-LOOP on altitude ERROR: z_c - z_meas
#    gamma_c = Kz * (z_c - z_meas)
sys_z_cl = control.feedback(-Kz * Gz_measured, 1)

# 6) State-space form (for analysis & simulation)
sys_z_cl_ss = control.ss(sys_z_cl)

A_z = np.array(sys_z_cl_ss.A)
B_z = np.array(sys_z_cl_ss.B)
C_z = np.array(sys_z_cl_ss.C)
D_z = np.array(sys_z_cl_ss.D)

# 7) Poles and damping
wn_z, zeta_z, poles_z = control.damp(sys_z_cl_ss)

# --- 8) Step response z_c -> z -----------------------
t = np.linspace(0, 80, 2000)
t_step, z_resp = control.step_response(sys_z_cl_ss, t)

plt.figure(figsize=(8,4))
plt.plot(t_step, z_resp, lw=2)
plt.grid(True)
plt.xlabel("Time (s)")
plt.ylabel("z (m)")
plt.title("Closed-loop altitude response (z_c → z)")
plt.tight_layout()
plt.show()

# Performance indicators
OS_z, Tr_z, Ts_z = step_info(t_step, z_resp)
print(f"Kz = {Kz}")
print(f"Overshoot OS = {OS_z:.2f} %")
print(f"Settling time Ts(5%) = {Ts_z:.2f} s")
print(f"Rise time Tr = {Tr_z:.2f} s")
print("Poles:")
print(poles_z)
print("Damping ratios:")
print(zeta_z)


# %% Saturation

# State Space: gamma_c_sat -> alpha

# Output Matrix for alpha (3rd state)
C_alpha_cl = np.array([[0., 0., 1., 0.]])
D_alpha_cl = np.array([[0.]])

# System construction
sys_gamma_alpha = control.ss(A_gamma, B_gamma, C_alpha_cl, D_alpha_cl)
print("Closed Loop State Space: gamma_c_sat -> alpha")
print(sys_gamma_alpha)

t_alpha = np.linspace(0, 10, 1000)
T_alpha, y_alpha = control.step_response(sys_gamma_alpha, t_alpha)

plt.figure()
plt.plot(T_alpha, y_alpha, 'b', linewidth=2)
plt.grid(True)
plt.title(r'Response of $\alpha$ to a step in $\gamma_{c_{sat}}$')
plt.xlabel('Time (s)')
plt.ylabel(r'$\alpha$ (rad)')
plt.show()

# Alpha max
alphamax= 3*g*(alpha_eq-alpha0)+alpha_eq
print (alphamax)

#%%


def f(gamma_c):
    sys_gamma_alpha = control.ss(A_gamma, B_gamma*gamma_c, C_alpha_cl, D_alpha_cl*gamma_c)
    t_alpha = np.linspace(0, 10, 1000)
    T_alpha, y_alpha = control.step_response(sys_gamma_alpha, t_alpha)
    
    return max(y_alpha) - alphamax


diff = f(0.28)
print (diff) #



# %% Dichotomy algorithm 

def dicho(a,b, eps=1e-6 ):
    gamma_c=0
    it=0
    while (b-a)>eps:
        gamma_c=(a+b)/2
        # print() # verbose
        if f(a)*f(gamma_c)<=0:
            b=gamma_c
        else :
            a=gamma_c
        it+=1
    return gamma_c ,it

print (dicho(0,np.pi/2))
gammacmax,_=dicho(0,np.pi/2)


#%%




































#%%
# %% 9) Flight Management Simulation
print("\nFLIGHT MANAGEMENT SIMULATION")

# 1. System Preparation
sys_phase_gamma = control.ss(A_gamma6, B_gamma6, np.eye(6), np.zeros((6,1)))
sys_phase_z = sys_z_cl_ss 

# 2. Parameters
gamma_climb = gammacmax * (np.pi/180)

# Phase 1: Ascent
gamma_climb = gammacmax

# Phase 1: Ascent
duration_ascent = 40.0
# Phase 2: Cruise
duration_cruise = 100.0
# Phase 3: Descent
duration_descent = 40.0
# Phase 4: Flare
duration_flare = 30.0
z_ground = 0

# 3. Sequential Simulation

# Phase 1: Ascent (Gamma)
t1 = np.linspace(0, duration_ascent, 400)
u1 = np.full_like(t1, gamma_climb)
x0_1 = np.zeros(6) 
resp1 = control.forced_response(sys_phase_gamma, T=t1, U=u1, X0=x0_1)
x_end_1 = resp1.states[:, -1]

# Phase 2: Cruise (Z)
t2 = np.linspace(0, duration_cruise, 1000)
z_target_cruise = x_end_1[5]
u2 = np.full_like(t2, z_target_cruise)
x0_2 = np.append(x_end_1, x_end_1[5]) # 6 states + sensor
resp2 = control.forced_response(sys_phase_z, T=t2, U=u2, X0=x0_2)
x_end_2 = resp2.states[:, -1]

# Phase 3: Descent (Gamma)
t3 = np.linspace(0, duration_descent, 400)
u3 = np.full_like(t3, -gamma_climb) # Symetric descent
x0_3 = x_end_2[0:6]
resp3 = control.forced_response(sys_phase_gamma, T=t3, U=u3, X0=x0_3)
x_end_3 = resp3.states[:, -1]

# Phase 4: Flare (Z)
t4 = np.linspace(0, duration_flare, 300)
u4 = np.full_like(t4, x_end_3[5] - 20)
x0_4 = np.append(x_end_3, x_end_3[5])
resp4 = control.forced_response(sys_phase_z, T=t4, U=u4, X0=x0_4)

# 4. Plotting

# Global vector creation
# Time shift needed for continuity
time_all = np.concatenate([
    resp1.time, 
    resp2.time + resp1.time[-1], 
    resp3.time + resp2.time[-1] + resp1.time[-1],
    resp4.time + resp3.time[-1] + resp2.time[-1] + resp1.time[-1]
])

# Extract variables of interest
# Note: resp2/resp4 have 7 states, resp1/resp3 have 6.
# Concatenating only common physical states (indices 0 to 5)
states_1 = resp1.states
states_2 = resp2.states[0:6, :] # Ignore sensor state
states_3 = resp3.states
states_4 = resp4.states[0:6, :]

states_all = np.concatenate([states_1, states_2, states_3, states_4], axis=1)

# Get altitude (index 5) and gamma (index 1)
z_all = states_all[5, :]
gamma_all = states_all[1, :]
theta_all = states_all[4, :]

# PLOT GLOBAL
plt.figure(figsize=(10, 8))

# Altitude curve
plt.subplot(3, 1, 1)
plt.plot(time_all, z_all, 'b', lw=2)
plt.ylabel('Altitude z (m)')
plt.title('Flight Management Simulation')
plt.grid(True)
# Add colored zones for phases
plt.axvspan(0, t1[-1], color='green', alpha=0.1, label='Ascent (Gamma)')
plt.axvspan(t1[-1], t1[-1]+t2[-1], color='blue', alpha=0.1, label='Cruise (Z)')
plt.axvspan(t1[-1]+t2[-1], t1[-1]+t2[-1]+t3[-1], color='orange', alpha=0.1, label='Descent (Gamma)')
plt.legend()

# Gamma curve
plt.subplot(3, 1, 2)
plt.plot(time_all, gamma_all * 180/np.pi, 'r', lw=2)
plt.ylabel('Gamma (deg)')
plt.grid(True)

# Theta curve
plt.subplot(3, 1, 3)
plt.plot(time_all, theta_all * 180/np.pi, 'k', lw=2)
plt.ylabel('Theta (deg)')
plt.xlabel('Time (s)')
plt.grid(True)

plt.tight_layout()
plt.show()
# %% 10) Alternative Final Flare
# Alternative 1: Ramp
print("\n--- Alternative 1: Ramp Command ---")

# Parameters
z_start_flare = x_end_3[5]
z_target_flare = z_start_flare - 50
duration_alt_flare = 40.0

t_ramp = np.linspace(0, duration_alt_flare, 400)
z_ramp = z_start_flare + (z_target_flare - z_start_flare) * (t_ramp / duration_alt_flare)
x0_ramp = np.append(x_end_3, x_end_3[5])
resp_ramp = control.forced_response(sys_phase_z, T=t_ramp, U=z_ramp, X0=x0_ramp)

z_resp_ramp = resp_ramp.states[5, :]
gamma_resp_ramp = resp_ramp.states[1, :]

print(f"  Initial: {z_start_flare:.1f} m")
print(f"  Target: {z_target_flare:.1f} m")
print(f"  Response: {z_resp_ramp[-1]:.1f} m")


# --- Alternative 2: Exponential Command ---
print("\n--- Alternative 2: Exponential Command ---")

t_exp = np.linspace(0, duration_alt_flare, 400)
# Exp decay
tau_exp = 10.0
z_exp = z_target_flare + (z_start_flare - z_target_flare) * np.exp(-t_exp / tau_exp)

x0_exp = np.append(x_end_3, x_end_3[5])
resp_exp = control.forced_response(sys_phase_z, T=t_exp, U=z_exp, X0=x0_exp)

z_resp_exp = resp_exp.states[5, :]
gamma_resp_exp = resp_exp.states[1, :]

print(f"  tau: {tau_exp} s")
print(f"  Response: {z_resp_exp[-1]:.1f} m")


# Comparison Plot
fig, axes = plt.subplots(3, 2, figsize=(14, 10))

# Ramp Plot
axes[0, 0].plot(t_ramp, z_ramp, 'g--', lw=2, label='Command')
axes[0, 0].plot(t_ramp, z_resp_ramp, 'b-', lw=2, label='Response')
axes[0, 0].set_ylabel('Altitude z (m)')
axes[0, 0].set_title('Ramp Command')
axes[0, 0].legend()
axes[0, 0].grid(True)

axes[1, 0].plot(t_ramp, gamma_resp_ramp * 180/np.pi, 'r-', lw=2)
axes[1, 0].set_ylabel('Gamma (deg)')
axes[1, 0].grid(True)

axes[2, 0].plot(t_ramp, z_ramp - z_resp_ramp, 'm-', lw=2)
axes[2, 0].set_ylabel('Tracking Error (m)')
axes[2, 0].set_xlabel('Time (s)')
axes[2, 0].grid(True)

# Exp Plot
axes[0, 1].plot(t_exp, z_exp, 'g--', lw=2, label='Command')
axes[0, 1].plot(t_exp, z_resp_exp, 'b-', lw=2, label='Response')
axes[0, 1].set_ylabel('Altitude z (m)')
axes[0, 1].set_title(f'Exponential Command (tau={tau_exp}s)')
axes[0, 1].legend()
axes[0, 1].grid(True)

axes[1, 1].plot(t_exp, gamma_resp_exp * 180/np.pi, 'r-', lw=2)
axes[1, 1].set_ylabel('Gamma (deg)')
axes[1, 1].grid(True)

axes[2, 1].plot(t_exp, z_exp - z_resp_exp, 'm-', lw=2)
axes[2, 1].set_ylabel('Tracking Error (m)')
axes[2, 1].set_xlabel('Time (s)')
axes[2, 1].grid(True)

plt.suptitle('Alternative Final Flare Comparison: Ramp vs Exponential', fontsize=14)
plt.tight_layout()
plt.show()

print("\nTP COMPLETED - Subject 82")