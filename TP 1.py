
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
# -----------------------------
# CONSTANTS
# -----------------------------

g = 9.81                     # gravitational acceleration (m/s^2)
m = 8400                     # aircraft mass (kg)
S = 34                       # wing reference area (m^2)
alt_ft = 21230               # altitude (ft) - subject number 82
alt_m = alt_ft * 0.3048      # altitude converted to meters (m)
M = 0.95                     # Mach number (dimensionless) - subject number 82
R = 287.05                   # specific gas constant for air (J/(kg·K))
l_ref = 5.24
rg = 2.65
c = 0.52 

# -----------------------------
# GRAPHICAL DATA (from aerodynamic plots)
# -----------------------------

Cx0 = 0.018                  # zero-lift drag coefficient Cx0 (-) [M=0.95]
k = 0.28                     # polar coefficient k (-) [M=0.95]
Cz_alpha = 2.85              # lift curve slope w.r.t incidence Cz_alpha (1/rad) [M=0.95]
Cz_delta_m = 0.75            # lift curve slope w.r.t elevator Cz_delta_m (1/rad) [M=0.95]
delta_m0 = 0.002             # equilibrium elevator deflection for null lift δm0 (rad) [M=0.95]
alpha0 = 0.013               # incidence for zero lift at δm = 0 α0 (rad) [M=0.95]
f = 0.545                    # aerodynamic center coefficient of body and wings (−) [M=0.95]
f_delta = 0.82               # aerodynamic center coefficient of fins (−) [M=0.95]
Cmq = -0.55                  # damping coefficient (s/rad) [M=0.95]


# -----------------------------
# COMPUTED DATA 
# -----------------------------

# Values will be computed from US Standard Atmosphere 1976 model
# For alt = 21230 ft = 6471 m, using interpolation from atm_std:
hgeo, rho, a = get_cte_atm(alt_m)
T = 249.2 - (alt_m - 6000) * (249.2 - 236.2) / 2000  # linear interpolation ~245 K
P = rho * R * T               # ideal gas law
V_eq = M * a                 # true airspeed (m/s)
Q = 0.5 * rho * V_eq**2      # dynamic pressure (Pa = N/m^2)
l_t = (3/2) * l_ref          # total aircraft length lt = 3/2 * l_ref (m)
XF = -f * l_t
XG = -c * l_t
X = XF - XG                  # aerodynamic center position of wings/body (m)
XF_delta = -f_delta * l_t
Y = XF_delta - XG           # aerodynamic center position of fins (m)


# -----------------------------
# ALGORITHM FOR COMPUTING THE EQUILIBRIUM POINT
# -----------------------------
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
# -----------------------------
# Additional needed definitions
# -----------------------------

gamma_eq = 0.0                     # level flight (rad)
F_tau = 0.0                        # no thrust variation in simplified model

F_eq = Q * S * Cxeq / np.cos(alpha_eq)


# Aerodynamic derivatives
Cx_alpha = 2 * k * Czeq * Cz_alpha
CN_alpha = Cx_alpha * np.sin(alpha_eq) + Cz_alpha * np.cos(alpha_eq)
CN_delta_m = Cx_delta_m * np.sin(alpha_eq) + Cz_delta_m * np.cos(alpha_eq)

Cm_alpha = (X / l_ref) * CN_alpha
Cm_delta_m = (Y / l_ref) * CN_delta_m

Iyy = m * rg**2 


# -----------------------------
# X-coefficients of simplified longitudinal model
# -----------------------------
X_V = (2 * Q * S * Cxeq) / (m * V_eq)
X_alpha = (F_eq * np.sin(alpha_eq)) / (m * V_eq) + (Q * S * Cx_alpha) / (m * V_eq)
X_gamma = (g * np.cos(gamma_eq)) / V_eq
X_delta_m = (Q * S * Cx_delta_m) / (m * V_eq)
X_tau = -(F_tau * np.cos(alpha_eq)) / (m * V_eq)

# -----------------------------
# m-coefficients of simplified longitudinal model
# -----------------------------
m_V = 0
m_alpha = (Q * S * l_ref * Cm_alpha) / Iyy
m_q = (Q * S * l_ref**2 * Cmq) / (V_eq * Iyy)
m_delta_m = (Q * S * l_ref * Cm_delta_m) / Iyy


# -----------------------------
# Z-coefficients of simplified longitudinal model
# -----------------------------
Z_V = (2 * Q * S * Czeq) / (m * V_eq)
Z_alpha = (F_eq * np.cos(alpha_eq) / (m * V_eq)) + (Q * S * Cz_alpha) / (m * V_eq)
Z_gamma = (g * np.sin(gamma_eq)) / V_eq
Z_delta_m = (Q * S * Cz_delta_m) / (m * V_eq)
Z_tau = (F_tau * np.sin(alpha_eq)) / (m * V_eq)

# -----------------------------------------
# MATRIX 
# -----------------------------------------
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


title(r"Step response $\alpha/\delta m$ et $q/\delta m$")
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

sisotool(TqDm_ss)

#%% 5) CLOSED LOOP WITH GAIN Kr page 26

# Extract the q output: C_q selects only the q state (2nd state)
# For the 2x2 system Ai, the output is the 2nd state (theta or q depending on formulation)
Cq = np.matrix([[0, 1]])  # Output matrix for q only (1x2 for 2x2 system)
# Choose the damping ratio: xi = 0.65
xi_desired = 0.65
TqDm_tf = control.ss2tf(control.ss(Ai, Bi, Cq, 0))
# Kr = sisotool(TqDm_tf, kmin=0.0001, kmax=40.0, kdefault=0.010, xispec=xi_desired)
Kr = -0.10  # Subject 82: gives damping ratio xi ≈ 0.65
# ==========================================
# CLOSED LOOP STATE SPACE REPRESENTATION
# ==========================================
# According to slide: A_k = A - Kr*B*C_q
#                     B_k = Kr*B
#                     C_k = C_out
#                     D_k = 0

# A_k = A - Kr * B * C_q (from feedback: X = A_k*X + B_k*Y_c)
Ak = Ai - Kr * Bi @ Cq
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
print("\n" + "="*70)
print("Testing closed-loop step response")
print("="*70)

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

tau = 1                                     # tau > 0.5 s
W = control.tf([tau, 0], [tau, 1])          # Washout filter
G_fwd = Kr * TaDm_tf                        # Forward path from qc to alpha

# Open-loop system (from qc to alpha) - no feedback
T_open = G_fwd

# Closed-loop without filter (direct q feedback)
# Feedback H = TqDm_tf (q path). Closed loop: G_fwd / (1 + G_fwd * H)
T_cl_no_filter = control.feedback(G_fwd, TqDm_tf)

# Closed-loop with washout: feedback uses W * TqDm_tf
H_wash = control.series(W, TqDm_tf)
T_cl_wash = control.feedback(G_fwd, H_wash)

# Compute DC gains to compare steady-state alpha gains
dc_open = control.dcgain(T_open)
dc_no_filter = control.dcgain(T_cl_no_filter)
dc_wash = control.dcgain(T_cl_wash)


# Note: W(0)=0, so the washout removes the DC contribution of the q feedback
# therefore T_cl_wash DC gain should be equal to T_open DC gain (within numerical precision)

# Simulate step responses (unit step in qc)
t = np.arange(0, 10, 0.01)
Y_open, T_open_t = control.matlab.step(T_open, t)
Y_no_filter, T_nf_t = control.matlab.step(T_cl_no_filter, t)
Y_wash, T_w_t = control.matlab.step(T_cl_wash, t)

# Plot comparison
plt.figure(figsize=(8,5))
plt.plot(T_open_t, Y_open, 'k--', lw=2, label='Open-loop (Kr*Ta)')
plt.plot(T_nf_t, Y_no_filter, 'r-', lw=2, label='Closed-loop (no washout)')
plt.plot(T_w_t, Y_wash, 'b-', lw=2, label=f'Closed-loop (washout tau={tau}s)')
plt.xlabel('Time (s)')
plt.ylabel(r'\alpha (rad) per unit $q_c$')
plt.title('Comparison: Open-loop / Closed-loop without filter / with washout')
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()

# %% 7) Gamma feedback loop – page 30

# ---- Modèle 4x4 : états = [V, γ, α, q]^T ----
# Ar et Br sont déjà définis plus haut dans le code.

# Matrices de sortie pour γ et q
C_gamma4 = np.array([[0., 1., 0., 0.]])  # γ
C_q4     = np.array([[0., 0., 0., 1.]])  # q

A_q = Ar - Br @ (Kr * C_q4)   # dynamique avec retour q
B_q = Br*Kr                      # entrée = commande δm_qc

# Plante vue par la boucle γ : de δm_qc à γ
sys_gamma_plant = control.ss(A_q, B_q, C_gamma4, 0)

# choix de K avec siso
#Kgamma = sisotool(sys_gamma_plant, kmin=1e-4, kmax=50, kdefault=1.0, xispec=0.7)
Kgamma = 5.0  # Subject 82: gives damping ratio xi ≈ 0.5, OS < 15%

# --- Closed-loop with gamma feedback (γc -> γ) ---

# Matrices d'état du système avec boucle q déjà fermée : A_q, B_q
# Boucle γ : δm_qc = Kgamma * (γc - γ)

A_gamma = A_q - B_q @ (Kgamma * C_gamma4)   # A_cl
B_gamma = B_q * Kgamma                      # B_cl
C_gamma = C_gamma4                          # C_cl (sortie = γ)
D_gamma = np.array([[0.]])                  # D_cl

# 1) Représentation d'état fermée
sys_gamma_cl = control.ss(A_gamma, B_gamma, C_gamma, D_gamma)

# 2) Fonction de transfert fermée (γc -> γ)
T_gamma_cl = control.ss2tf(sys_gamma_cl)

# 3) Pôles, amortissement, pulsations propres
wn_g, zeta_g, poles_g = control.damp(sys_gamma_cl)

print("Closed-loop poles (gamma loop):")
print(poles_g)
print("Damping ratios:")
print(zeta_g)
print("Natural frequencies (rad/s):")
print(wn_g)

# 4) Réponse indicielle γc → γ
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

# (optionnel mais pratique pour vérifier les specs)
OS_g, Tr_g, Ts_g = step_info(t, y_gamma)
print(f"Overshoot OS = {OS_g:.2f} %")
print(f"Settling time Ts(5%) = {Ts_g:.2f} s")
print(f"Rise time Tr = {Tr_g:.2f} s")

# %% 8) Z feedback loop – page 31

# Modèle complet 6x6 : X = [V, γ, α, q, θ, z]^T
# A, B, Kr, Kgamma déjà définis

# Sélection des sorties
C_q6      = np.array([[0., 0., 0., 1., 0., 0.]])  # q
C_gamma6  = np.array([[0., 1., 0., 0., 0., 0.]])  # γ
C_z6      = np.array([[0., 0., 0., 0., 0., 1.]])  # z

# --- 1) Avion + boucle q + boucle γ (plante vue de γc) ---

A_q6 = A - B @ (Kr * C_q6)
B_q6 = B * Kr

A_gamma6 = A_q6 - B_q6 @ (Kgamma * C_gamma6)
B_gamma6 = B_q6 * Kgamma      # entrée = γc

# Plante de γc vers z (sans capteur altitude)
Gz = control.ss(A_gamma6, B_gamma6, C_z6, 0)

# Capteur d'altitude : 1 / (1 + τ s), τ = 0.65 s
tau_alt = 0.65
H_alt = control.tf([1.0], [tau_alt, 1.0])

# Plante utilisée dans sisotool : z_c -> z_mesuré (Gz suivi du capteur)
Gz_for_siso = control.series(Gz, H_alt)

# --- 2) Choix de Kz avec sisotool -------------------------
# Kz = sisotool(Gz_for_siso, kmin=1e-4, kmax=50.0, kdefault=1.0, xispec=0.7)
Kz = 0.001  # Subject 82: Ts ≈ 5s, OS ≈ 2.5% < 15%

# --- 3) Boucle z fermée : z_c -> z ------------------------
# γc = Kz (z_c - z_mesuré), z_mesuré = H_alt(z)
# => Tz(s) = Kz Gz(s) / (1 + Kz Gz(s) H_alt(s))

G_cl = Kz * Gz
sys_z_cl = control.feedback(G_cl, H_alt)   # sortie = z

# Représentation d'état fermée (A_z, B_z, C_z, D_z)
sys_z_cl_ss = control.ss(sys_z_cl)
A_z = np.array(sys_z_cl_ss.A)
B_z = np.array(sys_z_cl_ss.B)
C_z = np.array(sys_z_cl_ss.C)
D_z = np.array(sys_z_cl_ss.D)

# --- 4) Pôles, amortissement, pulsations propres ----------
wn_z, zeta_z, poles_z = control.damp(sys_z_cl_ss)

# --- 5) Réponse indicielle z_c -> z -----------------------
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

# Indicateurs de performance
OS_z, Tr_z, Ts_z = step_info(t_step, z_resp)
print(f"Kz = {Kz}")
print(f"Overshoot OS = {OS_z:.2f} %")
print(f"Settling time Ts(5%) = {Ts_z:.2f} s")
print(f"Rise time Tr = {Tr_z:.2f} s")
print("Poles:")
print(poles_z)
print("Damping ratios:")
print(zeta_z)


#%% partie saturation 

# %% BUILD STATE SPACE: gamma_c_sat -> alpha

# 1. Define the Output Matrix for alpha
# The state vector is x = [V, gamma, alpha, q]
# alpha is the 3rd state (index 2)
C_alpha_cl = np.array([[0., 0., 1., 0.]]) # pour avoir alpha ( 3 eme position )
D_alpha_cl = np.array([[0.]]) # reste  a 0 

# 2. Construct the State Space System

sys_gamma_alpha = control.ss(A_gamma, B_gamma, C_alpha_cl, D_alpha_cl)

print("Closed Loop State Space: gamma_c_sat -> alpha")

print(sys_gamma_alpha)


t_alpha = np.linspace(0, 10, 1000)
T_alpha, y_alpha = control.step_response(sys_gamma_alpha, t_alpha)  # par défault gammac est a 1 quand on prend un step
# on peut multiplier par la valeur final voulu comme cest un state space donc on multiplie B et D 
plt.figure()
plt.plot(T_alpha, y_alpha, 'b', linewidth=2)
plt.grid(True)
plt.title(r'Response of $\alpha$ to a step in $\gamma_{c_{sat}}$')
plt.xlabel('Time (s)')
plt.ylabel(r'$\alpha$ (rad)')
plt.show()
#%%

# calcul de alpha max

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


#%% algo dichotomie 

def dicho(a,b, eps=1e-6 ):
    gamma_c=0
    it=0
    while (b-a)>eps:
        gamma_c=(a+b)/2
        print()
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
#%% 9) FLIGHT MANAGEMENT - MULTI-PHASE SIMULATION

print("\n" + "="*70)
print("FLIGHT MANAGEMENT SIMULATION")
print("="*70)

# --- 1. Préparation des systèmes ---

# Système pour le maintien de Pente (Gamma Hold)
# On utilise le modèle 6 états (V, gamma, alpha, q, theta, z) défini dans la section 8
# Entrée : gamma_c, Sortie : états complets
sys_phase_gamma = control.ss(A_gamma6, B_gamma6, np.eye(6), np.zeros((6,1)))

# Système pour le maintien d'Altitude (Z Hold)
# On utilise le système bouclé complet défini section 8
# Il a 7 états : 6 états avion + 1 état capteur (z_mesuré)
sys_phase_z = sys_z_cl_ss 


# --- 2. Définition des paramètres de vol ---

# Phase 1: Montée (Ascent)
gamma_climb_deg = gammacmax
gamma_climb = gamma_climb_deg * (np.pi/180)
duration_ascent = 40.0 # secondes

# Phase 2: Croisière (Cruise)
duration_cruise = 100.0 # secondes (demandé environ 100s)

# Phase 3: Descente (Descent)
gamma_desc_deg = -gammacmax
gamma_desc = gamma_desc_deg * (np.pi/180)
duration_descent = 40.0

# Phase 4: Flare / Palier
duration_flare = 30.0
z_ground = 0 # ou altitude finale souhaitée


# --- 3. Simulation Séquentielle ---

# --- PHASE 1 : ASCENT (Gamma Control) ---
t1 = np.linspace(0, duration_ascent, 400)
u1 = np.full_like(t1, gamma_climb) # Consigne constante

# Condition initiale : tout à 0 (point d'équilibre)
x0_1 = np.zeros(6) 

resp1 = control.forced_response(sys_phase_gamma, T=t1, U=u1, X0=x0_1)
x_end_1 = resp1.states[:, -1] # On récupère le vecteur d'état final (6 valeurs)


# --- PHASE 2 : CRUISE (Z Control) ---
t2 = np.linspace(0, duration_cruise, 1000)
z_target_cruise = x_end_1[5] # On vise l'altitude atteinte à la fin de la montée
u2 = np.full_like(t2, z_target_cruise)

# Transition X0 : Le système Z a 7 états. 
# On prend les 6 états de la fin de phase 1 + on initialise le capteur
# On suppose que le capteur lit l'altitude réelle instantanément : x_sensor = z
x0_2 = np.append(x_end_1, x_end_1[5]) 

resp2 = control.forced_response(sys_phase_z, T=t2, U=u2, X0=x0_2)
x_end_2 = resp2.states[:, -1] # On récupère le vecteur d'état final (7 valeurs)


# --- PHASE 3 : DESCENT (Gamma Control) ---
t3 = np.linspace(0, duration_descent, 400)
u3 = np.full_like(t3, gamma_desc)

# Transition X0 : Le système Gamma a 6 états.
# On prend les 6 premiers états de la phase 2 (on ignore l'état du capteur)
x0_3 = x_end_2[0:6]

resp3 = control.forced_response(sys_phase_gamma, T=t3, U=u3, X0=x0_3)
x_end_3 = resp3.states[:, -1]


# --- PHASE 4 : FLARE / LEVEL (Z Control) ---
t4 = np.linspace(0, duration_flare, 300)
u4 = np.full_like(t4, x_end_3[5] - 20) # Exemple : on demande un palier un peu plus bas ou au sol

# Transition X0 : Retour vers 7 états
x0_4 = np.append(x_end_3, x_end_3[5])

resp4 = control.forced_response(sys_phase_z, T=t4, U=u4, X0=x0_4)


# --- 4. Concaténation et Affichage ---

# Création des vecteurs globaux
# Il faut décaler le temps pour qu'il soit continu
time_all = np.concatenate([
    resp1.time, 
    resp2.time + resp1.time[-1], 
    resp3.time + resp2.time[-1] + resp1.time[-1],
    resp4.time + resp3.time[-1] + resp2.time[-1] + resp1.time[-1]
])

# Extraction des variables d'intérêt
# Attention : resp2 et resp4 ont 7 états, resp1 et resp3 ont 6 états.
# On ne concatène que les états physiques communs (indices 0 à 5)
states_1 = resp1.states
states_2 = resp2.states[0:6, :] # On ignore l'état du capteur
states_3 = resp3.states
states_4 = resp4.states[0:6, :]

states_all = np.concatenate([states_1, states_2, states_3, states_4], axis=1)

# Récupération de l'altitude (index 5) et de la pente (index 1)
z_all = states_all[5, :]
gamma_all = states_all[1, :]
theta_all = states_all[4, :]

# PLOT GLOBAL
plt.figure(figsize=(10, 8))

# Courbe Altitude
plt.subplot(3, 1, 1)
plt.plot(time_all, z_all, 'b', lw=2)
plt.ylabel('Altitude z (m)')
plt.title('Flight Management Simulation')
plt.grid(True)
# Ajout de zones colorées pour les phases
plt.axvspan(0, t1[-1], color='green', alpha=0.1, label='Ascent (Gamma)')
plt.axvspan(t1[-1], t1[-1]+t2[-1], color='blue', alpha=0.1, label='Cruise (Z)')
plt.axvspan(t1[-1]+t2[-1], t1[-1]+t2[-1]+t3[-1], color='orange', alpha=0.1, label='Descent (Gamma)')
plt.legend()

# Courbe Pente (Gamma)
plt.subplot(3, 1, 2)
plt.plot(time_all, gamma_all * 180/np.pi, 'r', lw=2)
plt.ylabel('Gamma (deg)')
plt.grid(True)

# Courbe Theta (Assiette)
plt.subplot(3, 1, 3)
plt.plot(time_all, theta_all * 180/np.pi, 'k', lw=2)
plt.ylabel('Theta (deg)')
plt.xlabel('Time (s)')
plt.grid(True)

plt.tight_layout()
plt.show()