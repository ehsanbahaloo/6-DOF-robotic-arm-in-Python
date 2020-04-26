import sympy as sp
import Functions as fn
import matplotlib.pyplot as plt
pii = sp.pi
D_H = fn.D_H
q = fn.q
T0_6 = fn.H_joints_to_base(D_H)[5]
#initial, via and final points
P0 = (T0_6.subs({q[0]:   pii/3, q[1]:pii/3, q[2]:pii/4, q[3]:pii/3, q[4]:pii/3, q[5]:pii/3})).evalf()
Pv = (T0_6.subs({q[0]:   pii/2, q[1]:pii/4, q[2]:pii/3, q[3]:pii/4, q[4]:pii/4, q[5]:pii/3})).evalf()
Pf = (T0_6.subs({q[0]:2**pii/3, q[1]:pii/6, q[2]:pii/4, q[3]:pii/4, q[4]:pii/4, q[5]:pii/3})).evalf()

#extracting position information
x0, y0, z0,xv, yv, zv, xf, yf, zf = P0[0,3], P0[1,3], P0[2,3],\
 Pv[0,3], Pv[1,3], Pv[2,3], Pf[0,3], Pf[1,3], Pf[2,3]
 
#calculating trajectories
x, v_x, a_x = fn.quintic_trajectory_with_via_point(x0, xv, xf)
y, v_y, a_y = fn.quintic_trajectory_with_via_point(y0, yv, yf)
z, v_z, a_z = fn.quintic_trajectory_with_via_point(z0, zv, zf)

#plottng the results
plt.figure(1)
fn.plot_quintic_traj_with_via_point(x, v_x, a_x, 'x')
plt.figure(2)
fn.plot_quintic_traj_with_via_point(y, v_y, a_y, 'y')
plt.figure(3)
fn.plot_quintic_traj_with_via_point(z, v_z, a_z, 'z')
