import sympy as sp
import numpy as np
from DenavHart import T0_6, q
import matplotlib.pyplot as plt
import Functions as fn
pii = sp.pi
#initial and final transformation matrix
P0 = (T0_6.subs({q[0]:-pii/3, q[1]:pii/4, q[2]:pii/4, q[3]:pii/3, q[4]:pii/3, q[5]:pii/3})).evalf()
Pf = (T0_6.subs({q[0]: pii/3, q[1]:pii/3, q[2]:pii/6, q[3]:pii/3, q[4]:pii/3, q[5]:pii/3})).evalf()
t = np.linspace(0,10,100)
x, velx, accx, y, vely, accy, z, velz, accz = fn.traj_quintic(P0, Pf)
xx, yy, zz = fn.traj_line(P0, Pf)
plt.figure(1)
fn.plot_line(xx,'x')
plt.figure(2)
fn.plot_line(yy,'y')
plt.figure(3)
fn.plot_line(zz,'z')
plt.figure(4)
fn.plot_quintic(x, velx, accx, 'x')
plt.figure(5)
fn.plot_quintic(y, vely, accy, 'y')
plt.figure(6)
fn.plot_quintic(z, velz, accz, 'z')



