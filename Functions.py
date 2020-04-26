import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

alpha = sp.Matrix(sp.symbols('alpha0:6'))
a = sp.Matrix(sp.symbols('a0:6'))
d = sp.Matrix(sp.symbols('d1:7'))
q = sp.Matrix(sp.symbols('q1:7'))
qDot = sp.Matrix(sp.symbols('qDot1:7'))
t = sp.Symbol('t')


D_H = { alpha[0]:        0, a[0]:      0,   d[0]:  0.75,   q[0]:          q[0],
        alpha[1]: -sp.pi/2, a[1]:   0.35,   d[1]:     0,   q[1]: -sp.pi/2+q[1],
        alpha[2]:        0, a[2]:   1.25,   d[2]:     0,   q[2]:          q[2],
        alpha[3]: -sp.pi/2, a[3]: -0.055,   d[3]:   1.1,   q[3]:          q[3],
        alpha[4]:  sp.pi/2, a[4]:      0,   d[4]:     0,   q[4]:          q[4],
        alpha[5]: -sp.pi/2, a[5]:      0,   d[5]:     0,   q[5]:          q[5]}

Oc = [[0, 0, 0.4, 1], [0, 0, 0.44751, 1], [0.18842, 0.18344, -0.042799, 1],\
      [0.27146, -0.007326, 0, 1], [0, 0, 0, 1], [0, 0 ,0, 1]]#mass centroids of each links

mass = [150, 138, 95, 71, 17, 7, 0.5]

Inertia = {'ixx':[ 60, 30, 40,  1, 0.18,     0,     0],
           'iyy':[ 70, 50, 40, 10, 0.55, 0.068, 0.003],
           'izz':[100, 50, 10, 10, 0.64, 0.068, 0.003]} 

V = sp.Matrix(sp.symbols('V1:7'))
U = sp.Matrix(sp.symbols('U1:7'))

Jm = [ 0.014,  0.014,  0.014, 0.010, 0.004,  0.002]
r =  [   150,    150,    100,    70,    70,     70]
km = [     1,      1,      1,     1,     1,      1]
R =  [     1,      1,      1,     5,    14,     18]



def homogeneous_transform(q, d, a, alpha):
    """ Creates Homogeneous Transform Matrix from DH parameters (based on Craige's book)"""
    T = sp.Matrix([[                sp.cos(q),              -sp.sin(q),              0,                a],
                  [ sp.sin(q) * sp.cos(alpha), sp.cos(q)*sp.cos(alpha), -sp.sin(alpha), -sp.sin(alpha)*d],
                  [   sp.sin(q)*sp.sin(alpha), sp.cos(q)*sp.sin(alpha),  sp.cos(alpha),  sp.cos(alpha)*d],
                  [                         0,                       0,              0,                1]])
    return T

def H_joints_to_base(D_H):
    """This function gets Denavit_Hartenberg parameters and returns T0_6, T3_6 and T0_3"""
    T = []
    alpha = sp.Matrix(sp.symbols('alpha0:6'))
    a = sp.Matrix(sp.symbols('a0:6'))
    d = sp.Matrix(sp.symbols('d1:7'))
    q = sp.Matrix(sp.symbols('q1:7'))
    for i in range(6):
        T.append(homogeneous_transform(q[i], d[i], a[i], alpha[i]).subs(D_H))
    T0_1 = T[0]
    T0_2 = T0_1*T[1]
    T0_3 = T0_2*T[2]
    T0_4 = T0_3*T[3]
    T0_5 = T0_4*T[4]
    T0_6 = T0_5*T[5]
    T3_6 = T[3]*T[4]*T[5]
    return [T0_1, T0_2, T0_3, T0_4, T0_5, T0_6, T3_6]

def reach_wspace(Q_range):
    """This function gets range of motion of each joint and returns workspace of the robot."""
    n = 6
    l = []
    u = []
    for i in range(n):
        l.append(Q_range[0,i])
        u.append(Q_range[1,i])
    WS = []
    for k in range(1):
        for i in range(l[1],u[1],5):
            for j in range(l[2], u[2],5):
                O6 = (H_joints_to_base(D_H)[5])[0:3,3] 
                WS.append(O6.subs({q[0]:np.deg2rad(0),\
                  q[1]:np.deg2rad(i), q[2]:np.deg2rad(j),\
                  q[3]:np.deg2rad(0), q[4]:np.deg2rad(0),\
                  q[5]:np.deg2rad(0)}))
    return(WS)

def geo_jacobian(k):
    """This function caculates the geometric jacobian of joint k_th"""
    z, o = [], []
    for i in range(6):
        z.append(sp.trigsimp(H_joints_to_base(D_H)[i][0:3,2]))
        o.append(H_joints_to_base(D_H)[i][0:3,3])
    J = sp.zeros(6,6)
    if k==6:
        o.append(o[5]) 
    for j in range(k):
        J[j] = sp.trigsimp(sp.Matrix.vstack(z[j].cross(o[k]-o[j]), z[j]))
    return J

def anal_jacobian(DH):
    """based on roll-pitch-yaw angles"""
    zero = sp.zeros(3)
    I = sp.eye(3)
    r, p = sp.symbols('r, p') #r--> roll & p-->pitch 
    B = sp.Matrix([[1, 0, sp.sin(p)],
                   [0, sp.cos(r), -sp.cos(p)*sp.sin(r)],
                   [0, sp.sin(r), sp.cos(p)*sp.cos(r)]])
    J_a = (sp.Matrix.vstack(sp.Matrix.hstack(I, zero), sp.Matrix.hstack(zero, sp.trigsimp(B.inv()))) * geo_jacobian(6))
    return(J_a)
    
def Matrix_Equality(M1, M2):
    """This Function checks whwther two Matrix are equal"""
    for i in range(3):
        for j in range(3):
            if not (M1[i,j]-M2[i,j] < 10**(-6)):
                return False
    return True

def first_three_angles(Xc, Yc, Zc):  
    """This function gets position information of end effector an return the first three anglse"""
    theta_1 = sp.atan2(Yc, Xc)
    a = ((1.1)**2 + (0.055)**2)**0.5 
    b = ((((Xc**2 + Yc**2)**0.5 - 0.350)**2) + ((Zc - 0.750)**2))**0.5 
    c = 1.25
#    angle_c = acos((a**2 + b**2 - c**2)/(2*a*b))
    angle_a = sp.acos((c**2 + b**2 - a**2)/(2*c*b))
    angle_b = sp.acos((a**2 + c**2 - b**2)/(2*a*c))
    theta_2 = (sp.pi/2-(angle_a + sp.atan2(Zc - 0.750, (Xc*2 + Yc**2)**0.5 - 0.350))).evalf()
    theta_3 = (sp.pi/2-(angle_b + sp.atan2(0.055,1.1))).evalf()
    return [theta_1, theta_2, theta_3]

def last_three_angles(R_3_6):
    """This function gets orientation information of end effector an return the first three anglse"""
    R3_6 = H_joints_to_base(D_H)[6]
    if (R_3_6[2,2])!=0 or (R_3_6[0,2])!=0:
        theta_5 = (sp.atan2(((1-R_3_6[1,2]**2)**0.5 ),R_3_6[1,2])).evalf()
        theta_4 = (sp.atan2(R_3_6[2,2], -R_3_6[0,2])).evalf()
        theta_6 = (sp.atan2(-R_3_6[1,1], R_3_6[1,0])).evalf()
        if Matrix_Equality(R3_6.subs({q[3]: theta_4, q[4]: theta_5, q[5]: theta_6}), R_3_6.evalf()):
            theta_4 = theta_4
            theta_5 = theta_5
            theta_6 = theta_6
        else:
            theta_5 = (sp.atan2((-(1-R_3_6[1,2]**2)**0.5 ) ,R_3_6[1,2])).evalf()
            theta_4 = (sp.atan2(-R_3_6[2,2], R_3_6[0,2])).evalf()
            theta_6 = (sp.atan2(R_3_6[1,1], -R_3_6[1,0])).evalf()
            if Matrix_Equality(R3_6.subs({q[3]: theta_4, q[4]: theta_5, q[5]: theta_6}), R_3_6.evalf()):
                theta_4 = theta_4
                theta_5 = theta_5
                theta_6 = theta_6        
    if (R_3_6[2,0])==0 and (R_3_6[2,1])==0 and R_3_6[2,2]==1:
        theta_5 = 0
        theta_4 = (sp.atan2(R_3_6[0,0], -R_3_6[2,0])).evalf()
        theta_6 = 0
    if R_3_6[2,0]==0 and R_3_6[2,1]==0 and R_3_6[2,2]==-1:
        theta_5 = (sp.pi).evalf()
        theta_4 = (sp.atan2(-R_3_6[0,0], -R_3_6[2,0])).evalf()
        theta_6 = 0
    return [theta_4, theta_5, theta_6]

def inertia_of_links():
    """This function returnss inertial matrix of each joint."""
    Inertia = {'ixx':[ 60, 30, 40,  1, 0.18,     0,     0],
               'iyy':[ 70, 50, 40, 10, 0.55, 0.068, 0.003],
               'izz':[100, 50, 10, 10, 0.64, 0.068, 0.003]}#information extracted from kuka_kr210 urdf   
    I = []
    for i in range(6):
        I.append(sp.diag(Inertia['ixx'][i], Inertia['iyy'][i], Inertia['izz'][i]))
    return I
def mass_matrix(): 
    I = inertia_of_links()
    D = sp.zeros(6,6)
    for i in range(6):
        J = geo_jacobian(i+1)
        Jv = J[0:3,:]
        Jw = J[3:6,:]
        D = (mass[i]*sp.Transpose(Jv)*Jv + sp.Transpose(Jw)*I[i]*Jw)
    return D
#    q = sp.symbols('q:'+str(7))

def coriolis_and_centrifugal_matrix():
    C = sp.zeros(6,6)  
    D = mass_matrix()
    for k in range(6):
        for j in range(6):
            c_jk = 0
            for i in range(6):
                c_jk += 1/2*sp.trigsimp((sp.diff(D[k,j],q[i])+\
                sp.diff(D[k,i],q[j])-sp.diff(D[i,j],q[k])))*qDot[i]
            C[k,j] = c_jk
    return C
def gravity_matrix(): 
    H = H_joints_to_base(D_H)
    P = sp.Matrix([0])
    g = 9.81*sp.Matrix([0, 0, -1])
    for i in range(6):
        oc = sp.Matrix(Oc[i])
        rr = H[i]*oc
        r = sp.Matrix(rr[0:3])
        P += mass[i]*sp.Transpose(g)*r
    G = sp.zeros(6,1)
    for i in range(6):
        G[i,0] = sp.diff(P,q[i])
    return G



def points_line(x0,xf, t0, tf):
        """This function gets both the initial and final x(or y or z) and returns points"""
        """via this points based on the unimportance of velocity and acceleration on initial and final points"""   
        X = sp.Matrix([x0, xf])#vector of initial and final point
        Coef = sp.Matrix([[1, t0],
                          [1, tf]])#matrix of coefficients
        a = Coef.inv()*X #calculation of unknown coefficients of polynominal
        t = np.linspace(t0,tf,100)#initial and final time
        x = []#defining a list to append Xs via initial and final points
        for i in range(100):#for loop for calculation of aforementioned points
            x_t = a[0] + a[1]*t[i]
            x.append(x_t)
        return x

def traj_line(P0, Pf):
    """This function gets two transformation matrices and returns a trajectory without any"""
    """importance on velocity and acceleration on initial and final points"""
    x0, y0, z0, xf, yf, zf = P0[0,3], P0[1,3], P0[2,3], Pf[0,3], Pf[1,3], Pf[2,3]#getting initial and final x, y, z
    t0, tf = 0, 10 #initial and final time
    x = points_line(x0,xf,t0,tf)#calculation of Xs
    y = points_line(y0,yf,t0,tf)#calculation of Ys
    z = points_line(z0,zf,t0,tf)#calculation of Zs
    return x, y, z

def points_quintic(x0, xf, t0, tf):
        """This function gets both the initial and final x(or y or z) and returns points"""
        """via this points based on the zero velocity and zero acceleration on initial and """
        """final points (quintic polynominal)""" 
        x0dot, x0ddot, xfdot, xfddot = 0, 0, 0, 0 #zero velocity and acceleration of initial and final points
        X = sp.Matrix([x0, x0dot, x0ddot, xf, xfdot, xfddot])#vector of initial and final point
        t0, tf = 0, 10 #initial and final time 
        Coef = sp.Matrix([[1, t0, t0**2,   t0**3,    t0**4,    t0**5],
                          [0,  1,  2*t0, 3*t0**2,  4*t0**3,  5*t0**4],
                          [0,  0,     2,    6*t0, 12*t0**2, 20*t0**3],
                          [1, tf, tf**2,   tf**3,    tf**4,    tf**5],
                          [0,  1,  2*tf, 3*tf**2,  4*tf**3,  5*tf**4],
                          [0,  0,     2,    6*tf, 12*tf**2, 20*tf**3]])#matrix of coefficients
        a = Coef.inv()*X #calculation of unknown coefficients of polynominal
        t = np.linspace(t0,tf,100)
        x, vel, acc = [], [], []#defining blank lists to append X's, Velocity and acceleration of x's via initial and final points
        for i in range(100):#caculation of via X's due to a cubic polynominal
            x_t = a[0] + a[1]*t[i] + a[2]*t[i]**2 + a[3]*t[i]**3 + a[4]*t[i]**4 + a[5]*t[i]**5 
            v_t = a[1] + 2*a[2]*t[i] + 3*a[3]*t[i]**2 + 4*a[4]*t[i]**3 + 5*a[5]*t[i]**4
            a_t = 2*a[2] + 6*a[3]*t[i] + 12*a[4]*t[i]**2 + 20*a[5]*t[i]**3
            x.append(x_t)#appending via X's to the aforementioned blank list
            vel.append(v_t)
            acc.append(a_t)
        return x, vel, acc
def traj_quintic(P0, Pf):
    """This function gets two transformation matrices and returns a trajectory with"""
    """zero velocity and zero acceleration on initial and final points"""
    x0, y0, z0, xf, yf, zf = P0[0,3], P0[1,3], P0[2,3], Pf[0,3], Pf[1,3], Pf[2,3]#getting initial and final x, y, z
    t0, tf = 0, 10
    x, vel_x, acc_x = points_quintic(x0, xf, t0, tf)#calculation of X's(X, velocity and accelerationof X)
    y, vel_y, acc_y = points_quintic(y0, yf, t0, tf)
    z, vel_z, acc_z = points_quintic(z0, zf, t0, tf)
    return x, vel_x, acc_x, y, vel_y, acc_y, z, vel_z, acc_z

def plot_line(x, name):
    """plotting due to traj_wva function"""
    time = np.linspace(0,10,100)
    plt.plot(time, x, label=name)
    plt.xlabel('time')
    plt.ylabel(name)
    plt.legend()
    plt.show()
    
def plot_quintic(x, v_x, a_x, name):
    """plotting due to traj_va function"""
    time = np.linspace(0,10,100)
    plt.plot(time, x, label=name)
    plt.plot(time,v_x, label=name+'_velocity')
    plt.plot(time,a_x, label=name+'_acceleration')
    plt.xlabel('time')
    plt.ylabel(name)
    plt.legend()
    plt.show()         

def Joint_Space_Traj(Xc, Yc, Zc):
    """This function gets points a task space trajectory and returns joint space trajectory."""
    theta = []
    for i in range(len(Xc)):
        theta.append(first_three_angles(Xc[i], Yc[i], Zc[i]))
    q = [[],[],[]]
    for i in range(100):
        t = theta[i]
        for j in range(3):
            q[j].append(t[j])
    return q

def points_quintic_with_via_point(x1_t, x2_t):
    time = np.linspace(0,5,500)
    x = []
    for i in range(len(time)):
        x.append(x1_t.subs({t:time[i]}))
    for j in range(len(time)):
        x.append(x2_t.subs({t:time[j]}))
    return x
def quintic_trajectory_with_via_point(x0, xv, xf, t0=0, tf=5):
    X = sp.Matrix([x0, xv, xv, xf, 0, 0, 0, 0, 0, 0, 0, 0])
    Coef = sp.Matrix([[1,  0,     0,       0,        0,        0, 0,  0,     0,       0,        0,        0],
                      [1, tf, tf**2,   tf**3,    tf**4,    tf**5, 0,  0,     0,       0,        0,        0],
                      [0,  0,     0,       0,        0,        0, 1,  0,     0,       0,        0,        0],
                      [0,  0,     0,       0,        0,        0, 1, tf, tf**2,   tf**3,    tf**4,    tf**5],
                      [0,  1,     0,       0,        0,        0, 0,  0,     0,       0,        0,        0],
                      [0,  0,     0,       0,        0,        0, 0,  1,  2*tf, 3*tf**2,  4*tf**3,  5*tf**4],
                      [0,  0,     1,       0,        0,        0, 0,  0,     0,       0,        0,        0],
                      [0,  0,     0,       0,        0,        0, 0,  0,     2,    6*tf, 12*tf**2, 20*tf**3],
                      [0,  1,  2*tf, 3*tf**2,  4*tf**3,  5*tf**4, 0, -1,     0,       0,        0,        0],
                      [0,  0,     2,    6*tf, 12*tf**2, 20*tf**3, 0,  0,    -2,       0,        0,        0],
                      [0,  0,     0,       6,    24*tf, 60*tf**2, 0,  0,     0,      -6,        0,        0],
                      [0,  0,     0,       0,       24,   120*tf, 0,  0,     0,       0,      -24,        0]])
    a = Coef.inv()*X
    x1_t = sp.Matrix(a[0:6]).T*sp.Matrix([1, t, t**2,   t**3,    t**4,    t**5])
    x2_t = sp.Matrix(a[6:12]).T*sp.Matrix([1, t, t**2,   t**3,    t**4,    t**5])
    v1_t = sp.diff(x1_t, t)
    v2_t = sp.diff(x2_t, t)
    a1_t = sp.diff(v1_t, t)
    a2_t = sp.diff(v2_t, t)
    x = points_quintic_with_via_point(x1_t, x2_t)
    v = points_quintic_with_via_point(v1_t, v2_t)
    a = points_quintic_with_via_point(a1_t, a2_t)
    return x, v, a


def plot_quintic_traj_with_via_point(x, v_x, a_x, name):
    time = np.linspace(0,10,1000)
    plt.plot(time, x, label=name)
    plt.plot(time,v_x, label=name+'_velocity')
    plt.plot(time,a_x, label=name+'_acceleration')
    plt.xlabel('time')
    plt.ylabel(name)
    plt.legend()
    plt.show()