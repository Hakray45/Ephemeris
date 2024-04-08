import numpy as np
from math import sqrt,pow

mu = 398600.4418e+9;    #- earth gravity
J = 1082625.75e-9;      #- Ephi zonal harmonic decomposition of the geopotential in a series of spherical functions;
a = 6378136;            #- the equatorial radius of the earth;
w = 7.292115e-5;        # Rad/s is the angular velocity of the earth's rotation;
n = 6

def length(x:list)->float:
    return sqrt(x[0]**2 + x[1]**2 + x[2]**2)

def f(coordinate_before:list)->list:
    r=length(coordinate_before)
    Ak = -mu / (r**3);
    Ck = 1 - 5 * coordinate_before[2] * coordinate_before[2] / (r * r)
    Bk = (-3.0 / 2.0) * J * mu * a * a / pow(r, 5) * Ck
    intermidiate = np.zeros(6)
    for i in range(len(coordinate_before)-3):
        intermidiate[i] = coordinate_before[i+int(len(coordinate_before)/2)]

    intermidiate[3] = Ak * coordinate_before[0] + Bk * coordinate_before[0] + w * w * coordinate_before[0] + 2 * w * coordinate_before[4]
    intermidiate[4] = Ak * coordinate_before[1] + Bk * coordinate_before[1] + w * w * coordinate_before[1] - 2 * w * coordinate_before[3]
    intermidiate[5] = Ak * coordinate_before[2] - (3. / 2.) * J * mu * a * a / pow(r, 5) * coordinate_before[2] * (3 - 5 * coordinate_before[2] * coordinate_before[2] / (r * r))
    return intermidiate

def Runge_Kutta(coordinate_before:list, step:int, te:float, ti:float)->list:
    t = lambda te,ti: (ti - te) * 3600
    h=t(te,ti)/step

    k1 = np.array(f(coordinate_before))
    k2 = np.array(f((k1*0.5*h)+coordinate_before))
    k3 = np.array(f((k2*0.5*h)+coordinate_before))
    k4 = np.array(f((k3*h)+coordinate_before))

    return coordinate_before+(k1 + 2*k2 + 2*k3 + k4) * h/6

def  PredictionEphemeris(te:float, ti:float, step:int, coordinate_before:list)->list:

    value = Runge_Kutta(coordinate_before, step, te, ti)

    for i in range(step):
        value = Runge_Kutta(value, step, te, ti)

    return value

if __name__ == '__main__':
    te = 6.25
    ti = 6.5
    step = 100

    coordinate_before = np.array((-14081.752701e+3, 18358.958252e+3,10861.302124e+3,-1.02576358e+3,1.08672147e+3, -3.15732343e+3))

    coordinate_after = PredictionEphemeris(te, ti, step, coordinate_before)
    print(coordinate_after)
