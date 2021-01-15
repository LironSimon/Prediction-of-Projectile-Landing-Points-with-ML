#----------------------------------------
# Import Environment
#----------------------------------------
import pandas as pd
import random
import numpy as np
from scipy.integrate import odeint

#----------------------------------------
# Setting projectile values
#----------------------------------------
'''Creating a dictionary with all the features.
Each feature directs to a random array of 500 elements.'''
N = 500
feature_dict = dict(v0 = [random.uniform(0,10**4) for i in range(N)],   # initial velocity [m/s]
                    alpha0 = [random.uniform(0,361) for i in range(N)], # angle of v0 compared to horizon [deg]
                    x0 = [random.uniform(-500,500) for i in range(N)],  # initial horizontal positioning [m]
                    y0 = [random.uniform(0,2000) for i in range(N)],    # initial height relative to ground [m], therefore y0>=0
                    g = [9.81 for i in range(N)],                       # gravitational acceleration [m/s2]
                    m = [random.uniform(0.001,10**4) for i in range(N)],# projectile's mass [kg]
                    )

#----------------------------------------
# Calculate landing location
#----------------------------------------

# Vx = horizontal velocity component [m/s]
# Vy0 = initial horizontal velocity component [m/s]
# h = max distance from the ground [m]
# t = total time of flight [s]
# x_end = landing point [m]

# Calculate: Vx = V0 * cos(alpha), Vy0 = V0 * sin(alpha)
cos_alpha = np.cos( np.deg2rad(feature_dict['alpha0']) )
sin_alpha = np.sin( np.deg2rad(feature_dict['alpha0']) )
Vx,Vy0 = [],[]

for i,V0 in enumerate(feature_dict['v0']):
    Vx.append(V0 * cos_alpha[i])
    Vy0.append(V0 * sin_alpha[i])

# solve: 0.5g*t^2+Vy0*t+y0=0
t = []
for i in range(N):
    coeff = [0.5*feature_dict['g'][i], Vy0[i], feature_dict['y0'][i]]
    options = sorted(np.roots(coeff))
    if options[-1]<0: t.append(0)
    if options[0]>0: t.append(options[0])

# Calculate: x_end = x0 + Vx * t, and add as target feature to dict
'''CHECK! Is it a problem that there's no two targets (x_end) alike?'''
x_end, mult = [], np.multiply(Vx, t).tolist()
for i in range(N):
    x_end.append( feature_dict['x0'][i] * mult[i])

feature_dict["x_end"] = x_end

df = pd.DataFrame(data = feature_dict)

print(df)




