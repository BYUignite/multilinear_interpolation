import numpy as np
from   scipy.interpolate         import RegularGridInterpolator

x = np.array([1.,2,5,9])
y = np.array([1.,2,5,9,10])

f = np.empty((len(x),len(y)))

for i in range(len(x)):
    for j in range(len(y)):
        f[i,j] = np.sqrt(i+j)

interp = RegularGridInterpolator((x,y), f)
    
xP = 7.0
yP = 2.2

print(f'{interp(np.array([xP,yP]))[0]:.5f}')

