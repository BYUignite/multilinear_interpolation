import numpy as np
from   scipy.interpolate         import RegularGridInterpolator

x = np.array([1.,2,5,9])
y = np.array([1.,2,5,9,10])
z = np.array([1.,2,5,9,10,11])

f = np.empty((len(x),len(y),len(z)))

for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            f[i,j,k] = np.sqrt(i+j+k)

interp = RegularGridInterpolator((x,y,z), f)
    
xP = 7.0
yP = 2.2
zP = 5.5

print(f'{interp(np.array([xP,yP,zP]))[0]:.5f}')

