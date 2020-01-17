import numpy as np
from   scipy.interpolate import interp1d

x = np.array([1.,2,5,9])
f = np.array([1.,2,3,4])

interp = interp1d(x, f)
    
xP = 7.0

print(interp(xP))

