##Dynamic Time Warping Project

[Examples](http://nbviewer.ipython.org/github/honeyext/cdtw/blob/master/examples.ipynb)

This module implements:

Distance functions:
 * manhattan 
 * euclidean 
 * squared euclidean

Local constraints(step patterns, step functions):
 * [well known step patterns dp1, dp2, dp3][1]
 * [local constraints classified by Sakoe-Chiba][2]
impo
Global constraints(windows):
 * Itakura parallelogram
 * [Sakoe-chiba band, Palival adjustment window][3]
 
```python
import numpy as np
from cdtw import pydtw
r = np.array([1,2,3,4])
q = np.array([2,3,4,5])
d = pydtw.dtw(r,q,pydtw.Settings(step = 'p0sym',     #Sakoe-Chiba symmetric step with slope constraint p = 0
                                window = 'palival', #type of the window
                                param = 2.0,        #window parameter
                                norm = False,       #normalization
                                compute_path = True))

d.get_dist()
#2.0
d.get_cost()
#array([[  1.,   3.,   6.,  inf],
#       [  1.,   2.,   4.,   7.],
#       [  2.,   1.,   2.,   4.],
#       [ inf,   2.,   1.,   2.]])
d.get_path()
#[(0, 0), (1, 0), (2, 1), (3, 2), (3, 3)]

  


```

[1]: http://cyber.felk.cvut.cz/gerstner/teaching/zbd/dtw.pdf
[2]: http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=1163055&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D1163055
[3]: https://maxwell.ict.griffith.edu.au/spl/publications/papers/sigpro82_kkp_dtw.pdf
