## Edit Distance in Cython

> Originally forked from https://github.com/honeyext/cdtw
> Extended functionalities include:
> - Edit Distance on Real Sequence (EDR) [4]
> - Sequences can be multi-dimensional

## Installation

```
pip install git+https://github.com/fzyukio/ced
```

[Examples](https://github.com/fzyukio/ced/blob/master/ced/examples.ipynb)

This module implements:

Distance functions:
 * manhattan (DTW)
 * euclidean (DTW)
 * squared euclidean (DTW)
 * binary cost (0 or 1) for EDR

Local constraints(step patterns, step functions):
 * [well known step patterns dp1, dp2, dp3][1]
 * [local constraints classified by Sakoe-Chiba][2]
impo
Global constraints(windows):
 * Itakura parallelogram
 * [Sakoe-chiba band, Palival adjustment window][3]

## Example

### For DTW, provide a query sequence, a reference sequence

```python
import numpy as np
from ced import Settings, Edr, Dtw

# Three dimensional arrays, with the second and third dimensions are dummy
# E.g. the result will be identical to the distance between `r=[4,4,2,4]` and `q=[5,4,5,6,4`
r = np.array([[4, 0, 0], [4, 0, 0], [2, 0, 0], [4, 0, 0]])
q = np.array([[5, 0, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0], [4, 0, 0]])
sigma = np.array([1, 0, 0])

d = Dtw(r, q, settings=Settings(step='p0sym',  # Sakoe-Chiba symmetric step with slope constraint p = 0
                      window='palival',  # type of the window
                      param=2.0,  # window parameter
                      norm=False,  # normalization
                      compute_path=True))
d.get_dist()
#6.0
d.get_cost()
# array([[  1.   1.   2.   4.  inf]
#        [  2.   1.   2.   4.   4.]
#        [  5.   3.   5.   8.   6.]
#        [ inf  inf   5.   7.   6.]])
d.get_path()
#[(0, 0), (0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)]
```

### For EDR, you need to provide the sigmas
 > The recommended values of sigmas are 0.25 * standard deviations (for each dimension), see [4]

```python
import numpy as np
from ced import Settings, Edr, Dtw

# Three dimensional arrays, with the second and third dimensions are dummy
# E.g. the result will be identical to the distance between `r=[4,4,2,4]` and `q=[5,4,5,6,4`
r = np.array([[4, 0, 0], [4, 0, 0], [2, 0, 0], [4, 0, 0]])
q = np.array([[5, 0, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0], [4, 0, 0]])
sigma = np.array([1, 0, 0])

d = Edr(r, q, sigma, settings=Settings(dist='euclid', norm=False, compute_path=True))
d.get_dist()
#1.0
d.get_cost()
# array([[ 0.  0.  0.  1.  1.]
#        [ 0.  0.  0.  1.  1.]
#        [ 1.  1.  1.  1.  2.]
#        [ 1.  1.  1.  2.  2.]])
d.get_path()
#[(0, 0), (0, 1), (1, 2), (2, 3), (3, 4)]
```


[1]: http://cyber.felk.cvut.cz/gerstner/teaching/zbd/dtw.pdf
[2]: http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=1163055&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D1163055
[3]: https://maxwell.ict.griffith.edu.au/spl/publications/papers/sigpro82_kkp_dtw.pdf
[4]: https://cs.uwaterloo.ca/~tozsu/publications/spatial/sigmod05-leichen.pdf
