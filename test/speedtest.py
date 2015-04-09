import numpy as np
import time

from cdtw.src.pydtw import*

r = np.arange(6000)
q = np.arange(6000)

s = Settings()
s.step.set_type('dp2')

# s.global_constraint.set_type('itakura')
# s.global_constraint.set_param(0.2)
# 

s.compute_path = False


t1 = time.time()
d = dtw(r, q, s)
t2 = time.time()

print s
print "python total time: " + str(t2 - t1)


import mlpy
t1 = time.time()
d = mlpy.dtw_std(r, q, dist_only = False)
t2 = time.time()
print "mlpy total time: " + str(t2 - t1)