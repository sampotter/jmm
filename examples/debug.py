import numpy as np
import pickle

from jmm.multiple_arrivals import ReflectedField
np.seterr(all='raise')

with open('reflection_BCs.pickle', 'rb') as f:
    domain, bd_faces, bd_T, bd_grad_T = pickle.load(f)

# field = ReflectedField(domain, bd_faces, bd_T, bd_grad_T)
