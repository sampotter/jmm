import logging
import numpy as np
import pickle

np.seterr(all='raise')

logging.basicConfig(level=logging.INFO)

log = logging.getLogger('debug.py')

with open('field.pickle', 'rb') as f:
    field = pickle.load(f)

log.info('solving *PARENT*')
field.parent.solve()
