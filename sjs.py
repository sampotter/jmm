'''
TODO: add some documentation
'''

from _jmm import *

from _jmm import _ind2l
from _jmm import _ind2lc
from _jmm import _indc2l
from _jmm import _indc2lc
from _jmm import _l2ind
from _jmm import _l2indc
from _jmm import _lc2ind
from _jmm import _lc2indc
from _jmm import _l2lc
from _jmm import _lc2l
from _jmm import _xy_to_lc_and_cc

def get_constant_slowness_field2():
    '''Get a Field2 instance corresponding to a the slowness function $s
\equiv = 1$.

    '''
    return Field2(lambda x, y: 1.0, lambda x, y: (0.0, 0.0))

def get_linear_speed_field2(vx, vy):
    '''Get a Field2 instance corresponding to a slowness function
representing the linear speed function $c(x, y) = 1 + \texttt{vx}
\cdot x + \texttt{vy} \cdot y$ (note that $s \equiv 1/c$.

    '''
    s = lambda x, y: 1.0/(1.0 + vx*x + vy*y)
    return Field2(s, lambda x, y: (-vx*s(x, y)**2, -vy*s(x, y)**2))
