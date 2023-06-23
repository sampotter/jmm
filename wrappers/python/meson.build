python = import('python').find_installation()

python_dep = python.dependency()

python.extension_module(
  'jmm',
  'jmm.pyx',
  dependencies : [python_dep, jmm_dep]
)
