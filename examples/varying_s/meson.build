executable('mesh', 'mesh.c', dependencies : [argp_dep, jmm_dep])
executable('varying_s', 'varying_s.c', dependencies : [argp_dep, jmm_dep])
executable('varying_s_convergence_test', 'varying_s_convergence_test.cpp', dependencies : [argp_dep, jmm_dep, tetgen_dep])

files = ['make_linear_speed_func.py']

foreach file : files
  fs.copyfile(file, file)
endforeach
