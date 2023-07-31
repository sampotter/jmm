executable('mesh', 'mesh.c', dependencies : [argp_dep, jmm_dep])
executable('simulate', 'simulate.c', dependencies : [argp_dep, jmm_dep])
executable('render', 'render.c', dependencies : [argp_dep, jmm_dep])
executable('slice', 'slice.c', dependencies : [argp_dep, jmm_dep])

files = ['render_bin_to_img.py', 'simulate.py']

foreach file : files
  fs.copyfile(file, file)
endforeach
