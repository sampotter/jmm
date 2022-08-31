executable('mesh', 'mesh.c', dependencies : jmm_dep)
executable('simulate', 'simulate.c', dependencies : jmm_dep)
executable('render', 'render.c', dependencies : jmm_dep)
executable('slice', 'slice.c', dependencies : jmm_dep)

custom_target('room.off',
  input : 'room.off',
  output : 'room.off',
  command : ['cp', '@INPUT@', '@OUTPUT@'],
  install : false,
  build_by_default : true)

custom_target('render_bin_to_img.py',
  input : 'render_bin_to_img.py',
  output : 'render_bin_to_img.py',
  command : ['cp', '@INPUT@', '@OUTPUT@'],
  install : false,
  build_by_default : true)
