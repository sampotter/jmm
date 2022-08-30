executable('itd', 'itd.c', dependencies : jmm_dep)

custom_target('HUTUB_pp2_in_cube_50k.off',
  input : 'HUTUB_pp2_in_cube_50k.off',
  output : 'HUTUB_pp2_in_cube_50k.off',
  command : ['cp', '@INPUT@', '@OUTPUT@'],
  install : false,
  build_by_default : true)
