executable(
  'simulate',
  'simulate.c',
  dependencies : jmm_dep
)

custom_target('room.off',
  input : 'room.off',
  output : 'room.off',
  command : ['cp', '@INPUT@', '@OUTPUT@'],
  install : false,
  build_by_default : true)
