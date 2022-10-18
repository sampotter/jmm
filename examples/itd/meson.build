executable('itd', 'itd.c', dependencies : jmm_dep)

files = ['HUTUB_pp2_in_cube_50k.off', 'scratch.py']

foreach file : files
  custom_target(
    file, input : file, output : file,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true)
endforeach
