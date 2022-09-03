executable('itd', 'itd.c', dependencies : jmm_dep)

files = ['HUTUB_pp2_in_cube_50k_high_res_ear.off',
         'HUTUB_pp2_in_cube_50k.off',
         'HUTUB_pp2_in_cube_50k_origin_fixed.off',
         'HUTUB_pp2_in_cube_50k_origin_fixed_bust_only_BROKEN.off',
         'scratch.py',
         'fliege64.txt']

foreach file : files
  custom_target(
    file, input : file, output : file,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true)
endforeach
