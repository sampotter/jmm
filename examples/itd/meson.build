executable('itd', 'itd.c', dependencies : [jmm_dep, openmp_dep])

files = ['scratch.py']

foreach file : files
  fs.copyfile(file, file)
endforeach
