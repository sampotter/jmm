executable('curvy', ['curvy.c'], dependencies : [jmm_dep])

scripts = [
  'generate_BaseSnowTest.py',
]

foreach script : scripts
  custom_target(
    script,
    input : script, output : script,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true
  )
endforeach
