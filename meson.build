project('QP', 'c')

cc = meson.get_compiler('c')

py_mod = import('python')
py = py_mod.find_installation('/home/rootacess/qp/.venv/bin/python3')
py_dep = py.dependency(required: true)

incdir = include_directories('/usr/include/python3.10')
libdir = '/usr/lib/x86_64-linux-gnu'

main = executable('main', 'main.c',
           dependencies: py_dep,
           include_directories: incdir,
           link_args: ['-L' + libdir, '-lpython3.10'])

test('basic', main)