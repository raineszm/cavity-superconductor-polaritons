project('bardasis_schrieffer', 'cpp',
        default_options : ['cpp_std=c++17',
        'buildtype=release',
        'default_library=static',
        'b_lto=true'])

cc = meson.get_compiler('cpp')

# Library Dependencies
m_dep = cc.find_library('m', required : false)
#fmt = cc.find_library('fmt')
eigen = dependency('eigen3')
nlopt = dependency('nlopt')
boost = dependency('boost', modules : ['math'])
gsl = dependency('gsl')

#pybind11 stuff
py3 = import('python3')
pybind_dirs = run_command('python', '-c',
      'import pybind11; print(pybind11.get_include())').stdout().strip()

if get_option('buildtype') == 'release'
  add_project_arguments('-fvisibility=hidden', language : 'cpp')
endif

include_dirs = include_directories('include', pybind_dirs)

add_project_arguments('-ffast-math', language : 'cpp')

subdir('src')

deps = [
    m_dep,
    nlopt,
    boost,
    eigen,
    gsl
    ]

py3.extension_module('bardasis_schrieffer',
    ext_files,
    include_directories: include_dirs,
    dependencies : deps,
    install : true,
    install_dir : join_paths(meson.build_root(), '..', 'python', 'lib'),
    cpp_args : '-DDOCTEST_CONFIG_DISABLE')

subdir('include')
cdata = configuration_data()
subdir('doc')