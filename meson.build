project('bardasis_schrieffer', 'cpp',
        default_options : ['cpp_std=c++17',
        'buildtype=release',
        'default_library=shared',
        'b_lto=false',
        'b_asneeded=false'])

cc = meson.get_compiler('cpp')

# Library Dependencies
m_dep = cc.find_library('m', required : false)
eigen = dependency('eigen3')
gsl = dependency('gsl')

#pybind11 stuff
pymod = import('python')
py = pymod.find_installation('python3', required: true)
pybind_dirs = run_command(py, '-c',
      'import pybind11; print(pybind11.get_include())').stdout().strip()

if get_option('buildtype') == 'release'
  add_project_arguments('-fvisibility=hidden', language : 'cpp')
endif

include_dirs = include_directories('include', pybind_dirs)

if meson.get_compiler('cpp').get_id() == 'gcc'
    add_project_arguments('-Wno-int-in-bool-context', language : 'cpp')
    if build_machine.system() == 'darwin'
        add_project_link_arguments('-undefined', 'dynamic_lookup', language : 'cpp')
    endif
endif

if build_machine.system() == 'darwin'
    pydep = declare_dependency(
    include_directories: include_directories([
    py.get_variable('INCLUDEPY'),
    py.get_path('include'),
    py.get_path('platinclude')
    ]))
elif build_machine.system() == 'linux'
    pydep = py.dependency()
endif

subdir('src')

deps = [
    m_dep,
    eigen,
    gsl,
    pydep
    ]


py.extension_module('_bardasis_schrieffer',
    ext_files,
    include_directories: include_dirs,
    dependencies : deps,
    install : true,
    install_dir : join_paths(meson.build_root(), '..', 'python', 'lib', 'bardasis_schrieffer'))

executable('debug', debug_files,
    include_directories: include_dirs,
    dependencies : deps,
    build_by_default:false)
