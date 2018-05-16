cdata.set('TOP_SRCDIR', meson.source_root())
cdata.set('TOP_BUILDDIR', meson.build_root())

doxygen = find_program('doxygen', required : false)

doxyfile = configure_file(input: 'Doxyfile.in',
                          output: 'Doxyfile',
                          configuration: cdata,
                          install: false)

if doxygen.found()
  doxygen_docs = custom_target('doxygen_docs',
      input: doxyfile,
      output: ['xml'],
      depend_files : includes,
      command: [doxygen, doxyfile],
      install: true,
      install_dir: join_paths(meson.build_root(), '..', 'python', 'doc'))
endif
