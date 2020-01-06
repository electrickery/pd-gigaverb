# Makefile to build class 'helloworld' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = gigaverb~

# input source file (class name == source file basename)
class.sources := gigaverb~.c 

suppress-wunused := yes

# all extra files to be included in binary distribution of the library
datafiles = gigaverb~-help.pd gigaverb~-meta.pd output~.pd voice.wav

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
include Makefile.pdlibbuilder
