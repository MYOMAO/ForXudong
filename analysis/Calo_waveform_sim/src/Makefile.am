AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libcalowaveformsim.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

pkginclude_HEADERS = \
  CaloWaveFormSim.h

if ! MAKEROOT6
  ROOT5_DICTS = \
    CaloWaveFormSim_Dict.cc
endif

libcalowaveformsim_la_SOURCES = \
  $(ROOT5_DICTS) \
  CaloWaveFormSim.cc 

libcalowaveformsim_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lfun4all \
  -lg4detectors_io \
  -lg4detectors \
  -lsph_onnx \
  -lcalo_reco \
  -lphg4hit


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libcalowaveformsim.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
