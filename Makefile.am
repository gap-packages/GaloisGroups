#
# GaloisGroups: Computing Galois Groups using GAP and PARI
#
# This file is part of the build system of a GAP kernel extension.
# Requires automake.
#
ACLOCAL_AMFLAGS = -I cnf/m4

BINARCHDIR = bin/$(GAPARCH)
GAPINSTALLLIB = $(BINARCHDIR)/GaloisGroups.so

lib_LTLIBRARIES = GaloisGroups.la

GaloisGroups_la_SOURCES = src/GaloisGroups.c
GaloisGroups_la_CPPFLAGS = $(GAP_CPPFLAGS) -DCONFIG_H
GaloisGroups_la_LDFLAGS = -module -avoid-version
if SYS_IS_CYGWIN
GaloisGroups_la_LDFLAGS += -no-undefined -version-info 0:0:0 
GaloisGroups_la_LDFLAGS += -Wl,$(GAPROOT)/bin/$(GAPARCH)/gap.dll
endif

all-local: $(GAPINSTALLLIB)

$(GAPINSTALLLIB): GaloisGroups.la
	$(mkdir_p) $(BINARCHDIR)
if SYS_IS_CYGWIN
	cp .libs/GaloisGroups.dll $(GAPINSTALLLIB)
else
	cp .libs/GaloisGroups.so $(GAPINSTALLLIB)
endif

clean-local:
	rm -f $(GAPINSTALLLIB)

distclean-local:
	rm -rf bin/*

doc:
	($(GAPROOT)/bin/gap.sh -A makedoc.g)
