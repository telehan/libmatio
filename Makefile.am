##
## Copyright (C) 2005-2006   Christopher C. Hulbert
##
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Lesser General Public
##  License as published by the Free Software Foundation; either
##  version 2.1 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##  Lesser General Public License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License along with this library; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
##
AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
VPATH = @srcdir@

ACLOCAL_AMFLAGS = -Iconfig

datadir = $(prefix)/share/matio
EXTRA_DIST = matio.pc.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = matio.pc
$(pkgconfig_DATA): config.status


if DOCS_INSTALL
DOCS=doxygen
endif

if BUILD_TEST
data_DATA  = share/test_file.mat share/test_file_v6.mat
TEST=test
endif

SUBDIRS = src $(TEST) $(DOCS)

dosubst = sed -e 's,[@]PACKAGE[@],$(PACKAGE),g' \
              -e 's,[@]VERSION[@],$(VERSION),g' \
              -e 's,\/,\\,g'

doc:
	

print:
	subdirs='$(SUBDIRS)'; \
	here=`pwd`;\
	for f in $$subdirs; do \
	    cd $$f; \
	    make print;\
	    cd $$here;\
	done
