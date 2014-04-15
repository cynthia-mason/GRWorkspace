AM_CFLAGS = -I $(top_srcdir)/include

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	hsl_ma97ss \
	hsl_ma97ds \
	hsl_ma97cs \
	hsl_ma97zs \
	hsl_ma97ss1 \
	hsl_ma97ds1 \
	hsl_ma97cs1 \
	hsl_ma97zs1 \
	hsl_ma97ss2 \
	hsl_ma97ds2 \
	hsl_ma97cs2 \
	hsl_ma97zs2

check_PROGRAMS = \
	hsl_ma97ss \
	hsl_ma97ds \
	hsl_ma97cs \
	hsl_ma97zs \
	hsl_ma97ss1 \
	hsl_ma97ds1 \
	hsl_ma97cs1 \
	hsl_ma97zs1 \
	hsl_ma97ss2 \
	hsl_ma97ds2 \
	hsl_ma97cs2 \
	hsl_ma97zs2

hsl_ma97ss_SOURCES = hsl_ma97ss.c ../../include/hsl_ma97s.h
hsl_ma97ss_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ss.output hsl_ma97ss.data
hsl_ma97ds_SOURCES = hsl_ma97ds.c ../../include/hsl_ma97d.h
hsl_ma97ds_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ds.output hsl_ma97ds.data
hsl_ma97cs_SOURCES = hsl_ma97cs.c ../../include/hsl_ma97c.h
hsl_ma97cs_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97cs.output hsl_ma97cs.data
hsl_ma97zs_SOURCES = hsl_ma97zs.c ../../include/hsl_ma97z.h
hsl_ma97zs_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97zs.output hsl_ma97zs.data
hsl_ma97ss1_SOURCES = hsl_ma97ss1.c ../../include/hsl_ma97s.h
hsl_ma97ss1_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ss1.output hsl_ma97ss1.data
hsl_ma97ds1_SOURCES = hsl_ma97ds1.c ../../include/hsl_ma97d.h
hsl_ma97ds1_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ds1.output hsl_ma97ds1.data
hsl_ma97cs1_SOURCES = hsl_ma97cs1.c ../../include/hsl_ma97c.h
hsl_ma97cs1_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97cs1.output hsl_ma97cs1.data
hsl_ma97zs1_SOURCES = hsl_ma97zs1.c ../../include/hsl_ma97z.h
hsl_ma97zs1_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97zs1.output hsl_ma97zs1.data
hsl_ma97ss2_SOURCES = hsl_ma97ss2.c ../../include/hsl_ma97s.h
hsl_ma97ss2_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ss2.output hsl_ma97ss2.data
hsl_ma97ds2_SOURCES = hsl_ma97ds2.c ../../include/hsl_ma97d.h
hsl_ma97ds2_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97ds2.output hsl_ma97ds2.data
hsl_ma97cs2_SOURCES = hsl_ma97cs2.c ../../include/hsl_ma97c.h
hsl_ma97cs2_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97cs2.output hsl_ma97cs2.data
hsl_ma97zs2_SOURCES = hsl_ma97zs2.c ../../include/hsl_ma97z.h
hsl_ma97zs2_LDADD = $(top_builddir)/src/libhsl_ma97.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma97zs2.output hsl_ma97zs2.data
