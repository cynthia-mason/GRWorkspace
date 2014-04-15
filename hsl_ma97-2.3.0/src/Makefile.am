lib_LIBRARIES = libhsl_ma97.a libfakemetis.a

modulesdir = modules
nodist_modules_HEADERS = \
	hsl_ma97_single.mod \
	hsl_ma97_double.mod \
	hsl_ma97_complex.mod \
	hsl_ma97_double_complex.mod

libhsl_ma97_a_SOURCES = \
	hsl_ma97s.f90 \
	hsl_ma97d.f90 \
	hsl_ma97c.f90 \
	hsl_ma97z.f90 \
	common.f \
	common90.f90 \
	sdeps90.f90 \
	ddeps90.f90 \
	cdeps.f \
	cdeps90.f90 \
	zdeps.f \
	zdeps90.f90 \
	hsl_ma97s_ciface.f90 \
	hsl_ma97d_ciface.f90 \
	hsl_ma97c_ciface.f90 \
	hsl_ma97z_ciface.f90

libfakemetis_a_SOURCES = fakemetis.f

hsl_ma97s.$(OBJEXT): sdeps90.$(OBJEXT)
hsl_ma97s.$(OBJEXT): common90.$(OBJEXT)
hsl_ma97d.$(OBJEXT): ddeps90.$(OBJEXT)
hsl_ma97d.$(OBJEXT): common90.$(OBJEXT)
hsl_ma97c.$(OBJEXT): cdeps90.$(OBJEXT)
hsl_ma97c.$(OBJEXT): common90.$(OBJEXT)
hsl_ma97z.$(OBJEXT): zdeps90.$(OBJEXT)
hsl_ma97z.$(OBJEXT): common90.$(OBJEXT)
sdeps90.$(OBJEXT): common90.$(OBJEXT)
ddeps90.$(OBJEXT): common90.$(OBJEXT)
cdeps90.$(OBJEXT): common90.$(OBJEXT)
zdeps90.$(OBJEXT): common90.$(OBJEXT)
hsl_ma97s_ciface.$(OBJEXT): hsl_ma97s.$(OBJEXT)
hsl_ma97d_ciface.$(OBJEXT): hsl_ma97d.$(OBJEXT)
hsl_ma97c_ciface.$(OBJEXT): hsl_ma97c.$(OBJEXT)
hsl_ma97z_ciface.$(OBJEXT): hsl_ma97z.$(OBJEXT)
