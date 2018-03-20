noinst_PROGRAMS = templatelapacktest templatelapacktest_threaded

noinst_LIBRARIES = libtemplatelapack.a

templatelapacktest_SOURCES = template_lapack_test.cc
templatelapacktest_LDADD   = 				\
	libtemplatelapack.a $(FLIBS)			\
	$(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a

templatelapacktest_threaded_SOURCES = template_lapack_test_threaded.cc
templatelapacktest_threaded_LDADD   = 				\
	libtemplatelapack.a $(FLIBS)			\
	$(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a


libtemplatelapack_a_SOURCES = 		\
template_lapack_common.cc		\
template_lapack_common.h		\
template_lapack_geqr2.h			\
template_lapack_geqrf.h			\
template_lapack_gesv.h			\
template_lapack_getf2.h			\
template_lapack_getrf.h			\
template_lapack_getrs.h			\
template_lapack_ggbak.h			\
template_lapack_ggbal.h			\
template_lapack_ggev.h			\
template_lapack_gghrd.h			\
template_lapack_hgeqz.h			\
template_lapack_labad.h			\
template_lapack_lacon.h			\
template_lapack_lacpy.h			\
template_lapack_ladiv.h			\
template_lapack_lae2.h			\
template_lapack_laebz.h			\
template_lapack_laev2.h			\
template_lapack_lag2.h			\
template_lapack_lagtf.h			\
template_lapack_lagts.h			\
template_lapack_laln2.h			\
template_lapack_lamch.h			\
template_lapack_lange.h			\
template_lapack_lanhs.h			\
template_lapack_lanst.h			\
template_lapack_lansy.h			\
template_lapack_lapy2.h			\
template_lapack_lapy3.h			\
template_lapack_larfb.h			\
template_lapack_larfg.h			\
template_lapack_larf.h			\
template_lapack_larft.h			\
template_lapack_larnv.h			\
template_lapack_lartg.h			\
template_lapack_laruv.h			\
template_lapack_lascl.h			\
template_lapack_laset.h			\
template_lapack_lasr.h			\
template_lapack_lasrt.h			\
template_lapack_lasv2.h			\
template_lapack_laswp.h			\
template_lapack_latrd.h			\
template_lapack_latrs.h			\
template_lapack_org2l.h			\
template_lapack_org2r.h			\
template_lapack_orgql.h			\
template_lapack_orgqr.h			\
template_lapack_orgtr.h			\
template_lapack_orm2r.h			\
template_lapack_ormqr.h			\
template_lapack_pocon.h			\
template_lapack_potf2.h			\
template_lapack_potrf.h			\
template_lapack_pptrf.h			\
template_lapack_rscl.h			\
template_lapack_spgst.h			\
template_lapack_stebz.h			\
template_lapack_stein.h			\
template_lapack_steqr.h			\
template_lapack_sterf.h			\
template_lapack_stevx.h			\
template_lapack_larra.h                 \
template_lapack_larrb.h                 \
template_lapack_larrc.h                 \
template_lapack_larrd.h                 \
template_lapack_larre.h                 \
template_lapack_larrf.h                 \
template_lapack_larrj.h                 \
template_lapack_larrk.h                 \
template_lapack_larrr.h                 \
template_lapack_larrv.h                 \
template_lapack_lar1v.h                 \
template_lapack_laneg.h                 \
template_lapack_isnan.h                 \
template_lapack_laisnan.h               \
template_lapack_lasq2.h			\
template_lapack_lasq3.h			\
template_lapack_lasq4.h			\
template_lapack_lasq5.h			\
template_lapack_lasq6.h			\
template_lapack_stemr.h			\
template_lapack_stevr.h			\
template_lapack_syev.h			\
template_lapack_sygs2.h			\
template_lapack_sygst.h			\
template_lapack_sygv.h			\
template_lapack_sytd2.h			\
template_lapack_sytrd.h			\
template_lapack_test.cc 		\
template_lapack_tgevc.h			\
template_lapack_tptri.h			\
template_lapack_trti2.h			\
template_lapack_trtri.h

EXTRA_DIST =  \
test.sh

AM_CPPFLAGS = \
	-I$(top_srcdir)/source/matrix/template_lapack/blas

TESTS_ENVIRONMENT = top_builddir=$(top_builddir) top_srcdir=$(top_srcdir)

TESTS = \
	test.sh


