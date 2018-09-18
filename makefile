######  Fortran Complier  ##################################
#
#
#

objects = \
kappa.o pmf.o init.o cpmdrate.o

cpmdrate : $(objects)
	${FC} ${FLAGS} -o cpmdrate $(objects)

$(objects): %.o : %.F90
	${FC} ${FLAGS} -c $<

#
#
############################################################
