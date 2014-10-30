FC = pgf90
#FFLAGS = -Mbounds -Minform=inform -g77libs
FFLAGS = -fastsse -Minform=inform -g77libs
#FFLAGS = -fast -Mstandard -g77libs
#FFLAGS = -fast -Mscalarsse -Mvect=sse -Mflushz -Mstandard -g77libs
#FFLAGS = -g -Mstandard -g77libs
LIBS = ${NETCDF}/lib/libnetcdf.a ${GEMLIB}/gemlib.a ${GEMLIB}/cgemlib.a

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $<

FOBJ = dateutil.o diagnostics.o gempak.o maths_utils.o maths.o kind_types.o current_kind.o partition_wind.o registry.o wrf2gem.o wrf2gemsubs.o 

balance: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(FOBJ) $(LIBS)

wrf2gem_parameters.o: wrf2gem_parameters.f90
dateutil.o: dateutil.f90
diagnostics.o: diagnostics.f90
gempak.o: gempak.f90 wrf2gem_parameters.o
maths_utils.o: maths_utils.f90
maths.o: maths.f90 maths_utils.o
kind_types.o: kind_types.f90
current_kind.o: current_kind.f90 kind_types.o
partition_wind.o: partition_wind.f90 current_kind.o
registry.o: registry.f90 diagnostics.o gempak.o wrf2gemsubs.o maths.o partition_wind.o current_kind.o
wrf2gem.o: wrf2gem.f90 registry.o gempak.o wrf2gemsubs.o 
wrf2gemsubs.o: wrf2gemsubs.f90 dateutil.o wrf2gem_parameters.o

clean:
	rm -f *.o *.mod
