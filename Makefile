## use ifort and f2py at NERSC
ifeq (${NERSC_HOST}, edison)
	FF = ifort
	FPY = f2py
	OPT = --opt=-O3 -lifcore
else ifeq (${NERSC_HOST}, cori)
	FF = ifort
	FF = gfortran
	FPY = f2py
	OPT = --opt=-O3
else
	FF = gfortran
	FPY = f2py-2.7
	OPT = --opt=-ffixed-line-length-none --opt=-O3
endif

all: s4cmb

s4cmb:
	${FPY} -c instrument/scanning_strategy_f.f90 -m scanning_strategy_f ${OPT}
	-mv scanning_strategy_f.so instrument/
	${FPY} -c time_ordered_data/detector_pointing_f.f90 -m detector_pointing_f ${OPT}
	-mv detector_pointing_f.so time_ordered_data/
	${FPY} -c time_ordered_data/tod_f.f90 -m tod_f ${OPT}
	-mv tod_f.so time_ordered_data/

clean:
	-rm *.o *.so
