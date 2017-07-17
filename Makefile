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
else ifeq (${USER}, julien)
	FF = gfortran
	FPY = f2py-2.7
	OPT = --opt=-ffixed-line-length-none --opt=-O3
else
	FF = gfortran
	FPY = f2py
	OPT = --opt=-ffixed-line-length-none --opt=-O3
endif

all: cmb

cmb:
	${FPY} -c s4cmb/scanning_strategy_f.f90 -m scanning_strategy_f ${OPT}
	${FPY} -c s4cmb/detector_pointing_f.f90 -m detector_pointing_f ${OPT}
	${FPY} -c s4cmb/tod_f.f90 -m tod_f ${OPT}
	${FPY} -c s4cmb/systematics_f.f90 -m systematics_f ${OPT}
	-mv *.so s4cmb/

clean:
	-rm s4cmb/*.o s4cmb/*.so
