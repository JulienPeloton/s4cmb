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

all: scanning

scanning:
	${FPY} -c scanning_strategy_f.f90 -m scanning_strategy_f ${OPT}

clean:
	-rm *.o *.so
