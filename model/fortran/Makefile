#FFLAGS := -fcheck=bounds -ffree-line-length-0
FFLAGS := -fcheck=bounds -ffree-line-length-0 -L /usr/x86_64-w64-mingw32/lib -static
F_SOURCES := $(wildcard *.F95)
F_EXECUTABLE := $(F_SOURCES:.F95=.exe)
FFC := x86_64-w64-mingw32-gfortran

all:$(F_EXECUTABLE)

$(F_EXECUTABLE):$(F_SOURCES)
	$(FFC) $< $(FFLAGS) -o ../compiled/$@

