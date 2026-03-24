FC = gfortran

# Detect OS
ifeq ($(OS),Windows_NT)
    RM = del /Q
    EXE_EXT = .exe
else
    RM = rm -f
    EXE_EXT =
endif

TARGET = waves$(EXE_EXT)

# Debug flags
FFLAGS_DEBUG = -g -O0 -fcheck=all -Wall \
  -fno-range-check -fno-automatic -ffixed-line-length-none

# Release flags
FFLAGS_RELEASE = -O3 \
  -fno-range-check -fno-automatic -ffixed-line-length-none

# Default
FFLAGS = $(FFLAGS_DEBUG)

OBJS = Waves.o WavesAux.o WavesPrint.o IRI.o \
       lapack.o lapack1.o lapack2.o lapack3.o fft.o \
       nrlmsise00_modified.o hwm93.o hwm07e_modified.o hwm14.o

.PHONY: all clean debug release

all: $(TARGET)

debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean $(TARGET)

release: FFLAGS = $(FFLAGS_RELEASE)
release: clean $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Fixed-form Fortran
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Free-form Fortran
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
ifeq ($(OS),Windows_NT)
	-@$(RM) *.exe *.x *.o *.mod 2>nul
else
	-@$(RM) *.exe *.x *.o *.mod
endif