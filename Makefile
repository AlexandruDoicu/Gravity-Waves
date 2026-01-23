gf = gfortran -fno-range-check -fno-automatic -ffixed-line-length-none

#
compile:
	$(gf) -c Waves.f
	$(gf) -c WavesAux.f
	$(gf) -c WavesPrint.f
	$(gf) -c IRI.f
	$(gf) -c lapack.f
	$(gf) -c lapack1.f
	$(gf) -c lapack2.f
	$(gf) -c lapack3.f
	$(gf) -c fft.f
	$(gf) -c nrlmsise00_modified.f
	$(gf) -c hwm93.f
	$(gf) -c hwm07e_modified.f90
	$(gf) -c hwm14.f90
	$(gf) -o waves.x Waves.o WavesAux.o WavesPrint.o IRI.o\
		 lapack.o lapack1.o lapack2.o lapack3.o fft.o\
		 nrlmsise00_modified.o hwm93.o hwm07e_modified.o  hwm14.o

#
clean:
	rm *.x *.o *.mod
