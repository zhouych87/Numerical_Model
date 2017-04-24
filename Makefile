# ***********************************************************************
f90comp=mpif90
# ***********************************************************************
solar.out: subroutine.o commondat.mod main.o
	$(f90comp) subroutine.o commondat.o main.o -o solar.out 
commondat.mod:  commondat.f90
	$(f90comp) -c commondat.f90
subroutine.o:  subroutine.f90 commondat.mod
	$(f90comp) -c subroutine.f90
main.o:  main.f90 subroutine.f90 commondat.f90 
	$(f90comp) -c main.f90
# ***********************************************************************
clean:
	rm solar.out *.o *.mod
# ***********************************************************************

