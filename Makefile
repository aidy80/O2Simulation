nveSimO2: NVEO2.f90
	gfortran -o nveSimO2 NVEO2.f90

nvtSimO2: NVTO2.f90
	gfortran -o nvtSimO2 NVTO2.f90
	
sim: sim.f90
	gfortran -o sim sim.f90

O2Box: O2Box.f90
	gfortran -o O2Box O2Box.f90
