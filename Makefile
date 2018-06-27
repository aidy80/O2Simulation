COMPILER = gfortran
FLAGS = -O3

OBJS = constants.o mdSubroutines.o analysisSubroutines.o

nveSimO2: NVEO2.f90 $(OBJS)
	$(COMPILER) $(FLAGS) -o nveSimO2 NVEO2.f90 $(OBJS)

nvtSimO2: NVTO2.f90 testMod.o $(OBJS)
	$(COMPILER) $(FLAGS) -o nvtSimO2 NVTO2.f90 testMod.o $(OBJS)

O2Box: O2Box.f90
	$(COMPILER) $(FLAGS) -o O2Box O2Box.f90

drawVelDist: drawVelDist.f90
	$(COMPILER) $(FLAGS) -o drawVelDist drawVelDist.f90

constants.o: constants.f90
	$(COMPILER) $(FLAGS) -c constants.f90

mdSubroutines.o: mdSubroutines.f90 constants.o
	$(COMPILER) $(FLAGS) -c mdSubroutines.f90

analysisSubroutines.o: analysisSubroutines.f90 constants.o
	$(COMPILER) $(FLAGS) -c analysisSubroutines.f90

testMod.o: testMod.f90 constants.o
	$(COMPILER) $(FLAGS) -c testMod.f90

#sim: sim.f90
#	gfortran -o sim sim.f90
