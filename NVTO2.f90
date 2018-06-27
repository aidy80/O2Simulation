!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 100,000 4fs timesteps in a NVT
!ensemble. Temperature will begin at 320K and be gradually scaled to 77K

program nvtSim
    use Constants
    use mdSubroutines
    use analysisSubroutines
    use testMod

    implicit none

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]
    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel !oldPos: [m]
                                                                !oldVel: [m/s]
    real(dp), dimension(numDimensions) :: avgPos !Temp variable to measure 
                                                 !pos of an O2 molecule
    integer, parameter :: numBins = 10
    integer, dimension(numBins) :: bins

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 
    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure

    real(dp) :: cutoff ![m] Cutoff distance for short-range interaction
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared used for optimization

    real(dp) :: potentialEnergy ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system

    !Used to time simulation
    real :: start_time
    real :: end_time

    !Interator
    integer :: p

    real(dp) :: desiredTemperature = initTemperature

    CALL cpu_time(start_time)

    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldPos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldVel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(force(numMolecules, numAtomsPerMolecule, numDimensions))

    call openNVTFiles()
    call readIn(pos, vel, dim)

    !Set the cutoff for LJ interation to 1/2 box size
    cutoff = dim(1)/2.0
    cutoffSq = cutoff**2

    close (unit=11)

    call testLJ(sigma, epsilon)

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)

    do p = 1, numNVTSteps
        if(mod(p,100) == 0) then
            print *, p
        end if
        !Init force vectors to zero
        call initZero(force)

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        if (mod(p, zeroMomentTimeStep) == 0) then
            call zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)
        end if

        potentialEnergy = 0.0
        kineticEnergy = 0.0

        !Calculate forces at the current timestep
        call calcForces(potentialEnergy, pos, force, bins, numBins, dim, cutoffSq, "T")

        !Update positions based on forces
        call leapFrogAndShake(kineticEnergy, force, vel, pos, oldPos, oldVel)
        
        !Find the kinetic energy of the system and 
        !the total energy of the system
        totEnergy = potentialEnergy + kineticEnergy

        !Output the energy at the current timestep
        call writeEnergy(p, kineticEnergy, potentialEnergy, totEnergy)

        !Gradually scale down the temperature until it reaches a 
        !constant Temperature of 77K
        if (p.LE.(numNVTSteps - numNVTConstSteps*2).AND.p.GE.numNVTConstSteps&
                                        &.AND.mod(p, temperatureStep) == 0) then
            print *, "temp: ", desiredTemperature
            desiredTemperature = desiredTemperature - 1
        end if

        !Scale the temperature of the system down to the desired value
        call scaleTemp(desiredTemperature, kineticEnergy, vel, numMolecules, &
                                     &numAtomsPerMolecule, numDimensions, Bolz)

        !Write vel and position information every 100 timesteps
        call writeNVTTrajectory(p, dim(1), pos, vel)
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(oldPos)
    deallocate(oldVel)
    deallocate(force)

    call closeNVTFiles()

    CALL cpu_time(end_time)
    print *, "Time usage:", (end_time - start_time)/60.0, " minutes", &
                            &mod((end_time - start_time),60.0), " seconds"

end program nvtSim
