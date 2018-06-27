!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 50,000 4fs timesteps in an 
!NVE ensemble

program nveSim
    use Constants
    use mdSubroutines
    use analysisSubroutines

    implicit none

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]

    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel!oldPos: [m]
                                                               !oldVel: [m/s] 

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 

    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure

    real(dp) :: potentialEnergy ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system


    !Used to time simulation
    real :: start_time
    real :: end_time

    real(dp), dimension(CvvStep) :: Cvv
    real(dp), dimension(CvvStep) :: CvvRot
    real(dp), dimension(MSDStep) :: MSD
    real(dp), dimension(:, :, :), allocatable :: vStore
    real(dp), dimension(:, :, :), allocatable :: pStore
    real(dp), dimension(:, :, :), allocatable :: orientStore
    integer :: nCorr = 0!Used to count the number of times Cvv is added to 
    integer :: nCorr_MSD = 0!Used to count the number of times MSD is added to 
    integer :: nCorr_Rot = 0!Used to count the number of times CvvRot is added to 

    !Iterators
    integer :: p

    real(dp) :: cutoff ![m] Cutoff distance for short-range interaction
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared used for optimization

    !Set up data for finding g(r)
    integer :: numBins 
    integer :: numFullBins 
    integer, dimension(:), allocatable :: bins

    !Speed Distribution
    integer, dimension(numVelBox) :: velBoxes

    !Garbage
    character :: nullChar

    CALL cpu_time(start_time)

    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldPos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldVel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(force(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vStore(numMolecules, numDimensions, CvvStep))
    allocate(pStore(numMolecules, numDimensions, MSDStep))
    allocate(orientStore(numMolecules, numDimensions, CvvStep))

    call openNVEFiles()
    call readIn(pos, vel, dim)

    !Set up information about g(r)
    numBins = ANint((dim(1)*1.733)/(2.0 * delR)) + 1
    numFullBins = ANint(dim(1)/(2.0*delR))

    !Set the cutoff for LJ interation to 1/2 box size
    cutoff = dim(1)/2.0
    cutoffSq = cutoff**2

    allocate(bins(numBins))

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)

    !Init analysis arrays
    call initAnalysisArrays(Cvv, CvvRot, MSD, bins, numBins, velBoxes)

    do p = 1, numNVESteps
        if (mod(p,100) == 0) then
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
        call calcForces(potentialEnergy, pos, force, bins, numBins, dim, cutoffSq, "E")

        !Update positions based on forces
        call leapFrogAndShake(kineticEnergy, force, vel, pos, oldPos, oldVel)

        !Find the kinetic energy of the system and 
        !the total energy of the system
        totEnergy = potentialEnergy + kineticEnergy

        !Call analysis functions to calculate TCF for velocity, 
        !rotation, and distance
        call TCF(p, Cvv, nCorr, vStore, vel)
        call Calc_MSD(p, MSD, nCorr_MSD, pStore, pos)
        call TCFRot(p, CvvRot, nCorr_Rot, orientStore, pos)
        call calcSpeedDist(vel, velBoxes)

        call writeEnergy(p, kineticEnergy, potentialEnergy, totEnergy)
        call writeNVETrajectory(p, dim(1), pos, vel)
    end do

    call TCF_Analysis(Cvv, nCorr, CvvRot, nCorr_Rot, MSD, nCorr_MSD) 
    call Bins_Analysis(velBoxes, bins, numBins, numFullBins, dim(1))

    CALL cpu_time(end_time)
    print *, "Time usage:", (end_time - start_time)/60.0, " minutes", &
                            &mod((end_time - start_time),60.0), " seconds"

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(oldPos)
    deallocate(oldVel)
    deallocate(force)
    deallocate(vStore)
    deallocate(pStore)
    deallocate(orientStore)
    deallocate(bins)

    call closeNVEFiles()

end program nveSim
