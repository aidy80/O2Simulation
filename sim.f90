!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 150,000 4fs timesteps in an 
!NVT ensemble which scales from 320K-77K then enters an NVE ensemble

program nveSim
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)!Give reals double precision

    !Iterators
    integer :: i, j, k, m, l, n, p

    !Constants
    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space
    real(dp), parameter :: oMass = 2.66E-26 ![kg] Mass of an oxygen atom
    real(dp),parameter :: invOMass = 1.0 / oMass ![m^(-1)]
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: epsilon = 48.0 * Bolz![J] Minimum in the 
                                                 !    Lennard Jones equation
    real(dp), parameter :: desiredBondLength = 1.208E-10 ![m]Bond length betweem 
                                                         !oxygen atoms
    real(dp), parameter :: desiredBondLengthSq = desiredBondLength**2![m^2]
    real(dp), parameter :: allowedError = 1.0E-13 ![]Allowed error between actual 
                                                  !bond length and desired
    real(dp), parameter :: sigma = 3.006E-10 ![m]. Constant in the Lennard 
                                            !     Jones equation
    real(dp), parameter :: sigmaSq = sigma**2 ![m^2] Sigma squared
    real(dp) :: cutoff ![m] Cutoff distance for short-range interaction
    real(dp), parameter :: timeStep = 4.0E-15 ![s] Time between calculations
    real(dp), parameter :: fourEps = 4.0 * epsilon ![J] Epsilon*4. Used for optimization
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared
    real(dp), parameter :: twentyFourEps = 24.0 * epsilon ![J] Epsilon*24. 
                                                          !Used for optimization
    real(dp), parameter :: initTemperature = 320.0 ![K] initial const temperature 
                                                    !of the system
    real(dp), parameter :: finalTemperature = 77.0 ![K] final temperature of 
                                                    !the system
                                               
    integer, parameter :: CvvStep = 500 !Number of items in Cvv
    integer, parameter :: MSDStep = 500 !Number of items in MSD
    integer, parameter :: numTotSteps = 170000 !Number of timesteps in the program
    integer, parameter :: numNVTSteps = 150000 !Number of timesteps in the program
    integer, parameter :: numConstSteps = 50000 !Number of timesteps that the 
                                                !system is equilibrated at a 
                                                !constant temperature
    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                   !between momentum-zeroing
    integer, parameter :: numTrajSteps = 100 !Number of timesteps between 
                                             !trajectory outputs to .ttr file
    integer, parameter :: temperatureStep = 205 !Number of timesteps between 
                                                !temperature decrements
    
    integer, parameter :: numAtomsPerMolecule = 2 !Number of atoms in each Oxygen 
    integer :: numMolecules !Number of molecules in this simulation
    integer :: numAtoms !Number of atoms in the simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]
    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel!oldPos: [m]
                                                               !oldVel: [m/s] 

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 

    !Records distances between different atoms
    real(dp), dimension(:), allocatable :: distBins ![]

    real(dp), dimension(numDimensions) :: pastBondLength ![m]Holds bond length of an 
                                                         !O2 in its previous iteration
    real(dp), dimension(numDimensions) :: currBondLength ![m]Holds current bond length of an O2
    real(dp) :: currBondLengthSq ![m^2] Squared bondLength
    real(dp), dimension(numDimensions) :: avgPos !Used to hold median bond position
    real(dp) :: lambda!Coefficient which must be solved for through iteration
    real(dp) :: currError!Error between current bondLength and desired bondLength

    real(dp) :: Epot ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: potentialSum ![J] the total potential energy of the system
    real(dp) :: oldPotential ![J] the total potential energy of the system
    real(dp) :: vSQ ![(m/s)^2] Square velocity of a given atom
    real(dp) :: vOld ![m/s] Temp variable for velocity
    real(dp) :: Fmag ![N] Used to hold the magnitude of force btwn two atoms
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system at time t
    real(dp) :: kineticEnergyScale ![J] The total kinetic energy of the 
                                   !system at time t+0.5dt
    real(dp) :: addedKE ![J] Sum used to calc KE. Used to avoid     
                        !numerical percision issues
    real(dp) :: addedKEScale ![J] Sum used to calc KEScale. Used to avoid     
                        !numerical percision issues

    real(dp) :: desiredTemperature = initTemperature

    !Used to time simulation
    real :: start_time
    real :: end_time

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance ![m]
    real(dp) :: distanceSq ![m^2]
    real(dp) :: sigmaDistTwo![]. (sigma/(distance between two atoms))**2. 
    real(dp) :: sigmaDistSix ![]. (sigma/(distance between two atoms))**6. 
    real(dp) :: sigmaDistTwelve![]. (sigma/(distance between two atoms))**12.  
    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure

    !Variables to store information needed to calculate TCF, Cvv, TCFRotational
    real(dp), dimension(CvvStep) :: Cvv !Stores the sum of Cvvs calculations
    real(dp), dimension(CvvStep) :: CvvRot !Stores the sum of CvvRot calculations
    real(dp) :: CvvOne = 0 
    real(dp) :: CvvOneRot = 0
    real(dp), dimension(MSDStep) :: MSD
    real(dp), dimension(:, :, :), allocatable :: vStore
    real(dp), dimension(:, :, :), allocatable :: pStore
    real(dp), dimension(:, :, :), allocatable :: orientStore
    integer :: nCorr = 0!Used to count the number of times Cvv is added to 
    integer :: nCorr_MSD = 0!Used to count the number of times MSD is added to 
    integer :: nCorr_Rot = 0!Used to count the number of times CvvRot is added to 
    real(dp) :: diffusionCoeff = 0!Diffusion coefficient to be calcualted

    integer :: convergeCount = 0 !Counts number of iterations
    logical :: hasConverge = .false. !Used in iterative while loop

    !Set up data for finding g(r)
    real(dp), parameter :: delR = 0.003e-9
    integer :: numBins 
    integer :: numFullBins 
    integer, dimension(:), allocatable :: bins
    integer :: currIndex
    real(dp) :: rho 
    real(dp), parameter :: pi = 3.1415
    real(dp) :: sphereConst
    real(dp) :: rLower 
    real(dp) :: rUpper
    real(dp) :: shellVol

    !Garbage variables
    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1
    real(dp) :: dumr


    CALL cpu_time(start_time)

    !Open file containing position and velocity information 
    !about oxygen molecules and several files to save the temperature of
    !the system and analysis tests such as TCF, and MSD
    open(unit=11, file='oxygenInit.gro') 
    open(unit=90, file='Sim_temperature.dat')
    open(unit=91, file='Sim_totEnergy.dat')
    open(unit=92, file='Sim_potentialEnergy.dat')
    open(unit=93, file='Sim_kineticEnergy.dat')
    open(unit=94, file='Sim.gro')
    open(unit=95, file='Sim_final.gro')
    open(unit=96, file='TCFO2Norm.dat')
    open(unit=97, file='grO2.dat')
    open(unit=98, file='MSDO2.dat')
    open(unit=99, file='TCFRotO2.dat')
    open(unit=100, file='TCFO2.dat')
    
    !Read in header information from the file
    read(11, *) nullChar
    read(11, 30) numAtoms
    numMolecules = numAtoms / numAtomsPerMolecule
    
    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldPos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldVel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(force(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vStore(numMolecules, numDimensions, CvvStep))
    allocate(pStore(numMolecules, numDimensions, MSDStep))
    allocate(orientStore(numMolecules, numDimensions, CvvStep))

    !Read in position and velocity information from the .gro file
    do m = 1, numMolecules
        do l = 1, numAtomsPerMolecule
            read(11, 10) nullInt, nullChar, nullChar1, nullInt1, pos(m,l,1), &
                    & pos(m,l,2), pos(m,l,3), vel(m,l,1), vel(m,l,2), vel(m,l,3)
        end do
    end do

    !Read in the dimensions of the box and convert to [m]
    read(11,*) dim(1),dim(2),dim(3)
    do m = 1, numDimensions
        dim(m) = dim(m) * 1.0E-9
    end do

    !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
    do m = 1, numMolecules
        do l = 1, numAtomsPerMolecule
            do k = 1, numDimensions
                pos(m,l,k) = pos(m,l,k) * 1.0E-9
                vel(m,l,k) = vel(m,l,k) * 1.0E3
            end do
        end do
    end do
    
    !Set up information about g(r)
    rho = numAtoms / (dim(1)**3)
    sphereConst = 4.0D0 * pi * rho / 3.0D0
    numBins = ANint((dim(1)*1.733)/(2.0 * delR)) + 1
    numFullBins = ANint(dim(1)/(2.0*delR))

    !Set the cutoff for LJ interation to 1/2 box size
    cutoff = dim(1)/2.0
    cutoffSq = cutoff**2

    allocate(bins(numBins))

    close (unit=11)

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)

    !Initialize Arrays
    do i = 1, CvvStep
        Cvv(i) = 0.0D0
    end do
    
    do i = 1, CvvStep
        CvvRot(i) = 0.0D0
    end do

    do i = 1, MSDStep
        MSD(i) = 0.0D0
    end do

    do i = 1, numBins
        bins(i) = 0
    end do

    do p = 1, numTotSteps
        if (mod(p,500) == 0) then
            print *, p
        end if

        !Init force vectors to zero
        call initZero(force, numMolecules, numAtomsPerMolecule, numDimensions)

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        if (mod(p, zeroMomentTimeStep) == 0) then
            call zeroNetMomentum(vel, numAtoms, numAtomsPerMolecule, numDimensions)
        end if

        Epot = 0.0
        kineticEnergy = 0.0
        kineticEnergyScale = 0.0

        !Compute Force between each pair of atoms if their distance is below
        !the cutoff radius. Each pair is evaluated only once
        do i = 1, numMolecules - 1
            potentialSum = 0.0
            do m = 1, numAtomsPerMolecule
                do j = i + 1, numMolecules
                    do l = 1, numAtomsPerMolecule
                        !Calculate the distance between the current atoms, applying 
                        !periodic boundary conditions to find the minimum possible 
                        !distance for each coordinate
                        do k = 1, numDimensions
                            distance(k) = pos(i, m, k) - pos(j, l, k)
                            distance(k) = distance(k) - dim(k) * &
                                                  &anint(distance(k) / dim(k))
                        end do

                        distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2

                        !If the distance between the two atoms is below the cutoff, 
                        !calculate the force exerted on each of them based on the 
                        !lennard jones potential
                        if (distanceSq < cutoffSq) then
                            sigmaDistTwo = sigmaSq / distanceSq
                            sigmaDistSix = sigmaDistTwo**3
                            sigmaDistTwelve = sigmaDistSix**2

                            !Calc potential from lennard jones equation between the 
                            !current pair of atoms and add it to the current 
                            !potential energy sum
                            potential = fourEps * (sigmaDistTwelve - sigmaDistSix)
                            potentialSum = potentialSum + potential

                            !Calculate the resulting force on the current two atoms 
                            !based on the lennard jones potential between them. Calculated 
                            !using the negative gradient
                            Fmag = twentyFourEps * (2.0 * sigmaDistTwelve - sigmaDistSix)

                            !Apply calculated force the the net force for each atom
                            do k = 1, numDimensions
                                force(i, m, k) = force(i, m, k) + &
                                                    &Fmag * (distance(k) / distanceSq)
                                force(j, l, k) = force(j, l, k) - &
                                                    &Fmag * (distance(k) / distanceSq)
                            end do
                        end if

                        !Add current position information to g(r) data
                        if (p.GE.numNVTSteps) then
                            currIndex = Int(sqrt(distanceSq)/delR) + 1
                            bins(currIndex) = bins(currIndex) + 2
                        end if
                    end do
                end do
            end do
            Epot = Epot + potentialSum
        end do

        !Use the leap-frog verlet algorithm to calculate new position and 
        !velocity vectors for all atoms based on the forces 
        !calculated between them.
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions 
                    oldVel(m, l, k) = vel(m, l, k)
                    oldPos(m, l, k) = pos(m, l, k)
                    vel(m, l, k) = vel(m, l, k) + (force(m, l, k) / oMass) * timeStep
                    pos(m, l, k) = pos(m, l, k) + vel(m, l, k) * timeStep
                end do
            end do

            !SHAKE Alg to fix bond angles
            hasConverge = .false.
            convergeCount = 0
            !Use iteration to solve for lambda, a coefficient adjusting the amount
            !that positions/velocities must be adjusted to prevent bond lengths from 
            !changing
            do while(.NOT.hasConverge)
                do k = 1, numDimensions
                    pastBondLength(k) = oldPos(m, 1, k) - oldPos(m, 2, k)
                    currBondLength(k) = pos(m, 1, k) - pos(m, 2, k)
                end do
               
                lambda = (currBondLength(1)**2 + currBondLength(2)**2 + &
                        &currBondLength(3)**2 - desiredBondLengthSq) / &
                        &(4.0*invOMass*dot(currBondLength, pastBondLength))
                
                !Use lambda to adjust the position and velocities of atoms in the simulation
                currBondLengthSq = 0.0
                do k = 1, numDimensions
                    pos(m, 1, k) = pos(m, 1, k) - pastBondLength(k)*invOMass * lambda
                    pos(m, 2, k) = pos(m, 2, k) + pastBondLength(k)*invOMass * lambda
                    currBondLengthSq = currBondLengthSq + (pos(m, 1, k) - pos(m, 2, k))**2 
                end do

                !Find the error between the current bond 
                !length and desired bond length
                currError = ABS(currBondLengthSq - desiredBondLengthSq) / desiredBondLengthSq

                !Repeat the iteration until the bond is within the allowedError 
                !of the desired length or the iteration occurs 500 times
                if(currError.LE.allowedError) then
                    hasConverge = .true.
                else 
                    hasConverge = .false.
                    convergeCount = convergeCount + 1

                    if (convergeCount > 500) then
                        print *, "Did not converge"
                        hasConverge = .true.
                    end if
                end if
            end do

            !Update velocity baseed on changed velocity and calc KE 
            !post-bondlength adjustments
            addedKE = 0.0
            addedKEScale = 0.0
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m,l,k) = (pos(m,l,k) - oldPos(m,l,k)) / timeStep
                    addedKE = addedKE + 0.5 * (0.5*(oldVel(m,l,k) + vel(m,l,k)))**2 * oMass 
                    addedKEScale = addedKEScale + 0.5 * vel(m,l,k)**2 * oMass
                end do
            end do

            kineticEnergy = kineticEnergy + addedKE
            kineticEnergyScale = kineticEnergyScale + addedKEScale
        end do

        !Calculate the total energy of the system
        totEnergy = Epot + kineticEnergy

        !Gradually scale down the temperature until it reaches a 
        !constant Temperature of 77K from steps 50,000 to 100,000
        if (p.LE.(numNVTSteps - numConstSteps).AND.p.GE.numConstSteps&
                                        &.AND.mod(p, temperatureStep) == 0) then
            print *, "temp: ", desiredTemperature
            desiredTemperature = desiredTemperature - 1
        end if

        !Scale the temperature of the system down to the desired value
        if (p.LE.numNVTSteps) then
            call scaleTemp(desiredTemperature, kineticEnergyScale, vel, numMolecules, &
                                     &numAtomsPerMolecule, numDimensions, Bolz)
        end if

        !Write energy and temperature information to files
        write(90, *) timestep * real(p), &
                    &(kineticEnergy * 2.0D0) / (5.0 * real(numMolecules) * Bolz)
        write(91, *) (timestep * real(p)), totEnergy
        write(92, *) (timestep * real(p)), Epot
        write(93, *) (timestep * real(p)), kineticEnergy

        !Call analysis functions to calculate TCF for velocity, 
        !rotation, and for MSD
        if(p.GE.numNVTSteps) then
            call TCF(p, Cvv, nCorr, vStore, vel)
            call Calc_MSD(p, MSD, nCorr_MSD, pStore, pos)
            call TCFRot(p, CvvRot, nCorr_Rot, orientStore, pos)
        end if

        !Output position and velocities for the trajectory file
        if (mod(p, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for O2Simulation. 100ns total time"
            write(94, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules 
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(94, 10) m, "OXY", "O1", 2*(m-1) + 1, &
                    &(pos(m,1,1) * 1E9) - 1E9*dim(1) * floor(avgPos(1)/dim(1)), &
                    &(pos(m,1,2) * 1E9) - 1E9*dim(1) * floor(avgPos(2)/dim(1)), &
                    &(pos(m,1,3) * 1E9) - 1E9*dim(1) * floor(avgPos(3)/dim(1)), & 
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(94, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9) - 1E9*dim(1) * floor(avgPos(1)/dim(1)), &
                    &(pos(m,2,2) * 1E9) - 1E9*dim(1) * floor(avgPos(2)/dim(1)), &
                    &(pos(m,2,3) * 1E9) - 1E9*dim(1) * floor(avgPos(3)/dim(1)), & 
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(94, 20) dim(1) * 1E9, dim(1) * 1E9, dim(1) * 1E9
        end if

        !Output the final frame of the simulation
        if (p == numTotSteps) then 
            write(95, *) "Last frame of the simulation of 500 O2 molecules."
            write(95, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(95, 10) m, "OXY", "O1", 2*(m-1) + 1, &
                    &(pos(m,1,1) * 1E9),(pos(m,1,2) * 1E9), (pos(m,1,3) * 1E9),&
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(95, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9),(pos(m,2,2) * 1E9),(pos(m,2,3) * 1E9),&
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(95, 20) dim(1) * 1E9, dim(1) * 1E9, dim(1) * 1E9
        end if
    end do

    !Write information about the velocity TCF
    CvvOne = Cvv(1)/(real(nCorr))
    do m = 1, CvvStep
        write(96, *) (real(m) * timestep) * 1e12, (Cvv(m)/(real(nCorr))) / CvvOne
        write(100, *) (real(m) * timestep) * 1e12, (Cvv(m)/(real(nCorr))) 
        diffusionCoeff = diffusionCoeff + (Cvv(m) / (3.0 * real(nCorr))) * timestep
    end do
    print *, "Desired Cvv1", 3.0 * Bolz * desiredTemperature / (oMass * 2.0),&
            &"Actual Cvv1", Cvv(1) / real(nCorr)
    print *, "DiffusionCoefficient: ", diffusionCoeff * 1e4

    !Write information about the Rotational TCF
    CvvOneRot = CvvRot(1)/real(nCorr_Rot)
    do m = 1, CvvStep
        write(99, *) (real(m) * timestep) * 1e12, (CvvRot(m)/real(nCorr_Rot)) / CvvOneRot
    end do

    !Write information about the MSD
    do m = 1, MSDStep
        write(98, *) (real(m) * timestep) * 1e12, (MSD(m)/real(nCorr_MSD)) * 1e20
    end do

    !Write information about g(r), the radial distribution function
    do m = 1, numFullBins
        rLower = real(m - 1) * delR    
        rUpper = rLower + delR
        shellVol = sphereConst * (rUpper ** 3 - rLower ** 3)
        write(97, *) (rlower + delR * 0.5)*1E9, real(bins(m)) / (real(numTotSteps - numNVTSteps) * &
                                                    &real(numAtoms) * shellVol)
    end do

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

    close (unit=90)
    close (unit=91)
    close (unit=92)
    close (unit=93)
    close (unit=94)
    close (unit=95)
    close (unit=96)
    close (unit=97)
    close (unit=98)
    close (unit=99)

    CALL cpu_time(end_time)
    print *, "Time usage:", (end_time - start_time)/60.0, " minutes", &
                            &mod((end_time - start_time),60.0), " seconds"

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

contains
    !Helper subroutine to initialize a 3D array of dimension 
    !(numMolecules, numAtomsPerMolecule, numDimensions) with all zeros. 
    !"array" is the 3D array passed and returned with all zeros
    subroutine initZero(array, numMolecules, numAtomsPerMolecule, numDimensions)
        implicit none
        
        integer, intent(in) :: numMolecules
        integer, intent(in) :: numAtomsPerMolecule
        integer, intent(in) :: numDimensions
        real(dp), dimension(numMolecules,numAtomsPerMolecule, numDimensions), &
                                                        &intent(inout) :: array

        integer :: m,l,k

        do m = 1, numMolecules 
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    array(m,l,k) = 0.0
                end do
            end do
        end do

        return

    end subroutine initZero

    !Zeros the net magnitude of the vectors passed in. Used to prevent the 
    !Wandering ice cube problem
    subroutine zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)
        implicit none
        
        integer, intent(in) :: numMolecules
        integer, intent(in) :: numAtomsPerMolecule
        integer, intent(in) :: numDimensions
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                          &intent(inout) :: vel

        integer :: m, l, k
    
        real(dp), dimension(numDimensions) :: netVel
        real(dp), dimension(numDimensions) :: velScale

        do k = 1, numDimensions
            netVel(k) = 0.0
        end do

        !Store the netVelcoity in each direction of the system
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    netVel(k) = netVel(k) + vel(m,l,k) 
                end do
            end do
        end do

        !Calc how much to scale each velocity by such that net velocity is zero
        do k = 1, numDimensions
            velScale(k) = netVel(k) / real(numMolecules * numAtomsPerMolecule)
        end do

        !Scale the velocity of each atom in the system
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m, l, k) = vel(m, l, k) - velScale(k)
                end do
            end do
        end do
    end subroutine zeroNetMomentum

    !Compute the dot product of two vectors v1 and v2
    real function dot(v1, v2)
        implicit none

        real(dp), dimension(numDimensions),intent(in) :: v1
        real(dp), dimension(numDimensions),intent(in) :: v2

        dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
        return
    end function

    !Scale the velocities in the system such that the system has a 
    !temperature equal to desiredTemperature
    subroutine scaleTemp(desiredTemperature, kineticEnergy, vel, numMolecules, &
                                            &numAtomsPerMolecule, numDimensions,Bolz)
        implicit none 
    
        real(dp), intent(in) :: kineticEnergy
        real(dp), intent(in) :: Bolz
        real(dp), intent(in) :: desiredTemperature
        integer, intent(in) :: numDimensions
        integer, intent(in) :: numAtomsPerMolecule
        integer, intent(in) :: numMolecules
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout)::vel

        real(dp) :: currTemp 
        real(dp) :: tempScale
        real(dp) :: currVel
        real(dp) :: newKE
        real(dp) :: nextVSq
        real(dp) :: addedKE
        real(dp) :: sumTotalVel
        integer :: degreesFreedom = 5

        integer :: m,l,k

        !Use equipartition thm to find the current temperature of the system
        currTemp = (kineticEnergy * 2.0D0) / &
                        &(real(degreesFreedom) * real(numMolecules) * Bolz)
        tempScale = sqrt(desiredTemperature / currTemp)

        !Scale each velocity such that the kintic energy of the system can 
        !be applied to the equipartition theorem to find that the system has 
        !a temperature of "desiredTemperature"
        do m=1, numMolecules
            do l=1, numAtomsPerMolecule
                do k=1, numDimensions
                    vel(m,l,k) = vel(m,l,k) * tempScale
                end do
            end do
        end do
    end subroutine scaleTemp

    !Test that the Shake Algorithm successfully fixes bond length
    subroutine testBondLengths(pos, numMolecules, numAtomsPerMolecule, numDimensions,&
                              &desiredBondLengthSq) 
        implicit none
    
        integer, intent(in) :: numMolecules 
        integer, intent(in) :: numAtomsPerMolecule
        integer, intent(in) :: numDimensions
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                        &intent(in) :: pos
        real(dp), intent(in) :: desiredBondLengthSq

        integer :: m,l,k

        real(dp) :: currBondLengthSq

        do m = 1, numMolecules
            currBondLengthSq = (pos(m, 1, 1) - pos(m, 2, 1))**2 +&
                             & (pos(m, 1, 2) - pos(m, 2, 2))**2 +&
                             & (pos(m, 1, 3) - pos(m, 2, 3))**2
            if (ABS(sqrt(desiredBondLengthSq) - sqrt(currBondLengthSq)).GE.1.0E-14) then
                print *, "desiredBondLength: ", sqrt(desiredBondLengthSq), &
                        &"currBondLength: ", sqrt(currBondLengthSq)
            end if
        end do
        
    end subroutine testBondLengths
    
    !Calculate the velocity TCF for the simulation at the current step add it 
    !to the running sum array Cvv
    subroutine TCF(step, Cvv, nCorr, vStore, vel)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr
        real(dp), dimension(CvvStep), intent(inout) :: Cvv
        real(dp), dimension(numMolecules, numDimensions, CvvStep), intent(inout) :: vStore
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: vel

        real(dp), dimension(CvvStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot product
                                            !is successfully added, then
                                            !this sum is added to the Cvv at the end
        integer :: index0
        integer :: currIndex
        integer :: m,l,k,q

        do q = 1, CvvStep
            sum(q) = 0.0
        end do

        currIndex = mod(step - 1, CvvStep) + 1 

        !Store record of the current velocities at the current index
        do m = 1, numMolecules
            do k = 1, numDimensions
                vStore(m, k, currIndex) = 0.5 * (vel(m,1,k) + vel(m,2,k))
            end do
        end do

        !Go thorugh add add all dot products for the current step after all 
        !vStores have been filled once
        if (step.GE.CvvStep) then
            nCorr = nCorr + 1
            index0 = mod(step, CvvStep) + 1
            do m = 1, numMolecules
                do k = 1, numDimensions
                    do q = 1, CvvStep
                        currIndex = mod(q - 1 + step, CvvStep) + 1
                        sum(q) = sum(q) + vStore(m,k,index0)*vStore(m,k,currIndex)
                    end do
                end do
            end do
        end if

        !Add sums to Cvv
        do q = 1, CvvStep
            Cvv(q) = Cvv(q) + sum(q) / real(numMolecules)
        end do
    end subroutine TCF

    !Calculate the velocity TCF for the simulation at the current step add it 
    !to the running sum array Cvv
    subroutine TCFRot(step, CvvRot, nCorr_Rot, orientStore, pos)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr_Rot
        real(dp), dimension(CvvStep), intent(inout) :: CvvRot
        real(dp), dimension(numMolecules, numDimensions, CvvStep), intent(inout) :: orientStore
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: pos

        real(dp), dimension(CvvStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot product
                                            !is successfully added, then
                                            !this sum is added to the Cvv at the end

        real(dp), dimension(numDimensions) :: orientV0
        real(dp) :: orientMag
        real(dp), dimension(numDimensions) :: orientVt
        integer :: index0
        integer :: currIndex
        integer :: m,l,k,q

        do q = 1, CvvStep
            sum(q) = 0.0
        end do

        currIndex = mod(step - 1, CvvStep) + 1 

        !Store record of the current molecular orientations at the current index
        do m = 1, numMolecules
            do k = 1, numDimensions
                orientStore(m, k, currIndex) = pos(m, 1, k) - pos(m, 2, k)
            end do
            orientMag = orientStore(m,1,currIndex) * orientStore(m,1,currIndex) +&
                       &orientStore(m,2,currIndex) * orientStore(m,2,currIndex) +&
                       &orientStore(m,3,currIndex) * orientStore(m,3,currIndex)
            orientMag = sqrt(orientMag)
            orientStore(m,1,currIndex) = orientStore(m,1,currIndex) / orientMag
            orientStore(m,2,currIndex) = orientStore(m,2,currIndex) / orientMag
            orientStore(m,3,currIndex) = orientStore(m,3,currIndex) / orientMag
        end do

        !Go thorugh and add to the sum array the second Legendre polynomial of the 
        !dot product of unit vectors of orientation 
        if (step.GE.CvvStep) then
            nCorr_Rot = nCorr_Rot + 1
            index0 = mod(step, CvvStep) + 1
            do m = 1, numMolecules
                orientV0(1) = orientStore(m,1,index0)
                orientV0(2) = orientStore(m,2,index0)
                orientV0(3) = orientStore(m,3,index0)
                do q = 1, CvvStep
                    currIndex = mod(q - 1 + step, CvvStep) + 1
                    orientVt(1) = orientStore(m,1,currIndex)
                    orientVt(2) = orientStore(m,2,currIndex)
                    orientVt(3) = orientStore(m,3,currIndex)
                    sum(q) = sum(q) + 0.5 * (3.0 * dot(orientV0, orientVt)**2 - 1.0)
                end do
            end do
        end if

        !Add sums to Cvv
        do q = 1, CvvStep
            CvvRot(q) = CvvRot(q) + sum(q) / real(numMolecules)
        end do
    end subroutine TCFRot

    !Calculate the mean square displacement at the current timestep
    subroutine Calc_MSD(step, MSD, nCorr_MSD, pStore, pos)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr_MSD
        real(dp), dimension(MSDStep), intent(inout) :: MSD
        real(dp), dimension(numMolecules, numDimensions, MSDStep), intent(inout) :: pStore
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: pos

        integer :: index0
        integer :: currIndex
        integer :: m, l, k, q

        real(dp), dimension(MSDStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot product
                                            !is successfully added, then
                                            !this sum is added to the MSD at the end

        do q = 1, MSDStep
            sum(q) = 0.0
        end do

        !Index of the pStore array you will store position information to this 
        !step
        currIndex = mod(step - 1, MSDStep) + 1 

        !Store the position of each molecule
        do m = 1, numMolecules
            do k = 1, numDimensions
                pStore(m, k, currIndex) = 0.5 * (pos(m, 1, k) + pos(m, 2, k))
            end do
        end do

        !After all of pStores indicies have been filled, calculate the displacement 
        !at the current timestep
        if (step.GE.MSDStep) then
            nCorr_MSD = nCorr_MSD + 1
            index0 = mod(step, MSDStep) + 1
            do m = 1, numMolecules
                do k = 1, numDimensions
                    do q = 1, MSDStep
                        currIndex = mod(q - 1 + step, MSDStep) + 1
                        sum(q) = sum(q) + (pStore(m, k, index0) - pStore(m, k, currIndex))**2
                    end do
                end do
            end do
        end if

        !Add the summed value to the total MSD
        do q = 1, MSDStep
            MSD(q) = MSD(q) + sum(q) / real(numMolecules)
        end do
    end subroutine Calc_MSD
end program nveSim
