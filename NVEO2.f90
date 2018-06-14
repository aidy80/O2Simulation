!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 100,000 4fs timesteps in an 
!NVE ensemble

program nveSim
    implicit none

    integer, parameter :: dp = kind(1.0)!Give reals double precision

    integer, parameter :: numAtomsPerMolecule = 2 !Number of atoms in each Oxygen 
    integer :: numMolecules !Number of molecules in this simulation
    integer :: numAtoms !Number of atoms in the simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]
    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel !oldPos: [m]

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 
    real(dp), dimension(:), allocatable :: distBins ![]

    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space
    integer, parameter :: CvvStep = 500 !Number of items in Cvv
    integer, parameter :: MSDStep = 500 !Number of items in MSD

    real(dp), dimension(numDimensions) :: pastBondLength !Holds bond length of an 
                                                         !O2 in its previous iteration
    real(dp), dimension(numDimensions) :: currBondLength !Holds current bond length of an O2
    real(dp) :: currBondLengthSq
    real(dp), dimension(numDimensions) :: avgPos !Used to hold median bond position
    real(dp) :: lambda!Coefficient which must be solved for through iteration
    
    real(dp), parameter :: desiredBondLength = 1.208E-10 ![m]Bond length betweem 
                                                         !oxygen atoms
    real(dp), parameter :: desiredBondLengthSq = desiredBondLength**2![m^2]
    real(dp), parameter :: allowedError = 1.0E-13 ![m]Allowed error between actual 
                                                  !bond length and desired
    real(dp) :: currError

    real(dp) :: Epot ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: vSQ ![(m/s)^2] Square velocity of a given atom
    real(dp) :: vOld ![m/s] Temp variable for velocity
    real(dp) :: Fmag ![N] Used to hold the magnitude of force btwn two atoms
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system

    !Used to time simulation
    real :: start_time
    real :: end_time

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance ![m]
    real(dp) :: distanceSq ![m^2]

    real(dp), dimension(CvvStep) :: Cvv
    real(dp), dimension(CvvStep) :: CvvRot
    real(dp) :: CvvOne = 0
    real(dp) :: CvvOneRot = 0
    real(dp), dimension(MSDStep) :: MSD
    real(dp), dimension(:, :, :), allocatable :: vStore
    real(dp), dimension(:, :, :), allocatable :: pStore
    real(dp), dimension(:, :, :), allocatable :: orientStore
    integer :: nCorr = 0!Used to count the number of times Cvv is added to 
    integer :: nCorr_MSD = 0!USed to count the number of times MSD is added to 
    integer :: nCorr_Rot = 0!USed to count the number of times CvvRot is added to 
    real(dp) :: diffusionCoeff = 0!Diffusion coefficient to be calcualted

    !Iterators
    integer :: i, j, k, m, l, n, p
    character :: oxygenChar*2

    !Constants
    real(dp), parameter :: oMass = 2.66E-26 ![kg] Mass of an oxygen atom
    real(dp),parameter :: invOMass = 1 / oMass
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: epsilon = 48 * Bolz![J] Minimum in the 
                                                 !    Lennard Jones equation
    real(dp), parameter :: sigma = 3.006E-10 ![m]. Constant in the Lennard 
                                            !     Jones equation
    real(dp), parameter :: sigmaSq = sigma**2 ![m^2] Sigma squared
    real(dp) :: sigmaDistTwo![]. (sigma/(distance between two atoms))**2. 
    real(dp) :: sigmaDistSix ![]. (sigma/(distance between two atoms))**6. 
    real(dp) :: sigmaDistTwelve![]. (sigma/(distance between two atoms))**12.  
    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure
    real(dp), parameter :: temperature = 77.0 ![K] final temperature of 

    real :: cutoff ![m] Cutoff distance for short-range interaction
    real, parameter :: timeStep = 4.0E-15 ![s] Time between calculations
    real(dp) :: fourEps = 4.0 * epsilon ![J] Epsilon*4. Used for optimization
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared
    real(dp) :: twentyFourEps = 24.0 * epsilon ![J] Epsilon*24. 
                                               !        Used for optimization
    integer, parameter :: numSteps = 12000 !Number of timesteps in the program

    integer :: convergeCount = 0 !Counts number of iterations
    logical :: hasConverge = .false. !Used in iterative while loop

    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                   !between momentum-zeroing
    integer, parameter :: numTrajSteps = 10 !Number of timesteps between 
                                            !trajectory outputs to .ttr file

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
    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1

    CALL cpu_time(start_time)

    !Open file containing position and velocity information 
    !about argon atoms in the cubic boundary from argon.gro
    open(unit=11, file='NVTO2_final.gro') 
    open(unit=91, file='NVEO2_totEnergy.dat')
    open(unit=92, file='NVEO2_potentialEnergy.dat')
    open(unit=93, file='NVEO2_kineticEnergy.dat')
    open(unit=94, file='NVEO2.gro')
    open(unit=95, file='NVEO2_final.gro')
    open(unit=96, file='TCFO2.dat')
    open(unit=97, file='grO2.dat')
    open(unit=98, file='MSDO2.dat')
    open(unit=99, file='TCFRotO2.dat')
    
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
    rho = numAtoms / (dim(1) ** 3)
    sphereConst = 4.0D0 * pi * rho / 3.0D0
    numBins = ANint((dim(1)*1.733)/(2.0 * delR)) + 1
    numFullBins = ANint(dim(1)/(2.0*delR))

    !Set the cutoff for LJ interation to 1/2 box size
    cutoff = dim(1)/2
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

    do i = 1, MSDStep
        MSD(i) = 0.0D0
    end do

    do i = 1, numBins
        bins(i) = 0
    end do

    do p = 1, numSteps
        !Init force vectors to zero
        call initZero(force, numAtoms, numDimensions)

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        if (mod(p, zeroMomentTimeStep) == 0) then
            call zeroNetMomentum(vel, numAtoms, numAtomsPerMolecule, numDimensions)
        end if

        Epot = 0.0
        kineticEnergy = 0.0

        !Compute Force between each pair of atoms if their distance is below
        !the cutoff radius. Each pair is evaluated only once
        do i = 1, numMolecules - 1
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

                        sigmaDistTwo = sigmaSq / distanceSq
                        sigmaDistSix = sigmaDistTwo**3
                        sigmaDistTwelve = sigmaDistSix**2

                        !Calc potential from lennard jones equation between the 
                        !current pair of atoms and add it to the current 
                        !potential energy sum
                        potential = fourEps * (sigmaDistTwelve - sigmaDistSix)
                        Epot = Epot + potential

                        !Calculate the resulting force on the current two atoms 
                        !based on the lennard jones potential between them. Calculated 
                        !using the negative gradient
                        Fmag = twentyFourEps * (2 * sigmaDistTwelve - sigmaDistSix)

                        !If the distance between the two atoms is below the cutoff, 
                        !calculate the force exerted on each of them based on the 
                        !lennard jones potential
                        if (distanceSq < cutoffSq) then
                            do k = 1, numDimensions
                                force(i, m, k) = force(i, m, k) + &
                                                    &Fmag * (distance(k) / distanceSq)
                                force(j, l, k) = force(j, l, k) - &
                                                    &Fmag * (distance(k) / distanceSq)
                            end do
                        end if

                        currIndex = Int(sqrt(distanceSq)/delR) + 1
                        bins(currIndex) = bins(currIndex) + 2
                    end do
                end do
            end do
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
                    !vel(m, 1, k) = vel(m, 1, k) - pastBondLength(k)*invOMass/timestep * lambda
                    !vel(m, 2, k) = vel(m, 2, k) + pastBondLength(k)*invOMass/timestep * lambda
                    currBondLengthSq = currBondLengthSq + (pos(m, 1, k) - pos(m, 2, k))**2 
                end do

                !Find the error between the current bond 
                !length and desired bond length
                currError = ABS(currBondLengthSq - desiredBondLengthSq) 

                !Repeat the iteration until the bond is within 1E-18 of the desired 
                !length or the iteration occurs 500 times
                if(currError.LE.allowedError) then
                    hasConverge = .true.
                else 
                    hasConverge = .false.
                    convergeCount = convergeCount + 1

                    if (convergeCount > 500) then
                        print *, "Did not converge"
                        hasConverge = .true.
                    else
                        !print *, "Did converge"
                    end if
                end if
            end do

            !call testBondLengths(pos, numMolecules, numAtomsPerMolecule, &
            !                                &numDimensions, desiredBondLengthSq)

            !Update velocity baseed on changed velocity and calc KE 
            !post-bondlength adjustments
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m,l,k) = (pos(m,l,k) - oldPos(m,l,k)) / timeStep
                    kineticEnergy = kineticEnergy + 0.5 * &
                              &(0.5 * (oldVel(m,l,k) + vel(m,l,k)))**2 * oMass 
                end do
            end do
        end do

        !Find the kinetic energy of the system and 
        !the total energy of the system
        totEnergy = Epot + kineticEnergy

        !Write energy information to files
        write(91, *) (timestep * real(p)), totEnergy
        write(92, *) (timestep * real(p)), Epot
        write(93, *) (timestep * real(p)), kineticEnergy

        call TCF(p, Cvv, nCorr, vStore, vel)
        call Calc_MSD(p, MSD, nCorr_MSD, pStore, pos)
        call TCFRot(p, Cvv, nCorr_Rot, orientStore, pos)

        !Output position and velocities for the trajectory file
        if (mod(p, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for NVE ensemble. 100ns total time"
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

        !Output the final frame of the NVE simulation
        if (p == numSteps) then 
            write(95, *) "Last frame of the NVE simulation of 500 O2 molecules."
            write(95, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(95, 10) m, "OXY", "O1", 2*(m-1) + 1, &
                    &(pos(m,1,1) * 1E9),& 
                    &(pos(m,1,2) * 1E9),&
                    &(pos(m,1,3) * 1E9),&
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(95, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9),&
                    &(pos(m,2,2) * 1E9),&
                    &(pos(m,2,3) * 1E9),&
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(95, 20) dim(1) * 1E9, dim(1) * 1E9, dim(1) * 1E9
        end if
    end do

    !Write information about the TCF
    CvvOne = Cvv(1)/real(nCorr)
    do m = 1, CvvStep
        if (m < 200) then
            write(96, *) (real(m) * timestep) * 1e12, (Cvv(m)/real(nCorr)) / CvvOne
        end if
        diffusionCoeff = diffusionCoeff + (Cvv(m) / real(nCorr)) * timestep
    end do
    print *, "DiffusionCoefficient: ", diffusionCoeff * 1e4

    !Write information about the TCF
    CvvOneRot = CvvRot(1)/real(nCorr_Rot)
    do m = 1, CvvStep
        if (m < 200) then
            write(99, *) (real(m) * timestep) * 1e12, (CvvRot(m)/real(nCorr_Rot)) / CvvOneRot
        end if
    end do

    !Write information about the MSD
    do m = 1, MSDStep
        write(98, *) (real(m) * timestep) * 1e12, &
                                &(MSD(m)/real(nCorr_MSD)) * 1e20
    end do

    !Write information about g(r), the radial distribution function
    do m = 1, numFullBins
        rLower = real(m - 1) * delR    
        rUpper = rLower + delR
        shellVol = sphereConst * (rUpper ** 3 - rLower ** 3)
        write(97, *) (rlower + delR * 0.5) / sigma, real(bins(m)) / (real(numSteps) * &
                                                    &real(numAtoms) * shellVol)
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(force)
    deallocate(vStore)
    deallocate(pStore)
    deallocate(orientStore)
    deallocate(bins)

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
    print *, "Time usage:", (end_time - start_time)/60, " minutes"
    print *, "Time usage:", (end_time - start_time) - &
                                       &(end_time - start_time)/60, " seconds"

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

contains
    !Helper subroutine to initialize a 2D array of dimension 
    !(arraySize, numDimensions) with all zeros. "array" is the 2D array passed 
    !and returned with all zeros
    subroutine initZero(array, arraySize, numDimensions)
        implicit none
        
        integer, intent(in) :: arraySize
        integer, intent(in) :: numDimensions
        real(dp), dimension(arraySize, numDimensions), intent(inout) :: array

        integer :: i
        integer :: k

        do i = 1, arraySize 
            do k = 1, numDimensions
                array(i,k) = 0.0
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

        !Find the current velocity of the system
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
            print *, "desiredBondLength: ",sqrt(desiredBondLengthSq), &
                                    &"currBondLength: ",sqrt(currBondLengthSq)
        end do
        
    end subroutine testBondLengths
    
    !Calculate the velocity TCF for the simulation at the current step add it 
    !to the running sum array Cvv
    subroutine TCF(step, Cvv, nCorr, vStore, vel)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr
        real, dimension(CvvStep), intent(inout) :: Cvv
        real, dimension(numMolecules, numDimensions, CvvStep), intent(inout) :: vStore
        real, dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: vel

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
            Cvv(q) = Cvv(q) + sum(q) / real(numMolecules*numDimensions)
        end do
    end subroutine TCF

    !Calculate the velocity TCF for the simulation at the current step add it 
    !to the running sum array Cvv
    subroutine TCFRot(step, CvvRot, nCorr, orientStore, pos)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr
        real, dimension(CvvStep), intent(inout) :: CvvRot
        real, dimension(numMolecules, numDimensions, CvvStep), intent(inout) :: orientStore
        real, dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: pos

        real(dp), dimension(CvvStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot product
                                            !is successfully added, then
                                            !this sum is added to the Cvv at the end

        real(dp), dimension(numDimensions) :: orientV0
        real(dp), dimension(numDimensions) :: orientVt
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
                orientStore(m, k, currIndex) = pos(m, 1, k) - pos(m, 2, k)
            end do
        end do

        !Go thorugh add add all dot products for the current step after all 
        !vStores have been filled once
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
                    sum(q) = sum(q) + 0.5 * (3 * dot(orientV0, orientVt)**2 - 1)
                end do
            end do
        end if

        !Add sums to Cvv
        do q = 1, CvvStep
            CvvRot(q) = CvvRot(q) + sum(q) / real(numMolecules)
        end do
    end subroutine TCFRot

    subroutine Calc_MSD(step, MSD, nCorr_MSD, pStore, pos)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr_MSD
        real, dimension(MSDStep), intent(inout) :: MSD
        real, dimension(numMolecules, numDimensions, MSDStep), intent(inout) :: pStore
        real, dimension(numMolecules, numAtomsPerMolecule, numDimensions), intent(in) :: pos

        integer :: index0
        integer :: currIndex
        integer :: m, l, k, q

        real(dp), dimension(MSDStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot product
                                            !is successfully added, then
                                            !this sum is added to the Cvv at the end

        do q = 1, MSDStep
            sum(q) = 0.0
        end do

        currIndex = mod(step - 1, MSDStep) + 1 

        !Store the position of each molecule
        do m = 1, numAtoms
            do k = 1, numDimensions
                pStore(m,k,currIndex) = 0.5 * (pos(m,1,k) + pos(m,2,k))
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
