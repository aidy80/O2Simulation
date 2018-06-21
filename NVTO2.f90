!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 100,000 4fs timesteps in a NVT
!ensemble. Temperature will begin at 320K and be gradually scaled to 77K

program nvtSim
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)!Give reals double precision

    integer, parameter :: numAtomsPerMolecule = 2 !Number of atoms in each Oxygen 
    integer :: numMolecules !Number of molecules in this simulation
    integer :: numAtoms !Number of atoms in the simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]
    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel !oldPos: [m]
                                                                !oldVel: [m/s]

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 

    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space

    real(dp), dimension(numDimensions) :: pastBondLength !Holds bond length of an 
                                                         !O2 in its previous iteration
    real(dp), dimension(numDimensions) :: currBondLength !Holds current bond length of an O2
    real(dp) :: currBondLengthSq
    real(dp), dimension(numDimensions) :: avgPos !Temp variable to measure 
                                                 !pos of an O2 molecule
    real(dp) :: lambda
    
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
    real(dp) :: kineticEnergyScale ![J] The total kinetic energy of the system
    real(dp) :: addedKE ![J] Kinetic Energy sum to prevent numerical 
                        !percision issues
    real(dp) :: addedKEScale ![J] Kinetic Energy sum to prevent numerical 
                             !percision issues
    real(dp) :: oldKE

    !Used to time simulation
    real :: start_time
    real :: end_time

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance ![m]
    real(dp) :: distanceSq ![m^2]

    !Iterators
    integer :: i, j, k, m, l, n, p
    character :: oxygenChar*2

    !Constants
    real(dp), parameter :: oMass = 2.66E-26 ![kg] Mass of an oxygen atom
    real(dp),parameter :: invOMass = 1.0 / oMass
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: epsilon = 48.0 * Bolz![J] Minimum in the 
                                                 !    Lennard Jones equation
    real(dp), parameter :: sigma = 3.006E-10 ![m]. Constant in the Lennard 
                                            !     Jones equation
    real(dp), parameter :: sigmaSq = sigma**2 ![m^2] Sigma squared
    real(dp) :: sigmaDistTwo![]. (sigma/(distance between two atoms))**2. 
    real(dp) :: sigmaDistSix ![]. (sigma/(distance between two atoms))**6. 
    real(dp) :: sigmaDistTwelve![]. (sigma/(distance between two atoms))**12.  
    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure
    real(dp), parameter :: initTemperature = 320.0 ![K] initial const temperature 
                                                    !of the system
    real(dp), parameter :: finalTemperature = 77.0 ![K] final temperature of 
                                                    !the system
    real(dp) :: desiredTemperature = initTemperature

    real :: cutoff ![m] Cutoff distance for short-range interaction
    real, parameter :: timeStep = 4.0E-15 ![s] Time between calculations
    real(dp) :: fourEps = 4.0 * epsilon ![J] Epsilon*4. Used for optimization
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared
    real(dp) :: twentyFourEps = 24.0 * epsilon ![J] Epsilon*24. 
                                               !        Used for optimization
    integer, parameter :: numSteps = 150000 !Number of timesteps in the program
    integer, parameter :: numConstSteps = 50000 !Number of timesteps that the 
                                                !system is equilibrated at a 
                                                !constant temperature

    integer :: convergeCount = 0 !Counts number of iterations
    logical :: hasConverge = .false. !Used in iterative while loop

    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                   !between momentum-zeroing
    integer, parameter :: numTrajSteps = 100 !Number of timesteps between 
                                            !trajectory outputs to .ttr file
    integer, parameter :: temperatureStep = 205 !Number of timesteps between 
                                                !temperature decrements
    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1

    CALL cpu_time(start_time)

    !Open file containing position and velocity information 
    !about argon atoms in the cubic boundary from argon.gro
    open(unit=11, file='oxygenInit.gro') 
    open(unit=91, file='NVTO2_totEnergy.dat')
    open(unit=92, file='NVTO2_potentialEnergy.dat')
    open(unit=93, file='NVTO2_kineticEnergy.dat')
    open(unit=94, file='NVTO2.gro')
    open(unit=95, file='NVTO2_final.gro')
    open(unit=96, file='NVTO2_temperature.dat')
    
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

    !Set the cutoff for LJ interation to 1/2 box size
    cutoff = dim(1)/2.0
    cutoffSq = cutoff**2

    !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
    do m = 1, numMolecules
        do l = 1, numAtomsPerMolecule
            do k = 1, numDimensions
                pos(m,l,k) = pos(m,l,k) * 1.0E-9
                vel(m,l,k) = vel(m,l,k) * 1.0E3
            end do
        end do
    end do

    close (unit=11)

    !call testLJ(sigma, epsilon)

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel, numMolecules, numAtomsPerMolecule, numDimensions)

    do p = 1, numSteps
        if(mod(p,100) == 0) then
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
                            Epot = Epot + potential

                            !Calculate the resulting force on the current two atoms 
                            !based on the lennard jones potential between them. Calculated 
                            !using the negative gradient
                            Fmag = twentyFourEps * (2.0 * sigmaDistTwelve - sigmaDistSix)

                            do k = 1, numDimensions
                                force(i, m, k) = force(i, m, k) + &
                                                    &Fmag * (distance(k) / distanceSq)
                                force(j, l, k) = force(j, l, k) - &
                                                    &Fmag * (distance(k) / distanceSq)
                            end do
                        end if

                        !if ((i == 1).AND.(m == 1).AND.(j == 2).AND.(l == 1).AND.(p == 1)) then
                        !    print *, "Dist", distance(1), distance(2), distance(3)
                        !    print *, "FORCE", force(1,1,1), force(1,1,2), force(1,1,3)
                        !end if

                    end do
                end do
            end do
        end do


        !Use the leapfrog verlet algorithm to calculate new position and 
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
                currError = ABS(currBondLengthSq - desiredBondLengthSq) / desiredBondLengthSq

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
                    end if
                end if
            end do

            !if (m ==1) then
                !call testBondLengths(pos, numMolecules, numAtomsPerMolecule, &
            !                                &numDimensions, desiredBondLengthSq)
            !end if

            !Update velocity baseed on changed velocity and calc KE 
            !post-bondlength adjustments
            !oldKE = kineticEnergy
            addedKE = 0.0
            addedKEScale = 0.0
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m,l,k) = (pos(m,l,k) - oldPos(m,l,k)) / timeStep
                    addedKE = addedKE + 0.5 * (0.5 * (oldVel(m,l,k) + vel(m,l,k)))**2 * oMass 
                    addedKEScale = addedKEScale + 0.5 * vel(m,l,k)**2 * oMass
                end do
            end do

            kineticEnergyScale = kineticEnergyScale + addedKEScale
            kineticEnergy = kineticEnergy + addedKE

            !if ((kineticEnergy - oldKE) == 0) then
            !    print *, "uhoh", vel(m,l,k) *0.5*oMass
            !end if
        end do
        
        !Find the kinetic energy of the system and 
        !the total energy of the system
        totEnergy = Epot + kineticEnergy

        write(91, *) (timestep * real(p)), totEnergy
        write(92, *) (timestep * real(p)), Epot
        write(93, *) (timestep * real(p)), kineticEnergy
        write(96, *) (timestep * real(p)), (kineticEnergy * 2.0) / &
                                              &(5.0 * real(numMolecules) * Bolz)

        !Gradually scale down the temperature until it reaches a 
        !constant Temperature of 77K
        if (p.LE.(numSteps - numConstSteps).AND.p.GE.numConstSteps&
                                        &.AND.mod(p, temperatureStep) == 0) then
            print *, "temp: ", desiredTemperature
            desiredTemperature = desiredTemperature - 1
        end if

        !Scale the temperature of the system down to the desired value
        call scaleTemp(desiredTemperature, kineticEnergyScale, vel, numMolecules, &
                                     &numAtomsPerMolecule, numDimensions, Bolz)

        !Output position and velocities for the trajectory file
        if (mod(p-1, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for NVT ensemble. 100ns total time"
            write(94, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules 
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(94, 10) m, "OXY", "O1", 2*(m-1)+1, &
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

        !Output the final frame of the NVT simulation
        if (p == numSteps) then 
            write(95, *) "Final Frame of the NVT ensemble. 500 O2 molecules at 77K"
            write(95, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules
                write(95, 10) m, "OXY", "O1", 2*(m-1) + 1, &
                    &(pos(m,1,1) * 1E9), (pos(m,1,2) * 1E9),(pos(m,1,3) * 1E9),&
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(95, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9),(pos(m,2,2) * 1E9),(pos(m,2,3) * 1E9),&
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(95, 20) dim(1) * 1E9, dim(1) * 1E9, dim(1) * 1E9
        end if
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(oldPos)
    deallocate(oldVel)
    deallocate(force)

    !Close all files
    close (unit=91)
    close (unit=92)
    close (unit=93)
    close (unit=94)
    close (unit=95)
    close (unit=96)

    CALL cpu_time(end_time)
    print *, "Time usage:", (end_time - start_time)/60.0, " minutes", &
                            &mod((end_time - start_time),60.0), " seconds"

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

contains
    !Helper subroutine to initialize a 2D array of dimension 
    !(arraySize, numDimensions) with all zeros. "array" is the 2D array passed 
    !and returned with all zeros
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

        !Find the current velocity of the system
        do k = 1, numDimensions
            netVel(k) = 0.0
        end do

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

        !TESTING
        !print *, "currTemp: ", currTemp
        !newKE = 0.0
        !sumTotalVel=0.0

        !Scale each velocity such that the kintic energy of the system can 
        !be applied to the equipartition theorem to find that the system has 
        !a temperature of "desiredTemperature"
        do m=1, numMolecules
            do l=1, numAtomsPerMolecule
                !currVel = 0.0
                do k=1, numDimensions
                    vel(m,l,k) = vel(m,l,k) * tempScale
                    !currVel = currVel + vel(m,l,k)**2
                    !newKE = newKE + 0.5 * vel(m,l,k)**2 * oMass
                end do
                !sumTotalVel = sumTotalVel + currVel
            end do
        end do

        !print *, "Desired: ", 0.644374E3, "ACtual:",&
        !                &sqrt(sumTotalVel / real(numAtomsPerMolecule * numMolecules))

        !currTemp = (newKE * 2.0D0) / (degreesFreedom * numMolecules * Bolz)
        !print *, "New temp: ", currTemp
    end subroutine scaleTemp

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
            if (ABS(sqrt(desiredBondLengthSq) - sqrt(currBondLengthSq)).GE.1.0E-14) then
                print *, "desiredBondLength: ",sqrt(desiredBondLengthSq), &
                        &"currBondLength: ",sqrt(currBondLengthSq)
            else if (p==500) then
                print *, "desiredBondLength: ",sqrt(desiredBondLengthSq), &
                        &"currBondLength: ",sqrt(currBondLengthSq)
            end if
        end do
        
    end subroutine testBondLengths

    !Test that the LJ potential is behaving as predicted, outputs a graph of the LJ
    !Potential as well as the measured and theoretical sigma and epsilon values
    subroutine testLJ(sigma, epsilon)
        implicit none

        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: epsilon

        integer, parameter :: numDistances = 100000
        real, parameter :: distStep = 1.0E-14
        real, parameter :: avagadrosNum = 6.0221409E23
        real(dp) :: twoSix
        real(dp) :: theoMin
        real(dp) :: distanceSq
        real(dp) :: sigmaDistTwo
        real(dp) :: sigmaDistSix
        real(dp) :: sigmaDistTwelve

        real(dp) :: minPotential = 1000000.0
        real(dp) :: minDistance

        integer :: i

        open(unit=99, file='lennardJones.dat') 

        !For each distance, record the LJ potential and output it to the file
        do i = 30000, numDistances
            distanceSq = (real(i) * distStep)**2
            sigmaDistTwo = (sigma**2 / distanceSq)
            sigmaDistSix = sigmaDistTwo * sigmaDistTwo * sigmaDistTwo
            sigmaDistTwelve = sigmaDistSix**2

            potential = 4.0 * epsilon * (sigmaDistTwelve - sigmaDistSix)
            if (potential < minPotential) then
                minPotential = potential
                minDistance = real(i) * distStep
            end if

            !Write potential is KJ/mol
            write(99, *) (real(i) * distStep) * 1E9, potential / 1000 * avagadrosNum
        end do

        close(unit=99)

        twoSix = 2.0**(1.0/6.0)
        theoMin = twoSix * sigma
        
        print *, "Theoretical minimum distance, sigma * 2^(1/6):", theoMin
        print *, "Theoretical minimum potential, -epsilon:", epsilon * (-1.0)
        print *, "Measured (sigma, epsilon): (", minDistance, minPotential, ")"
    end subroutine
end program nvtSim
