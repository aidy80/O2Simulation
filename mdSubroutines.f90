!Subrouintes for md functions, force calculations, leapfrogs, initializations, etc

module mdSubroutines
    use Constants
    implicit none

contains
    !Helper subroutine to initialize a 2D array of dimension 
    !(arraySize, numDimensions) with all zeros. "array" is the 2D array passed 
    !and returned with all zeros
    subroutine initZero(array)
        implicit none
        
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
    end subroutine initZero

    subroutine openNVEFiles()
        integer :: numAtoms

        character :: nullChar

        !Open file containing position and velocity information 
        !about argon atoms in the cubic boundary from argon.gro
        open(unit=11, file='NVTO2_final.gro') 
        open(unit=90, file='NVEO2_temperature.dat')
        open(unit=91, file='NVEO2_totEnergy.dat')
        open(unit=92, file='NVEO2_potentialEnergy.dat')
        open(unit=93, file='NVEO2_kineticEnergy.dat')
        open(unit=94, file='NVEO2.gro')
        open(unit=95, file='NVEO2_final.gro')
        open(unit=96, file='TCFO2Norm.dat')
        open(unit=97, file='grO2.dat')
        open(unit=98, file='MSDO2.dat')
        open(unit=99, file='TCFRotO2.dat')
        open(unit=100, file='TCFO2.dat')
        open(unit=111, file="speedDistO2.dat") 
        open(unit=112, file="maxwellDistO2.dat")
        
        !Read in header information from the file
        read(11, *) nullChar
        read(11, 30) numAtoms
        if (numAtoms / numAtomsPerMolecule.NE.numMolecules) then
            print *, "Incorrect number of molecules read in"
        end if

        30 format(I5)
    end subroutine openNVEFiles

    subroutine openNVTFiles()
        character :: nullChar

        integer :: numAtoms

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
        if (numAtoms / numAtomsPerMolecule.NE.numMolecules) then
            print *, "Incorrect number of molecules read in"
        end if

        30 format(I5)
    end subroutine openNVTFiles
    
    !Read in position, velocity and dimension information from the initial file
    subroutine readIn(pos, vel, dim)
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: vel
        real(dp), dimension(numDimensions), intent(inout) :: dim

        !Garbage variables
        character :: nullChar
        character :: nullChar1
        integer :: nullInt
        integer :: nullInt1

        !Iterators
        integer :: m,l,k

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

        close (unit=11)

        !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    pos(m,l,k) = pos(m,l,k) * 1.0E-9
                    vel(m,l,k) = vel(m,l,k) * 1.0E3
                end do
            end do
        end do

        30 format(I5)
        20 format(F10.5, F10.5, F10.5)
        10 format(i5,2a5,i5,3f8.3,3f8.3)
    end subroutine readIn

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

    !Calculate the forces between each pair of atoms and record the net force
    !on each atom
    subroutine calcForces(potentialEnergy, pos, force, bins, numBins, dim, cutoffSq, ensemble)
        real(dp), intent(inout) :: potentialEnergy
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                               &intent(in) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                          &intent(inout) :: force
        integer, intent(in) :: numBins
        integer, dimension(numBins), intent(inout) :: bins
        real(dp), dimension(numDimensions), intent(in) :: dim
        real(dp), intent(in) :: cutoffSq
        character, intent(in) :: ensemble

        real(dp) :: sigmaDistTwo![]. (sigma/(distance between two atoms))**2. 
        real(dp) :: sigmaDistSix ![]. (sigma/(distance between two atoms))**6. 
        real(dp) :: sigmaDistTwelve![]. (sigma/(distance between two atoms))**12.  
        real(dp) :: potentialSum 
        real(dp) :: distanceSq
        real(dp) :: currPotential
        real(dp) :: Fmag
        real(dp), dimension(numDimensions) :: distance

        integer :: currIndex
        integer :: i,j,m,l,k

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
                            currPotential = fourEps * (sigmaDistTwelve - sigmaDistSix)
                            potentialSum = potentialSum + currPotential

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

                        !Add current position information to g(r) data
                        if (ensemble == "E") then
                            currIndex = Int(sqrt(distanceSq)/delR) + 1
                            bins(currIndex) = bins(currIndex) + 2
                        end if
                    end do
                end do
            end do
            potentialEnergy = potentialEnergy + potentialSum
        end do
    end subroutine calcForces

    !Update position and velocity using leapfrog algorithm and shake algorithm
    subroutine leapFrogAndShake(kineticEnergy, force, vel, pos, oldPos, oldVel)
        real(dp), intent(inout) :: kineticEnergy

        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: vel
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: oldPos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: oldVel
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: force

        real(dp) :: addedKE
        integer :: m, l, k

        !Use the leap-frog verlet algorithm to calculate new position and 
        !velocity vectors for all atoms based on the forces 
        !calculated between them.
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions 
                    oldVel(m, l, k) = vel(m, l, k)
                    oldPos(m, l, k) = pos(m, l, k)
                    vel(m, l, k) = vel(m, l, k) + (force(m, l, k) * invOMass) * timeStep
                    pos(m, l, k) = pos(m, l, k) + vel(m, l, k) * timeStep
                end do
            end do

            call Shake(m, oldPos, pos)

            !Update velocity baseed on changed velocity and calc KE 
            !post-bondlength adjustments
            addedKE = 0.0
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m,l,k) = (pos(m,l,k) - oldPos(m,l,k)) / timeStep
                    addedKE = addedKE + 0.5 * (0.5*(oldVel(m,l,k) + vel(m,l,k)))**2 * oMass 
                end do
            end do

            kineticEnergy = kineticEnergy + addedKE
        end do
    end subroutine leapFrogAndShake 

    !SHAKE Alg to fix bond angles
    subroutine Shake(m, oldPos, pos)
        integer, intent(in) :: m
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(in) :: oldPos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos

        logical :: hasConverge !Used in iterative while loop
        integer :: convergeCount

        real(dp), dimension(numDimensions) :: pastBondLength
        real(dp), dimension(numDimensions) :: currBondLength
        real(dp) :: currError
        real(dp) :: lambda
        real(dp) :: currBondLengthSq

        integer :: k
        
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
                pos(m, 1, k) = pos(m, 1, k) - pastBondLength(k) * invOMass * lambda
                pos(m, 2, k) = pos(m, 2, k) + pastBondLength(k) * invOMass * lambda
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

        if(m==1)then
            call testBondLengths(pos, numMolecules, numAtomsPerMolecule, &
                                        &numDimensions, desiredBondLengthSq)
        end if
    end subroutine Shake

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

    !Close all NVE files
    subroutine closeNVEFiles()
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
        close (unit=111)
        close (unit=112)
    end subroutine closeNVEFiles

    !Close all files
    subroutine closeNVTFiles()
        close (unit=91)
        close (unit=92)
        close (unit=93)
        close (unit=94)
        close (unit=95)
        close (unit=96)
    end subroutine closeNVTFiles

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

    !Compute the dot product of two vectors v1 and v2
    real function dot(v1, v2)
        implicit none

        real(dp), dimension(numDimensions),intent(in) :: v1
        real(dp), dimension(numDimensions),intent(in) :: v2

        dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
        return
    end function

end module mdSubroutines
