!Subroutines for analysis functions: g(r), TCF, MSD, etc

module analysisSubroutines
    use Constants
    implicit none

contains
    !Initialize several different analysis arrays
    subroutine initAnalysisArrays(Cvv, CvvRot, MSD, bins, numBins, velBoxes)
        real(dp), dimension(CvvStep), intent(inout) :: Cvv
        real(dp), dimension(CvvStep), intent(inout) :: CvvRot
        real(dp), dimension(MSDStep), intent(inout) :: MSD
        integer, intent(in) :: numBins
        integer, dimension(numBins), intent(inout) :: bins
        integer, dimension(numVelBox), intent(inout) :: velBoxes

        integer :: i

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

        do i = 1, numVelBox
            velBoxes(i) = 0
        end do
    end subroutine initAnalysisArrays

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

        !Store record of the current velocities at the current index
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

    !Add to MSD ensemble sum the displacement information of the current timestep
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
                        sum(q) = sum(q) + &
                                &(pStore(m, k, index0) - pStore(m, k, currIndex))**2
                    end do
                end do
            end do
        end if

        !Add the summed value to the total MSD
        do q = 1, MSDStep
            MSD(q) = MSD(q) + sum(q) / real(numMolecules)
        end do
    end subroutine Calc_MSD

    !Calculate the speed distribution of the oxygen simulation
    subroutine calcSpeedDist(vel, velBoxes)
        implicit none

        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions),&
                                                                &intent(in) :: vel
        integer, dimension(numVelBox), intent(inout) :: velBoxes

        integer :: index
        integer :: m,l,k
        
        real(dp) :: currVel

        !Calculate the current velocities and put them in appropriate boxes
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                currVel = 0.0
                do k = 1, numDimensions
                    currVel = currVel + vel(m,l,k)**2
                end do
                index = int(sqrt(currVel)*boxScale) + 1
                if (index > numVelBox) then
                    index = numVelBox
                end if
                velBoxes(index) = velBoxes(index) + 1
            end do
        end do
    end subroutine calcSpeedDist

    !Write energy and temperature information to files
    subroutine writeEnergy(p, kineticEnergy, potentialEnergy, totEnergy)
        integer, intent(in) :: p
        real(dp), intent(in) :: kineticEnergy
        real(dp), intent(in) :: potentialEnergy
        real(dp), intent(in) :: totEnergy

        write(90, *) timestep * real(p) * 1E12, &
                    &(kineticEnergy * 2.0D0) / (5.0 * real(numMolecules) * Bolz)
        write(91, *) (timestep * real(p)) * 1E12, totEnergy * avagadrosNum / &
                                                   &(real(numMolecules) * 1000.0) 
        write(92, *) (timestep * real(p)) * 1E12, potentialEnergy * avagadrosNum / &
                                                   &(real(numMolecules) * 1000.0) 
        write(93, *) (timestep * real(p)) * 1E12, kineticEnergy * avagadrosNum / &
                                                   &(real(numMolecules) * 1000.0) 
    end subroutine writeEnergy

    !Output position and velocities for the trajectory file
    subroutine writeNVETrajectory(p, dim, pos, vel)
        integer, intent(in) :: p
        real(dp), intent(in) :: dim
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                                &intent(in) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                                &intent(in) :: vel
        real(dp), dimension(numDimensions) :: avgPos

        integer :: m

        if (mod(p, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for NVE ensemble. 100ns total time"
            write(94, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules 
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(94, 10) m, "OXY", "O1", 2*(m-1) + 1, &
                    &(pos(m,1,1) * 1E9) - 1E9*dim * floor(avgPos(1)/dim), &
                    &(pos(m,1,2) * 1E9) - 1E9*dim * floor(avgPos(2)/dim), &
                    &(pos(m,1,3) * 1E9) - 1E9*dim * floor(avgPos(3)/dim), & 
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(94, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9) - 1E9*dim * floor(avgPos(1)/dim), &
                    &(pos(m,2,2) * 1E9) - 1E9*dim * floor(avgPos(2)/dim), &
                    &(pos(m,2,3) * 1E9) - 1E9*dim * floor(avgPos(3)/dim), & 
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(94, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if

        !Output the final frame of the NVE simulation
        if (p == numNVESteps) then 
            write(95, *) "Last frame of the NVE simulation of 500 O2 molecules."
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
            write(95, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if

        30 format(I5)
        20 format(F10.5, F10.5, F10.5)
        10 format(i5,2a5,i5,3f8.3,3f8.3)
    end subroutine writeNVETrajectory

    !Write trajectory files for NVT simulation
    subroutine writeNVTTrajectory(p, dim, pos, vel)
        integer, intent(in) :: p
        real(dp), intent(in) :: dim
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                                &intent(in) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                                &intent(in) :: vel
        real(dp), dimension(numDimensions) :: avgPos

        integer :: m

        !Output velocity and position information into a trajectory file
        if (mod(p, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for NVT ensemble. 100ns total time"
            write(94, 30) numMolecules * numAtomsPerMolecule
            do m = 1, numMolecules 
                avgPos(1) = 0.5 * (pos(m,1,1) + pos(m,2,1))
                avgPos(2) = 0.5 * (pos(m,1,2) + pos(m,2,2))
                avgPos(3) = 0.5 * (pos(m,1,3) + pos(m,2,3))
                write(94, 10) m, "OXY", "O1", 2*(m-1)+1, &
                    &(pos(m,1,1) * 1E9) - 1E9*dim * floor(avgPos(1)/dim), &
                    &(pos(m,1,2) * 1E9) - 1E9*dim * floor(avgPos(2)/dim), &
                    &(pos(m,1,3) * 1E9) - 1E9*dim * floor(avgPos(3)/dim), & 
                    &vel(m,1,1) * 1E-3, vel(m,1,2) * 1E-3, vel(m,1,3) * 1E-3

                write(94, 10) m, "OXY", "O2", 2*(m-1) + 2, &
                    &(pos(m,2,1) * 1E9) - 1E9*dim * floor(avgPos(1)/dim), &
                    &(pos(m,2,2) * 1E9) - 1E9*dim * floor(avgPos(2)/dim), &
                    &(pos(m,2,3) * 1E9) - 1E9*dim * floor(avgPos(3)/dim), & 
                    &vel(m,2,1) * 1E-3, vel(m,2,2) * 1E-3, vel(m,2,3) * 1E-3
            end do
            write(94, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if

        !Output the final frame of the NVT simulation
        if (p == numNVTSteps) then 
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
            write(95, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if

        30 format(I5)
        20 format(F10.5, F10.5, F10.5)
        10 format(i5,2a5,i5,3f8.3,3f8.3)
    end subroutine writeNVTTrajectory

    !Analyze TCFs including Velocity, Rotation, and Displacement
    subroutine TCF_Analysis(Cvv, nCorr, CvvRot, nCorr_Rot, MSD, nCorr_MSD)
        real(dp), dimension(CvvStep), intent(in) :: Cvv
        integer, intent(in) :: nCorr
        real(dp), dimension(CvvStep), intent(in) :: CvvRot
        integer, intent(in) :: nCorr_Rot
        real(dp), dimension(CvvStep), intent(in) :: MSD
        integer, intent(in) :: nCorr_MSD

        real(dp) :: CvvOne
        real(dp) :: CvvOneRot
        real(dp) :: diffusionCoeff = 0.0

        integer :: m

        !Write information about the TCF
        CvvOne = Cvv(1)/(real(nCorr))
        do m = 1, CvvStep
            write(96, *) (real(m) * timestep) * 1e12, (Cvv(m)/(real(nCorr))) / CvvOne
            write(100, *) (real(m) * timestep) * 1e12, (Cvv(m)/(real(nCorr))) 
            diffusionCoeff = diffusionCoeff + (Cvv(m) / (real(nCorr))) * timestep / 3.0
        end do
        print *, "Desired Cvv1", 3.0 * Bolz * finalTemperature / (oMass * 2.0),&
                &"Actual Cvv1", Cvv(1) / real(nCorr)
        print *, "DiffusionCoefficient: ", diffusionCoeff * 1e4

        !Write information about the TCF
        CvvOneRot = CvvRot(1)/real(nCorr_Rot)
        do m = 1, CvvStep
            write(99, *) (real(m) * timestep) * 1e12, &
                                        &(CvvRot(m)/real(nCorr_Rot)) / CvvOneRot
        end do

        !Write information about the MSD
        do m = 1, MSDStep
            write(98, *) (real(m) * timestep) * 1e12, (MSD(m)/real(nCorr_MSD)) * 1e20
        end do
    end subroutine TCF_Analysis

    !Complete analysis for g(r) and velocity distribution
    subroutine Bins_Analysis(velBoxes, bins, numBins, numFullBins, dim)
        integer, dimension(numVelBox), intent(in) :: velBoxes

        integer, intent(in) :: numBins
        integer, dimension(numBins), intent(in) :: bins
        real(dp), intent(in) :: dim

        integer, intent(in) :: numFullBins

        real(dp) :: rLower
        real(dp) :: rUpper
        real(dp) :: rho
        real(dp) :: sphereConst
        real(dp) :: shellVol

        integer :: i, m

        rho = numMolecules * numAtomsPerMolecule / (dim**3)
        sphereConst = 4.0D0 * pi * rho / 3.0D0

        
        !Write my speedDist
        do i = 1, numVelBox
            write(111, *) (real(i - 1) / boxScale), real(velBoxes(i)) / &
                            &real(numMolecules * numAtomsPerMolecule * numNVESteps)
        end do

        !Write Theoretical Maxwell Distribution
        do i = 1, numMaxwell
            write(112, *) real(i - 1), 4.0 * pi * &
                                &sqrt(oMass / (2.0*pi*Bolz*finalTemperature))**3 * &
                                &real(i - 1)**2 * e**(oMass*(real(i - 1))**2/ &
                                                        &(-2.0*Bolz*finalTemperature))
        end do        

        !Write information about g(r), the radial distribution function
        do m = 1, numFullBins
            rLower = real(m - 1) * delR    
            rUpper = rLower + delR
            shellVol = sphereConst * (rUpper ** 3 - rLower ** 3)
            write(97, *) (rLower + delR * 0.5)*1E9, real(bins(m)) / (real(numNVESteps) * &
                                &real(numMolecules * numAtomsPerMolecule) * shellVol)
        end do
    end subroutine Bins_Analysis

    !Compute the dot product of two vectors v1 and v2
    real function dot(v1, v2)
        implicit none

        real(dp), dimension(numDimensions),intent(in) :: v1
        real(dp), dimension(numDimensions),intent(in) :: v2

        dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
        return
    end function
end module analysisSubroutines
