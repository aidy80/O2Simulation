module testMod
    use Constants
    implicit none

contains
    !Test that the LJ potential is behaving as predicted, outputs a graph of the LJ
    !Potential as well as the measured and theoretical sigma and epsilon values
    subroutine testLJ(sigma, epsilon)
        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: epsilon

        integer, parameter :: numDistances = 120000
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
        real(dp) :: potential

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
            write(99, *) (real(i) * distStep) * 1E10, potential / 1000 * avagadrosNum
        end do

        close(unit=99)

        twoSix = 2.0**(1.0/6.0)
        theoMin = twoSix * sigma
        
        print *, "Theoretical minimum distance, sigma * 2^(1/6):", theoMin
        print *, "Theoretical minimum potential, -epsilon:", epsilon * (-1.0)
        print *, "Measured (sigma, epsilon): (", minDistance, minPotential, ")"
    end subroutine testLJ
end module testMod
