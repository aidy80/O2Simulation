!Reads in a gro files and outputs its information as positions in another gro file
!Given random velocity distribution, this should result in a circle appear with randomly
!Distributed points

program drawVelDist
    integer, parameter :: numMolecules = 500
    integer, parameter :: numAtomsPerMolecule = 2
    integer, parameter :: numDimensions = 3
    real, dimension(numMolecules, numAtomsPerMolecule, numDimensions) :: vel
    
    real :: zero = 0.0
    real :: dim

    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1
    real, dimension(numDimensions) :: nullReals

    integer :: m,l,k

    open(unit=11, file='oxygenInit.gro')
    open(unit=91, file='velCircle.gro')

    read (11, *) nullChar
    read (11, *) nullInt

    do m=1, numMolecules
        do l=1, numAtomsPerMolecule
            read(11, 10) nullInt, nullChar, nullChar1, nullInt1, &
                    & (nullReals(j), j=1,3), (vel(m,l,k), k=1,3)
        end do
    end do

    read(11, *) dim

    write(91, *) "Velocity circle"
    write(91, 30) numMolecules

    do m=1, numMolecules
        write(91, 10) m, "OXY", "O", m, (vel(m,1,k), k=1,3), (nullReals(j), j=1,3) 
    end do

    write(91, 20) dim, dim, dim

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)
end program drawVelDist
