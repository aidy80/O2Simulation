module Constants 
    integer, parameter :: dp = selected_real_kind(15, 307)!Give reals double precision

    integer, parameter :: numAtomsPerMolecule = 2 !Number of atoms in each Oxygen 
    integer, parameter :: numMolecules = 500 !Number of molecules in this simulation
    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space
    integer, parameter :: CvvStep = 500 !Number of items in Cvv
    integer, parameter :: MSDStep = 500 !Number of items in MSD

    real(dp), parameter :: desiredBondLength = 1.208E-10 ![m]Bond length betweem 
                                                         !oxygen atoms
    real(dp), parameter :: desiredBondLengthSq = desiredBondLength**2![m^2]
    real(dp), parameter :: allowedError = 1.0E-13 ![m]Allowed error between actual 
                                                  !bond length and desired

    !Constants
    real(dp), parameter :: oMass = 2.66E-26 ![kg] Mass of an oxygen atom
    real(dp), parameter :: e = 2.71828 !Mathematical constant
    real(dp), parameter :: invOMass = 1.0 / oMass
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: avagadrosNum = 6.0221409E23![molecules/mole]
    real(dp), parameter :: pi = 3.14159265359
    real(dp), parameter :: epsilon = 48.0 * Bolz![J] Minimum in the 
                                                 !    Lennard Jones equation
    real(dp), parameter :: sigma = 3.006E-10 ![m]. Constant in the Lennard 
                                            !     Jones equation
    real(dp), parameter :: sigmaSq = sigma**2 ![m^2] Sigma squared
    real(dp), parameter :: timeStep = 4.0E-15 ![s] Time between calculations
    real(dp), parameter :: delR = 0.003E-9 !Distance between shells used for g(r)

    real(dp) :: fourEps = 4.0 * epsilon ![J] Epsilon*4. Used for optimization
    real(dp) :: twentyFourEps = 24.0 * epsilon ![J] Epsilon*24. 
                                               !        Used for optimization

    !Speed Distribution
    integer, parameter :: numVelBox = 80
    integer, parameter :: numMaxwell = 800
    real(dp), parameter :: boxScale = 0.1


    real(dp), parameter :: initTemperature = 320.0 ![K] initial const temperature 
                                                    !of the system
    real(dp), parameter :: finalTemperature = 77.0 ![K] final temperature of 
                                                    !the system

    integer, parameter :: numNVESteps = 50000 !Number of timesteps in the program
    integer, parameter :: numNVTSteps = 200000 !Number of timesteps in the program
    integer, parameter :: numNVTConstSteps = 50000 !Number of timesteps that the 
                                                !system is equilibrated at a 
                                                !constant temperature

    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                   !between momentum-zeroing
    integer, parameter :: numTrajSteps = 100 !Number of timesteps between 
                                             !trajectory outputs to .ttr file
    integer, parameter :: temperatureStep = 205 !Number of timesteps between 
                                                !temperature decrements
end module Constants
