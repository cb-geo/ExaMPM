/* Main file to run Refactored ExaMPM */

#include <ExaMPM_ProblemManager.hpp>

// #include <ExaMPM_iowriter.hpp>

#include <Cabana_Core.hpp>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <array>
#include <cmath>
//#include <iostream>

//---------------------------------------------------------

/*int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);

  Kokkos::initialize(argc, arg);

 try {

    // Read input json and create IO object 
    auto io = std::make_shared<ExaMPM::IO>(argc, argv);

    // Initialize Problem Manager

    // Get analysis type
    const std::string analysis = io->analysis_type();

    // Create analysis (Initialize Solver)
    auto ExaMPM = 

    //Solve (Run the problem).                                                              
    ExaMPM->solve();
  
    }

  Kokkos::finalize();

  MPI_Finalize();

}  */
  
//-----------------------------------------------------------------------
