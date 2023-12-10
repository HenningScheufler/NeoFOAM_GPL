#include <t8.h>
#include "NeoFOAM/mesh/treeAMR/treeAMRMesh.hpp"
#include "NeoFOAM_GPL/treeAMR/t8_MeshModifier.hpp"
#include <memory>

int main (int argc, char **argv)
{
  int mpiret;
  /* The prefix for our output files. */
  // const char prefix[BUFSIZ] = "t8_step1_tetcube";
  // t8_locidx_t local_num_trees;
  // t8_gloidx_t global_num_trees;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  {
    // std::unique_ptr<NeoFOAM::treeAMRMeshModifier> t8mod = std::make_unique<t8_MeshModifier>(comm,5,8);
    NeoFOAM::treeAMRMesh tMesh(std::make_unique<t8_MeshModifier>(comm,5,8));

    tMesh.refine();

    tMesh.write();

    std::cout << "nElements: " << tMesh.nElements() << std::endl;
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
