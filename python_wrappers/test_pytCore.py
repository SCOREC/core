from mpi4py import MPI
import pyCore

# PCU initialization
pyCore.PCU_Comm_Init()

# gmi initialization
pyCore.gmi_register_mesh()

# load the mesh and model and write the initial mesh to vtk
mesh = pyCore.loadMdsMesh('../cube.dmg', '../cube.smb')
pyCore.writeASCIIVtkFiles('before', mesh);

# setup uniform refiner and call mesh adapt
ma_input = pyCore.configureUniformRefine(mesh, 2);
pyCore.adapt(ma_input);

# write the adapted mesh to vtk
pyCore.writeASCIIVtkFiles('after', mesh);

# gmi finalization
pyCore.PCU_Comm_Free()
