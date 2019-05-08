from mpi4py import MPI
import pyCore

# print MPI.COMM_WORLD.Get_rank()

pyCore.PCU_Comm_Init()
print pyCore.PCU_Time()
print pyCore.PCU_Comm_Self()


pyCore.gmi_register_mesh()

mesh = pyCore.loadMdsMesh('../cube.dmg', '../cube.smb')

ma_input = pyCore.configureUniformRefine(mesh, 2);

pyCore.PCU_Comm_Free()
