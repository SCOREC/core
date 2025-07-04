from mpi4py import MPI
import pyCore

import sys, getopt

def main(argv):
   model = ''
   mesh  = ''
   try:
      opts, args = getopt.getopt(argv,"hg:m:",["model=","mesh="])
   except getopt.GetoptError:
      print('test_pytCore.py -g <model> -m <mesh>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test_pytCore.py -g <model> -m <mesh>')
         sys.exit()
      elif opt in ("-g", "--model"):
         model = arg
      elif opt in ("-m", "--mesh"):
         mesh = arg
   print('Model file is "', model)
   print('Mesh  file is "', mesh)

   # PCU initialization
   PCUObj = pyCore.PCU(MPI.COMM_WORLD)

   # gmi initialization
   pyCore.gmi_register_mesh()

   # load the mesh and model and write the initial mesh to vtk
   mesh = pyCore.loadMdsMesh(model, mesh, PCUObj)
   pyCore.writeASCIIVtkFiles('before', mesh);

   # setup uniform refiner and call mesh adapt
   ma_input = pyCore.configureUniformRefine(mesh, 2);
   pyCore.adapt(ma_input);

   # write the adapted mesh to vtk
   pyCore.writeASCIIVtkFiles('after', mesh);

   # gmi finalization
   del PCUObj

if __name__ == "__main__":
   main(sys.argv[1:])

