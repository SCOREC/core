from mpi4py import MPI
import pyCore

import sys, getopt

def main(argv):
   model = ''
   mesh  = ''
   try:
      opts, args = getopt.getopt(argv,"hg:m:",["model=","mesh="])
   except getopt.GetoptError:
      print 'test_pytCore.py -g <model> -m <mesh>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test_pytCore.py -g <model> -m <mesh>'
         sys.exit()
      elif opt in ("-g", "--model"):
         model = arg
      elif opt in ("-m", "--mesh"):
         mesh = arg
   print 'Model file is "', model
   print 'Mesh  file is "', mesh

   # PCU initialization
   pyCore.PCU_Comm_Init()

   # SIMX initialization
   pyCore.start_sim('simlog.txt')

   # gmi initialization
   pyCore.gmi_register_mesh()

   # gmi_sim start
   pyCore.gmi_sim_start()
   pyCore.gmi_register_sim()

   # load the mesh and model and write the initial mesh to vtk
   mesh = pyCore.loadMdsMesh(model, mesh)
   print("num verts in the mesh is ", mesh.count(0))

   it = mesh.begin(0)
   while True:
     e = mesh.iterate(it)
     if (not e):
       break
     p = pyCore.Vector3()
     mesh.getPoint(e, 0, p)
     print(p.x(), p.y(), p.z())
     print(mesh.getModelTag(mesh.toModel(e)))
   mesh.end(it)



   pyCore.writeASCIIVtkFiles('before', mesh);

   # setup uniform refiner and call mesh adapt
   ma_input = pyCore.configureUniformRefine(mesh, 2);
   pyCore.adapt(ma_input);

   it = mesh.begin(0)
   while True:
     e = mesh.iterate(it)
     if (not e):
       break
     p = pyCore.Vector3()
     mesh.getPoint(e, 0, p)
     print(p.x(), p.y(), p.z())
     print(mesh.getModelTag(mesh.toModel(e)))
   mesh.end(it)




   # write the adapted mesh to vtk
   pyCore.writeASCIIVtkFiles('after', mesh);


   # gmi_sim stop
   pyCore.gmi_sim_stop()

   # SIMX finalization
   pyCore.stop_sim()

   # gmi finalization
   pyCore.PCU_Comm_Free()

if __name__ == "__main__":
   main(sys.argv[1:])

