mpif90 program.f90 -o solve
rm result
mpirun -np 4 solve
python2 refine_mesh.py