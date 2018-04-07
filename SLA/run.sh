#!/bin/bash
mpif90 program.f90 -o solve
mpirun -np 2 solve > debug.log