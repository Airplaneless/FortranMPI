#!/bin/bash
gfortran tests/split_matrix.f90 -o test_splitmatrix
./test_splitmatrix > test_splitmatrix.log