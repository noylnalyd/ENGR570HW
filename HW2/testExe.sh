#!/bin/bash
cmake ./
make

./SpMV.exe COO 2 lns_3937.mtx ex3937IN.txt ex3937OUT.txt