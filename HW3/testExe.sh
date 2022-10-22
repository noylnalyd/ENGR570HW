#!/bin/bash
#cmake ./
make qr

filePrefixes=("nnc666" "impcol_e")
for matrix in {0..1}; do
    filePrefix=${filePrefixes[$matrix]}
    echo $filePrefix
    ./qr.exe ${filePrefix}.mtx
done