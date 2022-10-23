#!/bin/bash
#cmake ./
make laplace

filePrefixes=("JI")
fileNs=(8 32 64 128 256)
for solver in {0..0}; do
    filePrefix=${filePrefixes[$solver]}
    for n in {0..4}; do
        echo $filePrefix $n
        ./laplace.exe ${filePrefix} $n
    done
done
