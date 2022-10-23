#!/bin/bash
#cmake ./
make laplace

filePrefixes=("JI" "GS" "RB")
fileNs=(8 32 64 128 256)
for solver in {0..2}; do
    filePrefix=${filePrefixes[$solver]}
    for n in {0..2}; do
        echo $filePrefix $n
        ./laplace.exe ${filePrefix} ${fileNs[$n]}
    done
done
