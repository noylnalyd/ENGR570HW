#!/bin/bash
#cmake ./
make


filePrefixes=("nnc666" "s1rmq4m1" "impcol_e" "fs_183_1" "lns_3937")
fileNs=(666 5489 225 183 3937)
for matrix in {0..4}; do
    n=${fileNs[$matrix]}
    filePrefix=${filePrefixes[$matrix]}
    # python3 makeVecFile.py $filePrefix $n
    echo $filePrefix
    for SpFmt in "COO" "DEN" "CSR" #"JDS" "DEN" "COO" "CSR" "ELL"
    do
        echo $SpFmt
        ./SpMV.exe ${SpFmt} 10 ${filePrefix}.mtx ${filePrefix}IN.txt ${filePrefix}${SpFmt}OUT.txt
    done
done