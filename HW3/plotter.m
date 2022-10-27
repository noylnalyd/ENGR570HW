myStabber = fopen('myStab.txt','r')
petscStabber = fopen('ksp15.txt','r')

myStab = fscanf(myStabber,'%e');
petscStab = fscanf(petscStabber,'%e');

%%
semilogy(1:length(myStab),myStab,1:length(petscStab),petscStab)
hold on
legend("Dylan","PETSC")
title("BiCGStab Residuals")
xlabel("Iteration #")
ylabel("Residual L2 Norm (Log)")