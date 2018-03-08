# First Step:
Please see the README file in https://github.com/RzgarHosseini/EMETNET/blob/master/README.md to see how to install the software for flux balance analysis that is required to run our "BIGG_deletion_phenotyping.cpp" program.

# BIGG_deletion_phenotyping.cpp:
This is the program that I used to compute the deletional robustness of bacterial genomes.
The syntax of the program is as follows:
BIGG_deletion_phenotyping  <universe.net> <genome.dat> <deletion.dat> <grrule.dat> <unique_reaction_set.dat> <output.dat> --num1 --num2 --num3 <fluxbounds.flx>
where:
## <universe.net>: 
is the universe of possible reactions in 55 prokaryotic genomes that we have considered. You can see an example in this folder : universe46.net
## <genome.dat>: 
is a binary vector whose length is equal to the number of reactions in the universe. Each element of this vecor is 1 if the corresponding reaction is present in the given genome and 0 otherwise.You can see an example in this folder for E.coli K12 MG1655: genome46.dat
## <deletion.dat>: 
is a file in each line of which the set of genes that are deleted is specified.  
## 

