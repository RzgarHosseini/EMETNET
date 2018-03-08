# First Step:
Please see the README file in https://github.com/RzgarHosseini/EMETNET/blob/master/README.md to see how to install the software for flux balance analysis that is required to run our "BIGG_deletion_phenotyping.cpp" program.

# BIGG_deletion_phenotyping.cpp:
This is the program that I used to compute the deletional robustness of bacterial genomes.
The syntax of the program is as follows:
BIGG_deletion_phenotyping  <universe.net> <genome.dat> <deletion.dat> <grrule.dat> <unique_reaction_set.dat> <output.dat> --num1 --num2  <fluxbounds.flx>
where:
## <universe.net>: 
is the universe of possible reactions in 55 prokaryotic genomes that we have considered. You can see an example in this folder : DATA/universe46.net
## <genome.dat>: 
is a binary vector whose length is equal to the number of reactions in the universe. Each element of this vecor is 1 if the corresponding reaction is present in the given genome and 0 otherwise.You can see an example in this folder for E.coli K12 MG1655: DATA/genome46.dat
## <deletion.dat>: 
is a file in each line of which the set of genes that are deleted is specified. The type of deletions is specified in this file. For example, in tandem deletion of n=3, each line contains 3 consecutive numbers, while in random deletions of n=3, each line contains 3 randomly chosen numbers between 1:N(number of metabolic genes in the given genome). You can see two examples in this folder: DATA/tandem_deletion3.dat and DATA/random_deletion3.dat.    
## <grrule.dat>:
In this file gene-reaction association rules are specified. The first number in each line is the index of the reaction specified in the reaction universe. The other numbers in each line are the set of genes encoding that enzyme(s) catalyzing that reaction. The logical AND applies to the set of genes in a given line. There might be multiple lines corresponding to a given reaction (i.e., multiple lines with the same first number). In this case logical OR applies to the different lines. You can see an example in this folder for E.coli K12 MG1655: DATA/grRule46.dat.
## <unique_reaction_set.dat>:
This file contains the unique set of reaction indices corresping to a given metabolism. You can download an example corresponding to E.coli K12 MG1655 in this folder (DATA/Reaction_set46.dat). 
## <output.dat>:
Is the output file, each line of which corresponds to a given deletional variant specified in the <deletion.dat> file. Each line of this file is a binary vector of length L (i.e. the number of carbon sources we have used in our study).
## num1
is the number of lines in the <grrule.dat> file
## num2
is the number of lines in <unique_reaction_set.dat> file
## <fluxbounds.flx>
is the name of the file(s) corresponding to minimal environment(s) with a given carbon source. In the subfolder ENVS you can see 101 different minimal environments with distinct carbon sources.

# Example:
In this example, we run this program twice to determine the tandem and random robustness to deletion of 3 metabolic genes in E.coli K12 MG1655:

i) tandem deletion:
BIGG_deletion_phenotyping  DATA/universe46.net DATA/genome46.dat DATA/tandem_deletion3.dat DATA/grrule46.dat DATA/Reaction_set46.dat DATA/tandem3_phenotype.dat --num1 3900 --num2 2572 ENVS/*

ii) random deletion:
BIGG_deletion_phenotyping  DATA/universe46.net DATA/genome46.dat DATA/random_deletion3.dat DATA/grrule46.dat DATA/Reaction_set46.dat DATA/random3_phenotype.dat --num1 3900 --num2 2572 ENVS/*

To quantify robustness to tandem deletion of n=3, we check what fraction of the deletional variants (each line of the output file (tandem3_phenotype.dat)) retains viability on all carbon sources (strict definition) or at least one carbon source (moderate definition) on which the wild-type genome is viable. Similarly, to compute robustness to random deletion of n=3 we can use random3_phenotype.dat.  














