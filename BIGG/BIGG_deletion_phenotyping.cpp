#include "enet.h"
#include "erandomwalk.h"
#include <eutils/logger.h>
#include <eutils/emain.h>
#include <signal.h>
#include <eutils/eheap.h>
#include <algorithm>
#include <eutils/ernd.h>
#include <eutils/estrarrayof.h>
#include <eutils/eregexp.h>
#include <vector>

using std::vector;
using namespace std;

// Function "deleted" will output a vec2 vector, whose non-zero elements corrspond to reactions to be deleted.

vector<int>  deleted(int unique_reaction_set[], int reaction_indices[], int grrule[][30],int gen_del[], const int num000, const int num00, const int del_size){

             int vec1[num00];
             vector<int> vec2(num000);
             for (int i=0;i<num00;i++){vec1[i]=1;}
             for (int i=0;i<num00;i++){
                  int alaki=0;
                  for (int j=0;j<30;j++){
                       if (grrule[i][j]==0){alaki++;}
                       for (int k=0;k<del_size;k++){
                           if (grrule[i][j]==gen_del[k]) {vec1[i]=0;break;}
                       }
                       if (vec1[i]==0){break;}
                       if (alaki>1){break;}
                  }
             }
             for (int i=0;i<num000;i++){
                  int count=0;
                  int flag=0;
                  for (int j=0;j<num00;j++){
                      if (unique_reaction_set[i]==reaction_indices[j]){count=count+vec1[j];flag=1;}
                      else {if (flag==1) {break;}}
                  }
                  if (count==0){vec2[i]=unique_reaction_set[i];}
                  else {vec2[i]=0;}
             }
             return vec2;           
}
estr solver="esolver_clp";
int netsize=-1;
int strict=0;
int periphery_only=0;
int mutate_transport=0;
int internal_secretion=0;
int only_viable=0;
double gendist=0.0;
enet net;
erandomWalk *prw=0x00;
int num1=-1;//length of grrule.dat file
int num2=-1;//length of unique_reaction_set.dat file
int num3=-1;//index of the carbon source

int emain()
{ ldieif(argvc<4,"syntax: ./BIGG_deletion_phenotyping  <universe.net> <genome.dat> <deletion.dat> <grrule.dat> <unique_reaction_set.dat> <output.dat> --num1 --num2 --num3 <fluxbounds.flx>"); 
////////////// ARGUMENTS /////////////////////////
  epregister(num1);
  epregister(num2);
  epregister(num3);
  eparseArgs(argvc,argv);

////////////// Universe loading //////////////////
  net.load(argv[1]); 
  net.correct_malformed();
  erandomWalk rw(net,solver,strict);
  prw=&rw; 
  rw.periphery_only=periphery_only;
  rw.mutate_transport=mutate_transport;
  rw.internal_secretion=internal_secretion;
  rw.only_viable=only_viable;
  rw.setRSize(netsize);
  rw.getEnv(argvc,argv);
  rw.load(net);

////////////////// Open files //////////////////////////////
  efile file_genome;
  file_genome.open(argv[2],"r");
  efile file_del;
  file_del.open(argv[3],"r");
  efile file_grrule;
  file_grrule.open(argv[4],"r");
  efile file_unique_reaction_set;
  file_unique_reaction_set.open(argv[5],"r");
  efile file_out;
  file_out.open(argv[6],"a");

/////////// Reading genome.dat file ////////////////////////
  estr str;
  estrarray parts;
  
  eintarray gen;
  int temp=0;

  while (file_genome.readln(str)) {
  	 parts=str.explode(" ");
         for (int i=0;i<5927;i++){gen.add(temp);}
         for (int k=0;k<5927;k++){gen[k]=parts[k].i();}
  }
  file_genome.close();
  eintarray gen2=gen;

/////////// Reading grrule.dat file ////////////////////////
   int grrule [num1][30];
   int reaction_indices [num1];
   for (int l=0;l<num1;l++){
       for (int p=0;p<30;p++){grrule[l][p]=0;}
       reaction_indices[l]=0;
   }

   int count=0;
   while (file_grrule.readln(str)){
          parts=str.explode(" ");
          int ss=parts.size();
          for (int k=1;k<ss;k++){grrule[count][k]=parts[k].i();}
          reaction_indices[count]=parts[0].i();
          count++;
   }
   file_grrule.close();
   for (int k=0;k<num1;k++){
        eintarray alaki;
        for (int kk=0;kk<30;kk++){
             alaki.add(grrule[k][kk]);
        }
        //cout<<intarr2str2(alaki)<<endl;
   }
/////////// Reading unique_reaction_set.dat file ////////////////////////
   int unique_reaction_set[num2];
   count=0;
   while (file_unique_reaction_set.readln(str)){
          parts=str.explode(" ");
          unique_reaction_set[count]=parts[0].i();
          count++;
   }
   file_unique_reaction_set.close();

////////////////////// Reading deletion.dat file ////////////////////////
   //int ccount=0;
   while (file_del.readln(str)){
          parts=str.explode(" ");
          int del_size=parts.size();
          int gen_del[del_size];
          for (int g=0;g<parts.size();g++){gen_del[g]=0;}              
          for (int g=0;g<parts.size();g++){gen_del[g]=parts[g].i();}              
          // in the line below the set of deleted reactions based on the set of deleted genes is determined.
          vector<int> rxn_del=deleted(unique_reaction_set,reaction_indices,grrule,gen_del,num2,num1,del_size);
          gen2=gen;
          for (int k=0;k<num2;k++){if (rxn_del[k]>0){gen2[(rxn_del[k]-1)]=0;}}// the genotype (gen2) is updated by excluding the deleted reactions
          for (int m=0; m<5927;m++){if (gen2[m]==0){rw.disable(m);}}
          rw.calcPhenotype(); //phenotype is determined
          eintarray phen = rw.phenotype;
	  for (int n=0; n<5927; n++) { rw.activate(n);}
          estr intstr = intarr2str2(phen);
          file_out.write(intstr+"\n");//phenotype vector is saved

   }
   file_del.close(); 
   file_out.close();        
   return(0);
}

