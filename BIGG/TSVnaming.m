function TSVnaming(num)
%%
nami1=sprintf('/home/rzgar/BIGG/OLD/TSV2_name_%d.dat',num);
nami2=sprintf('/home/rzgar/BIGG/OLD/FINAL_comp_%d.dat',num);
nami5=sprintf('/home/rzgar/BIGG/OLD/FINAL_REV_%d.dat',num);
namo=sprintf('/home/rzgar/BIGG/OLD/FINAL_name_%d.dat',num);
fclose('all');
count1=0;count2=0;count3=0;count4=0;count5=0;
fid1=fopen(nami1,'r');
fid2=fopen(nami2,'r');
fid5=fopen(nami5,'r');

while feof(fid1)==0
    fgetl(fid1);count1=count1+1;
end
while feof(fid2)==0
    fgetl(fid2);count2=count2+1;
end
while feof(fid5)==0
    fgetl(fid5);count5=count5+1;
end
fclose('all');


cel1=cell(count1,1);cel2=cell(count2,1);cel5=cell(count5,1);
count1=0;count2=0;count5=0;
fid1=fopen(nami1,'r');
fid2=fopen(nami2,'r');
fid5=fopen(nami5,'r');

while feof(fid1)==0
    count1=count1+1;cel1{count1,1}=fgetl(fid1);
end
while feof(fid2)==0
    count2=count2+1;cel2{count2,1}=fgetl(fid2);
end
while feof(fid5)==0
    count5=count5+1;cel5{count5,1}=fgetl(fid5);
end
fclose('all');

fod=fopen(namo,'a');
for i=1:count1
    string=cel1{i,1};
    if (length(strfind(cel2{i,1},'e'))+length(strfind(cel2{i,1},'p')))==0
        if (strcmp(cel5{i,1},'Reversible')==1)
             final=sprintf('RRR_%s_RV',string);
        else
             final=sprintf('RRR_%s_IRV',string);
        end
    else
        final=string;
    end
    fprintf(fod,'%s\n',final);
end
fclose('all');
end
