% from www.nature.com/nmeth/journal/v13/n10/extref/nmeth.3971-S2.zip
clear;
dirRoot = "your_project_root_dir"
load(dirRoot+"raw/nmeth.3971-S2.mat");

csvwrite(dirRoot+"raw/klein15_matrix_CellCylceCorrFromDPT.csv", datacc);

formatSpec = '%s\n';
fileIDgenes = fopen(dirRoot+"raw/klein15_geneIds_CellCylceCorrFromDPT.csv",'w');
ngenes = size( cc_genes );
for gene = 1:ngenes
    fprintf(fileIDgenes,formatSpec,cc_genes{gene});
end
fclose(fileIDgenes);

csvwrite(dirRoot+"raw/klein15_cellLabels_CellCylceCorrFromDPT.csv", labels);
