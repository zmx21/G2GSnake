#!/usr/bin/env python
import pandas as pd
import numpy as np
#Read pPCA file
ppca_file = pd.read_csv("../results/tmp/pPCA.txt",delim_whitespace=True)
#Read PCA file, rename headers
pca_file = pd.read_csv("../results/tmp/G2G_QC.eigenvec",delim_whitespace=True,header=None)        
pca_file = pca_file.iloc[:,1:len(pca_file.columns)]
pca_file.columns = ["#IID"] + ["PC" + str(x) for x in np.arange(len(pca_file.columns)-1)+1]
#Read provided covar file
covar_file = pd.read_table(snakemake.input[2])
#Join all files
merged_covar_file = covar_file.merge(ppca_file,on='PID')
merged_covar_file = merged_covar_file.merge(pca_file,on='#IID')
merged_covar_file = merged_covar_file.drop(columns=['PID'])
#Write to csv
merged_covar_file.to_csv("../results/tmp/merged_covar.txt",sep = ' ',index=False)
#Join PID and IID for AA table
aa_tbl = pd.read_table(snakemake.input[3])
map_file = covar_file[["#IID","PID"]]
aa_tbl_jned = aa_tbl.merge(map_file,left_index=True,right_on='PID')
aa_tbl_jned.to_csv("../results/tmp/aa_tbl_joined.txt",sep='\t',na_rep='NA',index=False)