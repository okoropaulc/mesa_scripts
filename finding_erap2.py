Python 3.7.2 (tags/v3.7.2:9a3ffc0492, Dec 23 2018, 23:09:28) [MSC v.1916 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> import numpy as np
>>> from sklearn.svm import SVR
>>> import pandas as pd
>>> from sklearn.model_selection import cross_val_score
>>> from statistics import mean
>>> import math
>>> from sklearn.linear_model import LinearRegression
>>> from sklearn.preprocessing import scale
>>> from pandas import DataFrame
>>> import pickle
>>> from sklearn.ensemble import RandomForestRegressor
>>> def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")
     #snpanot = snpanot[((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))]
     snpanot = snpanot[(((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))) & (snpanot["rsid"].notna())]
     snpanot = snpanot.drop_duplicates(["varID"])
     return snpanot

>>> def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) & (gene_anot["gene_type"].isin(gene_types))]
     return gene_anot

>>> def get_gene_type (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_type = gene_type.iloc[0,5]
     return gene_type

>>> def get_gene_name (gene_anot, gene):
     gene_name = gene_anot[gene_anot["gene_id"]==gene]
     gene_name = gene_name.iloc[0,2]
     return gene_name

>>> def get_gene_coords (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_coord = [gene_type.iloc[0,3], gene_type.iloc[0,4]]
     return gene_coord

>>> def get_covariates (cov_filepath):
      cov = pd.read_csv(cov_filepath, sep=" ")
      cov = cov.set_index("IID") #make IID to be the row names
      cov.index.names = [None] # remove the iid name from the row
      pc = ["PC1", "PC2", "PC3"] #a list of the PCs to retain
      cov = cov[pc]
      return cov

>>> def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	#print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df

>>> def adjust_for_covariates (expr_vec, cov_df):   
      reg = LinearRegression().fit(cov_df, expr_vec)
      ypred = reg.predict(cov_df)
      residuals = expr_vec - ypred
      residuals = scale(residuals)
      return residuals

>>> def get_maf_filtered_genotype(genotype_file_name,  maf):
	gt_df = pd.read_csv(genotype_file_name, 'r', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

>>> def get_cis_genotype (gt_df, snp_annot, coords, cis_window=1000000):
      snp_info = snpannot[(snpannot['pos'] >= (coords[0] - cis_window)) & (snpannot['rsid'].notna()) & (snpannot['pos'] <= (coords[1] + cis_window))]
      if len(snp_info) == 0:
          return 0
      else:
           gtdf_col = list(gt_df.columns)
           snpinfo_col = list(snp_info["varID"])
           intersect = snps_intersect(gtdf_col, snpinfo_col) #this function was defined earlier
           cis_gt = gt_df[intersect]
           return cis_gt

        
>>> def calc_R2 (y, y_pred):
    tss = 0
    rss = 0
    for i in range(len(y)):
        tss = tss + (y[i])**2
        rss = rss + (((y[i]) - (y_pred[i]))**2)
    tss = float(tss)
    rss = float(rss)
    r2 = 1 - (rss/tss)
    return r2

>>> 
>>> def calc_corr (y, y_pred): #pearson corr
    num = 0
    dem1 = 0
    dem2 = 0
    for i in range(len(y)):
        num = num + ((y[i]) * (y_pred[i]))
        dem1 = dem1 + (y[i])**2
        dem2 = dem2 + (y_pred[i])**2
    num = float(num)
    dem1 = math.sqrt(float(dem1))
    dem2 = math.sqrt(float(dem2))
    rho = num/(dem1*dem2)
    return rho

>>> def snps_intersect(list1, list2):
     return list(set(list1) & set(list2))

>>> afa_snp = "Z:data/mesa_models/AFA_"+str(chrom)+"_snp.txt"
Traceback (most recent call last):
  File "<pyshell#35>", line 1, in <module>
    afa_snp = "Z:data/mesa_models/AFA_"+str(chrom)+"_snp.txt"
NameError: name 'chrom' is not defined
>>> chrom = 5
>>> afa_snp = "Z:data/mesa_models/AFA_"+str(chrom)+"_snp.txt"
>>> gex = "/home/paul/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
>>> gex = "Z:/data/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
>>> afa_snp = "Z:/data/mesa_models/AFA_"+str(chrom)+"_snp.txt"
>>> geneanotfile = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt"
>>> snpfilepath = "Z:/data/mesa_models/AFA_"+str(chrom)+"_annot.txt"
>>> 
>>> test_snp = "Z:/data/METS_model/hg19/METS_"+str(chrom)+"_snp.txt"
>>> test_gex = "Z:/data/METS_model/hg19/METS_peer10_all_chr_prot_coding_gex.txt"
>>> test_covfile = "Z:/data/METS_model/hg19/METS_3_PCs.txt"
>>> test_snpfile = "Z:/data/METS_model/hg19/METS_"+str(chrom)+"_annot.txt"
>>> gencodev28 = "Z:/data/METS_model/hg19/gencode.v28_annotation.parsed.txt"
>>> snpannot = get_filtered_snp_annot(snpfilepath)
>>> geneannot = get_gene_annotation(geneanotfile, chrom)
>>> expr_df = get_gene_expression(gex, geneannot)
>>> annot_geneid = geneannot["gene_id"]
>>> annot_geneid = list(annot_geneid)
>>> agid = []
>>> for i in annot_geneid:
	agid.append(i[0:(i.find("."))])

	
>>> geneannot["gene_id"] = agid #replace with non decimal gene_id
>>> cov = get_covariates(cov_file)
Traceback (most recent call last):
  File "<pyshell#58>", line 1, in <module>
    cov = get_covariates(cov_file)
NameError: name 'cov_file' is not defined
>>> cov_file = "Z:/data/mesa_models/AFA_3_PCs.txt"
>>> cov = get_covariates(cov_file)
>>> geneannot.head(3)
      chr          gene_id gene_name   start     end       gene_type
14877   5  ENSG00000153404  PLEKHG4B  140373  190085  protein_coding
14879   5  ENSG00000185028   LRRC14B  191626  195468  protein_coding
14881   5  ENSG00000164366   CCDC127  196986  218330  protein_coding
>>> geneannot[geneannot['gene_name']=='ERRAP']
Empty DataFrame
Columns: [chr, gene_id, gene_name, start, end, gene_type]
Index: []
>>> geneannot[geneannot['gene_name']=='ERAP2']
      chr          gene_id gene_name     start       end       gene_type
16255   5  ENSG00000164308     ERAP2  96211643  96255420  protein_coding
>>> genes = list(expr_df.columns)
>>> gt_df = get_maf_filtered_genotype(afa_snp, 0.01)
>>> train_ids = list(gt_df.index)
>>> train_g = []
>>> for i in genes:
	train_g.append(i[0:(i.find("."))])

	
>>> expr_df.columns = train_g
>>> genes = list(expr_df.columns) #
>>> expr_df.head(3)
       ENSG00000172262  ENSG00000055163  ...  ENSG00000113638  ENSG00000113595
24779        -0.441862        -0.281812  ...        -0.354059        -0.287153
24795         0.071825         0.120564  ...         0.031082         0.295182
24824         0.008064        -0.180102  ...         0.098433         0.025786

[3 rows x 431 columns]
>>> test_snpannot = get_filtered_snp_annot(test_snpfile)
>>> test_geneannot = get_gene_annotation(gencodev28, chrom)
>>> test_expr_df = get_gene_expression(test_gex, test_geneannot)
>>> test_annot_geneid = test_geneannot["gene_id"]
>>> test_annot_geneid = list(test_annot_geneid)
>>> test_agid = []
>>> for i in test_annot_geneid:
	test_agid.append(i[0:(i.find("."))])

	
>>> test_geneannot["gene_id"] = test_agid
>>> train_geneannot[train_geneannot['gene_name']=='ERAP2']
Traceback (most recent call last):
  File "<pyshell#82>", line 1, in <module>
    train_geneannot[train_geneannot['gene_name']=='ERAP2']
NameError: name 'train_geneannot' is not defined
>>> test_geneannot[test_geneannot['gene_name']=='ERAP2']
      chr          gene_id gene_name     start       end       gene_type
16187   5  ENSG00000164308     ERAP2  96875939  96919716  protein_coding
>>> test_cov = get_covariates(test_covfile)
>>> test_genes = list(test_expr_df.columns)
>>> test_gt_df = get_maf_filtered_genotype(test_snp, 0.01)
>>> test_ids = list(test_gt_df.index)
>>> test_g = []
>>> for i in test_genes:
	test_g.append(i[0:(i.find("."))])

	
>>> test_expr_df.columns = test_g
>>> test_genes = list(test_expr_df.columns)
>>> test_expr_df.head(3)
         ENSG00000120705  ENSG00000170571  ...  ENSG00000145692  ENSG00000164244
LD1_LD1        -2.051153        -2.289820  ...         7.626484        -7.850946
LD2_LD2        -3.688696         0.199755  ...        -4.348942         1.735109
LD3_LD3         9.168135         0.657448  ...        -5.601410         3.169443

[3 rows x 735 columns]
>>> rf = RandomForestRegressor(max_depth=None, random_state=1234, n_estimators=100)
>>> svr = SVR(kernel="rbf", gamma="auto")
>>> gene = "ENSG00000164308"
>>> coords = get_gene_coords(geneannot, gene)
>>> test_coords = get_gene_coords(test_geneannot, gene)
>>> gene_name = get_gene_name(geneannot, gene)
>>> gene_name
'ERAP2'
>>> expr_vec = expr_df[gene]
>>> len(expr_vec)
233
>>> test_expr_vec = test_expr_df[gene]
>>> len(test_expr_vec)
77
>>> adj_exp = adjust_for_covariates(list(expr_vec), cov)
>>> test_adj_exp = adjust_for_covariates(list(test_expr_vec), test_cov)
>>> cis_gt = get_cis_genotype(gt_df, snpannot, coords)
>>> cis_gt.head(3)
id     5_97050559_G_A_b37  ...  5_96659728_T_C_b37
24779                 0.0  ...                 1.0
24795                 0.0  ...                 1.0
24824                 0.0  ...                 2.0

[3 rows x 8157 columns]
>>> test_cis_gt = get_cis_genotype(test_gt_df, test_snpannot, test_coords)
>>> test_cis_gt.head(3)
id       5_97050559_G_A_b37  ...  5_96971832_T_C_b37
LD1_LD1                 0.0  ...                 0.0
LD2_LD2                 0.0  ...                 0.0
LD3_LD3                 0.0  ...                 0.0

[3 rows x 7066 columns]
>>> train_snps = list(cis_gt.columns)
>>> test_snps = list(test_cis_gt.columns)
>>> snp_intersect = snps_intersect(train_snps, test_snps)
>>> cis_gt = cis_gt[snp_intersect]
>>> cis_gt.head(3)
id     5_97050559_G_A_b37  ...  5_96647141_A_G_b37
24779                 0.0  ...                 2.0
24795                 0.0  ...                 1.0
24824                 0.0  ...                 2.0

[3 rows x 4774 columns]
>>> test_cis_gt = test_cis_gt[snp_intersect]
>>> test_cis_gt(3)
Traceback (most recent call last):
  File "<pyshell#117>", line 1, in <module>
    test_cis_gt(3)
TypeError: 'DataFrame' object is not callable
>>> test_cis_gt.head(3)
id       5_97050559_G_A_b37  ...  5_96647141_A_G_b37
LD1_LD1                 0.0  ...                 2.0
LD2_LD2                 0.0  ...                 1.0
LD3_LD3                 0.0  ...                 2.0

[3 rows x 4774 columns]
>>> cis_gt = cis_gt.values
>>> test_cis_gt = test_cis_gt.values
>>> test_yobs = test_expr_vec.values
>>> rf.fit(cis_gt, adj_exp.ravel())
RandomForestRegressor(bootstrap=True, criterion='mse', max_depth=None,
           max_features='auto', max_leaf_nodes=None,
           min_impurity_decrease=0.0, min_impurity_split=None,
           min_samples_leaf=1, min_samples_split=2,
           min_weight_fraction_leaf=0.0, n_estimators=100, n_jobs=None,
           oob_score=False, random_state=1234, verbose=0, warm_start=False)
>>> ypred = rf.predict(test_cis_gt)
>>> stats.pearsonr(test_yobs, ypred)
Traceback (most recent call last):
  File "<pyshell#124>", line 1, in <module>
    stats.pearsonr(test_yobs, ypred)
NameError: name 'stats' is not defined
>>> from scipy import stats
>>> stats.pearsonr(test_yobs, ypred)
(0.7673698948281062, 3.9705761505773547e-16)
>>> svr.fit(cis_gt, adj_exp.ravel())
SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.1, gamma='auto',
  kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
>>> ypred = svr.predict(test_cis_gt)
>>> stats.pearsonr(test_yobs, ypred)
(0.8268051071099968, 2.0022068186125852e-20)
>>> stats.spearmanr(test_yobs, ypred)
SpearmanrResult(correlation=0.8181794282665455, pvalue=1.0479992021147932e-19)
>>> rf.fit(cis_gt, adj_exp.ravel())
RandomForestRegressor(bootstrap=True, criterion='mse', max_depth=None,
           max_features='auto', max_leaf_nodes=None,
           min_impurity_decrease=0.0, min_impurity_split=None,
           min_samples_leaf=1, min_samples_split=2,
           min_weight_fraction_leaf=0.0, n_estimators=100, n_jobs=None,
           oob_score=False, random_state=1234, verbose=0, warm_start=False)
>>> ypred = rf.predict(test_cis_gt)
>>> stats.spearmanr(test_yobs, ypred)
SpearmanrResult(correlation=0.7986197831777706, pvalue=3.2986058917975638e-18)
>>> 
