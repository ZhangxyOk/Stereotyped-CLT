# Stereotyped-CLT
The analysis of high-density lineage tree of lung progenitor cell differentiation model in vitro suggests that stable cell type composition among sub-CLTs and stereotypical development program support development robustness.
# about main figures 
1. all codes could be found according to the file name for each main figure.
2. all datas used by code could be found according to the file name for each main figure in the file folder "data".
# about mDELTA alignment 
1. Our team has developed mDELTA, a modified DELTA algorithm that can handle CLTs with multifurcating internal nodes and quantitatively labeled (i.e. with single-cell transcriptomes) terminal cells. Below are the algorithmic details of mDELTA. you could learn and used it (https://github.com/Chenjy0212/mdelta). 
2. For all lung progenitor differentiation samples,the parameters of mDELTA local alignment seperately are: top100(sub_CLT one-way),prune score:0.2,prune_percent:20%,merge score:100, diff:90, permutation:1000.
3. All input datas of mDELTA alignment of sample pair are located in the data path "data/fig5/mdelta_align_input".
4. The quantitatively score data of cell-cell type pair "df_euclidean_qualitative_30" was came from the input data "data/fig5/mdelta_prepare" and calculated with the python script "mDELTA_prepro"(https://github.com/Chenjy0212/mdelta_prepro) when the "cutoff" as "15 percentile and 55 percentile". 

