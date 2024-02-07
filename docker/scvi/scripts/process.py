import scanpy
import argparse
import pandas as pd

# Create the parser
parser = argparse.ArgumentParser(description='Normalize seurat objects')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')

parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
parser.add_argument('--adata-output', dest='adata_output', type=str, 
                    help='Output file to save AnnData object to')

parser.add_argument('--top-genes', dest='top_genes', type=str, 
                    help='csv containing top genes', default='top_genes.csv')


# Parse the arguments
args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# TODO: impliment cell cycle scoring 
# https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt

# cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
# # Split into 2 lists
# s_genes = cell_cycle_genes[:43]
# g2m_genes = cell_cycle_genes[43:]

# cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
# TODO load the top_genes from the qc plotting step and subset... to the top 8k genes
#  if we have memory issues, consider subsetting to a gene_list
# gene_list = pd.read_csv(top_genes, header=None)[0].tolist()
# adata = adata[:, gene_list.index]

# does this work with sparse uint8?
adata.layers['counts'] = adata.X.copy() # type: ignore

scanpy.pp.normalize_total(adata, target_sum=1e4)
scanpy.pp.log1p(adata)

# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
scanpy.pp.highly_variable_genes(
    adata, 
    batch_key='sample', 
    subset=True, 
    flavor='seurat_v3', 
    layer='counts', 
    n_top_genes=5000
)

adata.write_h5ad(filename=args.adata_output, compression='gzip') 