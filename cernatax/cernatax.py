import os
import pandas as pd


class CERNATAX():
    
    def __init__(self):
                
        wdr = os.path.dirname(os.path.abspath(__file__))        
        self.ref_db = pd.read_csv(os.path.join(wdr, 'db', 'ceRNA_database.csv'), index_col=0)
        
    
    def summarize_ref_db(self):
        ref_db = self.ref_db
        
        print('A total of {} miRNA-mRNA interaction and {} miRNA-lncRNA interaction'.format(
            ref_db[ref_db.type == 'miRNA-mRNA'].shape[0],
            ref_db[ref_db.type == 'miRNA-lncRNA'].shape[0],
        ))
        print(ref_db.type.value_counts())
        
        
    def find_ceRNA_axis_by_DEG(self, deg_df):
        ref_db = self.ref_db
        
        deg_miRNA = deg_df[deg_df.type == 'miRNA'].gene.unique()
        deg_ceRNA = deg_df[deg_df.type != 'miRNA'].gene.unique()
    
        df = ref_db[ref_db.miRNA.isin(deg_miRNA) & ref_db.ceRNA.isin(deg_ceRNA)]
        df = pd.merge(df, deg_df['log2FC'], left_on='miRNA', right_index=True)
        df = pd.merge(df, deg_df['log2FC'], left_on='ceRNA', right_index=True)
        df.columns = ['miRNA', 'ceRNA', 'species', 'database', 'type', 'miRNA_log2FC', 'ceRNA_log2FC']
        df['inference'] = 'strict'
        ceRNA_df = df[(df.miRNA_log2FC * df.ceRNA_log2FC < 0)]
    
        tmp = ceRNA_df.groupby('miRNA').agg({'type': lambda x: len(set(x))})
        axis_df = ceRNA_df[ceRNA_df.miRNA.isin(tmp[tmp.type > 1].index)]
    
        self.ceRNA_df = ceRNA_df
        self.axis_df = axis_df

        return ceRNA_df, axis_df        
    
    def expand_ceRNA_axis_by_loose_DEG(self, deg_strict_df, deg_df):
        # just expand ceRNA based on miRNA axis
        ref_db = self.ref_db
    
        axis_miRNA_list = self.axis_df.miRNA.unique()   
        deg_ceRNA = deg_df[deg_df.type != 'miRNA'].gene.unique()
        
        df = ref_db[ref_db.miRNA.isin(axis_miRNA_list) & ref_db.ceRNA.isin(deg_ceRNA)]
        df = pd.merge(df, deg_strict_df['log2FC'], left_on='miRNA', right_index=True)
        df = pd.merge(df, deg_df['log2FC'], left_on='ceRNA', right_index=True)
        df.columns = ['miRNA', 'ceRNA', 'species', 'database', 'type', 'miRNA_log2FC', 'ceRNA_log2FC']
        df['inference'] = 'loose'
        return df