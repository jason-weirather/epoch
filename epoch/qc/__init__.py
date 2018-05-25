from MulticoreTSNE import MulticoreTSNE as TSNE
from epoch import TCGA, GTEx
from iobio.explicitsemanticanalysis import ReferenceCorrelation
import pandas as pd
import numpy as np
import sys
class CohortQC:
    def __init__(self,expression_matrix,filter=None,limit=None,min_variance=1,verbose=False,use_TCGA=True,use_GTEx=True):
        self._exprs = expression_matrix
        self._verbose = verbose
        reference = []
        if use_TCGA:
            if verbose: sys.stderr.write("loading TCGA\n")
            tcga = TCGA.get_expression()
            if limit is not None:
                if verbose: sys.stderr.write("Filtering on best entropy\n")
                tindex = TCGA.best_entropy()[:limit]
                tcga = tcga[tindex]
            reference.append(tcga)
        if use_GTEx:
            if verbose: sys.stderr.write("loading GTEx\n")
            gtex = GTEx.get_expression()
            if limit is not None:
                if verbose: sys.stderr.write("Filtering on best entropy\n")
                gindex = GTEx.best_entropy()[:limit]
                gtex = gtex[gindex]
            reference.append(gtex)
        mutual = set()
        for r in reference:
        	mutual = mutual|set(r.index)
        # Collected all our names
        for r in reference:
        	mutual = mutual&set(r.index)
        # Limit them tot he mutual
        mutual = list(set(self._exprs.index)&mutual)
        if len(mutual) < 1: raise ValueError("Too few matching gene names with sufficient variance")
        reference = [x.loc[mutual] for x in reference]
        if verbose: sys.stderr.write("run reference correlation\n")
        rc = ReferenceCorrelation(self._exprs.loc[mutual],pd.concat(reference,1),min_variance=min_variance)
        self.features_selected = rc.features_selected
        self.reference_correlations = rc.reference_correlations
        self.data_correlations = rc.data_correlations
        self._tcga_count = 0
        self._gtex_count = 0
        if use_TCGA: self._tcga_count = tcga.columns.shape[0]
        if use_GTEx: self._gtex_count = gtex.columns.shape[0]
        self.rc = rc

    def tSNE(self,data_subset=None):
        # Get the TSNE for display
        if self._verbose: sys.stderr.write("run tSNE\n")
        if data_subset is None: data_subset = self._exprs.columns
        ecount = len(data_subset)
        df = pd.concat([self.rc.data_correlations[data_subset],self.reference_correlations],1)
        tsne = pd.DataFrame(TSNE().fit_transform(df.T))
        tsne.columns = ['x','y']
        tsne['case_submitter_id'] = df.columns
        tsne['type'] = 'observation'
        tcga_names = []
        if self._tcga_count > 0: tcga_names = tsne.iloc[ecount:ecount+self._tcga_count]['case_submitter_id']
        gtex_names = []
        if self._gtex_count > 0: gtex_names = tsne.iloc[-1*self._gtex_count:]['case_submitter_id']
        tsne.loc[tsne['case_submitter_id'].isin(tcga_names),'type'] = 'TCGA'
        tsne.loc[tsne['case_submitter_id'].isin(gtex_names),'type'] = 'GTEx'
        # label the tSNE
        tsne['label'] = np.nan
        tdata = tsne.merge(TCGA.get_meta()[['project']],left_on='case_submitter_id',right_index=True)[['case_submitter_id','project']]
        shared = list(set(tdata['case_submitter_id'])&set(tsne['case_submitter_id']))
        tsne = tsne.set_index('case_submitter_id')
        tdata = tdata.set_index('case_submitter_id')
        tsne.loc[shared,'label'] = tdata.loc[shared,'project']
        # label the GTEx
        tsne = tsne.reset_index()
        gdata = tsne.merge(GTEx.get_meta()[['SMTS']],left_on='case_submitter_id',right_index=True).reset_index()[['case_submitter_id','SMTS']]
        shared = list(set(gdata['case_submitter_id'])&set(tsne['case_submitter_id']))
        tsne = tsne.set_index('case_submitter_id')
        gdata = gdata.set_index('case_submitter_id')
        tsne.loc[shared,'label'] = gdata.loc[shared,'SMTS']
        return tsne

