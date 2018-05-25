import pandas as pd
import numpy as np
import os, math

_udir = os.path.dirname(os.path.realpath(__file__))

def _collect_best_entropy(projects,project_column,filename,myclass):
    eo = pd.read_hdf(filename)
    if projects is None: return pd.Series(list(eo.sort_values('redundancy',ascending=False)['case_submitter_id']))
    v = eo.merge(myclass.get_meta(),left_on='case_submitter_id',right_index=True).\
                 sort_values('redundancy',ascending=False)
    return pd.Series(list(v.loc[v[project_column].isin(projects)].sort_values('redundancy',ascending=False)['case_submitter_id']))

class TCGA:
    def __init__(self,*args,**kw):
        return
    @classmethod
    def get_expression(TCGA,projects=None,transform='plus1log2'):
        chosen_projects = TCGA.__select_projects(projects)
        data = pd.concat([pd.read_hdf(os.path.join(_udir,'data',x)) for x in chosen_projects['path']],1)
        if transform == 'plus1log2':
            scale = math.log(1.05)/math.log(2)
            data = data.multiply(scale)
        elif transform == 'raw':
            data = data.rpow(1.05).subtract(1)
        return data[sorted(data.columns)]
    @classmethod
    def get_meta(TCGA,projects=None): 
        chosen_projects = TCGA.__select_projects(projects)
        meta = pd.read_hdf(os.path.join(_udir,'data/TCGA/meta.h5'))
        meta = meta[meta['project'].isin(chosen_projects.index)]
        return meta.loc[sorted(meta.index)]
    @staticmethod
    def __select_projects(projects):
        chosen_projects = pd.read_hdf(os.path.join(_udir,'data/TCGA/locations.h5'))
        if projects is not None:
            missing = sorted(list(set(projects)-set(chosen_projects.index)))
            if len(missing) > 0: raise ValueError("Could not find projects named "+str(missing))
            chosen_projects = chosen_projects.loc[projects]
        return chosen_projects
    @classmethod
    def get_genes(TCGA): 
        return pd.read_hdf(os.path.join(_udir,'data/TCGA/genes.h5'))
    @classmethod
    def get_projects(TCGA): 
        return TCGA.__select_projects(None).index
    @classmethod
    def best_entropy(TCGA,projects=None):
        return _collect_best_entropy(projects,'project',os.path.join(_udir,'data/TCGA/entropy-order.h5'),TCGA)
        
class GTEx:
    def __init__(self,*args,**kw):
        return
    @classmethod
    def get_expression(GTEx,tissue_types=None,transform='plus1log2'):
        chosen_tissues = GTEx.__select_tissues(tissue_types)
        data = pd.concat([pd.read_hdf(os.path.join(_udir,'data',x)) for x in chosen_tissues['path']],1)
        if transform == 'plus1log2':
            scale = math.log(1.05)/math.log(2)
            data = data.multiply(scale)
        elif transform == 'raw':
            data = data.rpow(1.05).subtract(1)
        return data[sorted(data.columns)]
    @classmethod
    def get_meta(GTEx,tissue_types=None):
        chosen_tissues = GTEx.__select_tissues(tissue_types)
        meta = pd.read_hdf(os.path.join(_udir,'data/GTEx/meta.h5'))
        meta = meta[meta['SMTS'].isin(chosen_tissues.index)]
        return meta.loc[sorted(meta.index)]
    @staticmethod
    def __select_tissues(tissue_types):
        chosen_tissues = pd.read_hdf(os.path.join(_udir,'data/GTEx/locations.h5'))
        if tissue_types is not None:
            missing = sorted(list(set(tissue_types)-set(chosen_tissues.index)))
            if len(missing) > 0: raise ValueError("Could not find tissue named "+str(missing))
            chosen_tissues = chosen_tissues.loc[tissue_types]
        return chosen_tissues
    @classmethod
    def get_genes(GTEx): 
        return pd.read_hdf(os.path.join(_udir,'data/GTEx/genes.h5'))
    @classmethod
    def get_tissue_types(GTEx): 
        return GTEx.__select_tissues(None).index
    @classmethod
    def best_entropy(GTEx,tissue_types=None):
        return _collect_best_entropy(tissue_types,'SMTS',os.path.join(_udir,'data/GTEx/entropy-order.h5'),GTEx)
