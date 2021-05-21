import numpy as np
import seaborn as sns
from datetime import datetime
import pickle
from scipy import sparse
from sklearn.decomposition import PCA
from scipy.special import xlogy



import sys
sys.path.append('/tmp/FIt-SNE-master/')       ### ADAPT PATH TO FIt-SNE REPO AS NEEDED
from fast_tsne import fast_tsne
sys.path.append('../../libs/glmpca-py/')      ### ADAPT PATH TO glmpca-py REPO AS NEEDED
from glmpca import glmpca
sys.path.append('../../libs/rna-seq-tsne/')   ### ADAPT PATH TO rna-seq-tsne REPO AS NEEDED
import rnaseqTools




def add_labels(dataset, xdata,ydata,example_genes,textoffsets,lines,ax):
        
    for example_gene,textoffset,line in zip(example_genes,textoffsets,lines):
        gene_idx = dataset['genes'] == example_gene
        gene_position = np.array([xdata[gene_idx],ydata[gene_idx]]).flatten()
        text_position = np.array([gene_position[0]*10**textoffset[0],gene_position[1]+textoffset[1]]).flatten()
        ax.text(*text_position,example_gene.capitalize(),horizontalalignment='center')
        if line:            
            line_start = line[0]
            line_end =  line[1]
            line_x = [text_position[0]*10**line_start[0],gene_position[0]*10**line_end[0]]
            line_y = [text_position[1]+line_start[1],gene_position[1]+line_end[1]]
            ax.plot(line_x,line_y,'k',linewidth=2)

def add_largedot_legend(ax,loc,kwargs={}):
    lgnd = ax.legend(loc=loc,frameon=True,**kwargs)
    for l in lgnd.legendHandles:
        l._sizes = [30]
        
def compute_marginals(counts):
    '''compute depths per cell (ns) and relative expression fractions per gene (ps)'''
    ns = np.sum(counts,axis=1)
    ps = np.sum(counts,axis=0)
    ps = ps / np.sum(ps)    
    return np.squeeze(np.array(ns)), np.squeeze(np.array(ps))

def compute_means(counts):
    '''compute gene means and their min, max and 50 logspaced values covering this range'''
    means = np.squeeze(np.array(np.mean(counts,axis=0)))
    mean_min = min(means)
    mean_max = max(means)
    mean_range = np.logspace(np.log10(mean_min),np.log10(mean_max))    
    return means, mean_min, mean_max, mean_range

def get_glmpca_timestamps(date_string):    
    date_object = datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S.%f")
    return date_object

def kobak_tsne(data,name='',n_PCs=50,init=None,do_pca=True,seed=42,perplexities=None):    
    
    if do_pca:
        X,PCAinit  = PCA_sklearn(data,n_PCs=n_PCs,seed=seed)
        print('pca is done')
    else:
        X=data
    
    if init is None and do_pca:
        init = PCAinit    
    
    n = X.shape[0]
    
    if perplexities is None:
        perplexities = [30, int(n/100)]
    print(perplexities)

    learning_rate= n/12

    tsne = dict(perplexity_list=perplexities, initialization=init, learning_rate=learning_rate, seed=seed)
    print('tSNE input shape:',X.shape)
    tsne['coords'] = fast_tsne(X, perplexity_list=perplexities, initialization=init, learning_rate=learning_rate,seed=seed)
    
    filename = 'tsne/tsne_%s.pickle' % (name)
    with open(filename, "wb") as file:
        pickle.dump(tsne, file, protocol=4)
    
    
    print('tSNE output shape:',tsne['coords'].shape)
    
    return tsne

def kobak_tsne_w_exag(data,name='',n_PCs=50,init=None,do_pca=True,exag=4,start_late_exag_iter=250,seed=42,print_info=False):    
    
    pass
    if do_pca:
        X,PCAinit  = PCA_sklearn(data,n_PCs=n_PCs,seed=seed)
        print('pca is done')
    else:
        X=data
    
    if init is None:
        init = PCAinit
    
    perplexities = [30]
    learning_rate=X.shape[0]/12

    tsne = dict(perplexity_list=perplexities, initialization=init, learning_rate=learning_rate, seed=seed, exag=exag,start_late_exag_iter=start_late_exag_iter)
    
    if print_info:
        print('tSNE of X with',X.shape)
        print('init with ',init.shape)
        print('perplexity_list',perplexities)
        print('late_exag_coeff',exag)
        print('start_late_exag_iter',start_late_exag_iter)
        print('learning_rate',learning_rate)
        print('seed',seed)
            
    tsne['coords'] = fast_tsne(X, perplexity_list=perplexities, late_exag_coeff=exag, start_late_exag_iter=start_late_exag_iter,initialization=init, learning_rate=learning_rate,seed=seed)
    
    filename = 'tsne/tsne_%s.pickle' % (name)
    with open(filename, "wb") as file:
        pickle.dump(tsne, file, protocol=4)
        
    return tsne

def monitor_progress(gene_id,n_genes,print_time_every=1000,print_dot_every=25):
    if np.mod(gene_id,print_time_every)==0:
        print('')
        print('##',datetime.now().time(),'##',gene_id,'of',n_genes,'genes fit',end='')
    if np.mod(gene_id,print_dot_every)==0:
        print('.',end='')

def normalize_and_scale(counts,scale_mode='user_provided_scale',scale=1):
    '''Depth-normalizes and scales count matrix'''
    depths = np.squeeze(np.array(np.sum(counts,axis=1)))
    if scale_mode=='median':
        scale = np.median(depths)
    elif scale_mode=='user_provided_scale':
        #using provided scale value
        pass
    return (counts.T/depths).T * scale

def sqrt_lazy(counts):
    '''Depth-normalizes, scales by median depth and then takes the square root.'''    
    normalized = normalize_and_scale(counts,scale_mode='median')
    return np.sqrt(normalized)

def log_transform(counts,scale_mode,scale=1):
    '''Depth-normalizes, scales by median me'''
    normalized = normalize_and_scale(counts,scale_mode=scale_mode,scale=scale)
    return np.log1p(normalized)

def PCA_sklearn(data,n_PCs,seed):
    
    X = np.array(data)
    pca = PCA(n_components=n_PCs,random_state=seed)
    X = pca.fit_transform(X)
    print('variance explained by top%u PCs: %u %%' % (n_PCs, sum(pca.explained_variance_ratio_)*100))

    #scale values relative to absolute max value
    X = X / np.max(np.abs(X))
    #use top2 PCs for intialisation
    PCAinit = X[:,:2] / np.std(X[:,0]) * .0001
    
    return X, PCAinit


def PCA_sklearn_timing(data,n_PCs,seed):
    
    X = np.array(data)
    pca = PCA(n_components=n_PCs,random_state=seed)
    X = pca.fit_transform(X)

def pearson_residuals(counts, theta, clipping=True):
    '''Computes analytical residuals for NB model with a fixed theta, clipping outlier residuals to sqrt(N)'''
    counts_sum0 = np.sum(counts, axis=0, keepdims=True)
    counts_sum1 = np.sum(counts, axis=1, keepdims=True)
    counts_sum  = np.sum(counts)

    #get residuals
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    if clipping:
        n = counts.shape[0]
        z[z >  np.sqrt(n)] =  np.sqrt(n)
        z[z < -np.sqrt(n)] = -np.sqrt(n)
    
    return z


def deviance_residuals(x, theta,mu=None):
    '''Computes deviance residuals for NB model with a fixed theta'''

    if mu is None:
        counts_sum0 = np.sum(x, axis=0, keepdims=True)
        counts_sum1 = np.sum(x, axis=1, keepdims=True)
        counts_sum  = np.sum(x)
        #get residuals
        mu = counts_sum1 @ counts_sum0 / counts_sum
    
    
    
    
    def remove_negatives(sqrt_term):
        negatives_idx = sqrt_term < 0
        if np.any(negatives_idx):
            n_negatives = np.sum(negatives_idx)
            print('Setting %u negative sqrt term values to 0 (%f%%)' % (n_negatives,n_negatives/np.product(sqrt_term.shape)))
            sqrt_term[negatives_idx] = 0
    

    if np.isinf(theta): ### POISSON
        x_minus_mu = x-mu
        sqrt_term =                          2 * (xlogy(x,x/mu) - x_minus_mu   ) #xlogy(x,x/mu) computes xlog(x/mu) and returns 0 if x=0
        remove_negatives(sqrt_term)
        dev = np.sign(x_minus_mu) * np.sqrt(sqrt_term)
    else:               ### NEG BIN
        x_plus_theta = x+theta
        sqrt_term =                    2 * ( xlogy(x,x/mu)     -   (x_plus_theta)  * np.log(x_plus_theta/(mu+theta))     ) #xlogy(x,x/mu) computes xlog(x/mu) and returns 0 if x=0
        remove_negatives(sqrt_term)
        dev = np.sign(x-mu) * np.sqrt(sqrt_term)
    
    return dev


def prepare_largest_batch(dataset,extra_keys=[]):
    ###add counts and clusters for largest batch only to the dataset
    
    #find largest batch
    batch_ids, batch_counts = np.unique(dataset['batches'],return_counts=True)
    largest_batch_id_idx = np.argmax(batch_counts)
    largest_batch_id = batch_ids[largest_batch_id_idx]
    largest_batch_idx = dataset['batches'] == largest_batch_id
    
    #slice accordingly and clean out rare genes
    dataset['counts'],dataset['genes'] = remove_rare_genes(dataset['counts'][largest_batch_idx,:],dataset['genes'],5)
    dataset['clusters'] = dataset['clusters'][largest_batch_idx]
    dataset['labelshort'] = dataset['labelshort'] + '_largestBatch'
    dataset['batches'] = dataset['batches'][largest_batch_idx]
    for k in extra_keys:
        dataset[k] = dataset[k][largest_batch_idx]

def remove_rare_genes(counts,genes,minimum_detected_cells_per_gene):
    
    if type(counts) in [sparse.csr.csr_matrix, sparse.csc.csc_matrix]:
        
        #remove zero genes
        nonzero_genes_idx = np.array(counts.sum(axis=0)).flatten() > 0

        counts = counts[:,nonzero_genes_idx]
        genes = genes[nonzero_genes_idx]

        #count nonzero entries per gene
        nonzero_coords = counts.nonzero()
        n_nonzero = counts.count_nonzero()
        is_nonzero = sparse.csc_matrix((np.ones(n_nonzero),nonzero_coords))
        detected_cells_per_gene = np.array(is_nonzero.sum(axis=0)).flatten()

        keep_genes = detected_cells_per_gene >= minimum_detected_cells_per_gene    
        counts_kept = counts[:,keep_genes]
        genes_kept = genes[keep_genes]

        print('Of',len(detected_cells_per_gene),'total genes, returning',sum(keep_genes),'genes that are detected in %u or more cells.' % (minimum_detected_cells_per_gene))
        print('Output shape:', counts_kept.shape)

        return counts_kept,np.array(genes_kept)
    
    else:
        
        #remove zero genes
        nonzero_genes_idx = np.sum(counts,axis=0) > 0
        counts = counts[:,nonzero_genes_idx]
        genes = genes[nonzero_genes_idx]        
        
        #remove genes that are detected in less then n cells
        nonzero = counts > 0
        cells_per_gene = np.sum(nonzero,axis=0)
        include_genes = cells_per_gene >= minimum_detected_cells_per_gene
        counts_kept = counts[:,include_genes]
        genes_kept = genes[include_genes]
        print('Of',len(cells_per_gene),'total genes, returning',sum(include_genes),'genes that are detected in %u or more cells.' % (minimum_detected_cells_per_gene))
        print('Output shape:', counts_kept.shape)
        return counts_kept,genes_kept
    
def run_glmpca(counts,fam,theta = 100, penalty = 1, optimize_nb_theta=True, maxIter=1000, eps=0.0001, n_PCs=50, seed=42, dataset_label=''):
    '''Wrapper around GLM PCA by Will Townes: applies GLM PCA with given settings and saves results as pickle'''
    
    
    np.random.seed(seed)
    ctl = {"maxIter":maxIter, "eps":eps, "optimizeTheta":optimize_nb_theta}
    if maxIter==1000 and eps == 0.0001:
        ctl_str=''
    else:
        ctl_str='_maxIter%u_eps%s' % (maxIter,eps)

    starttime = str(datetime.now())
    res = glmpca.glmpca(counts.T,n_PCs,fam=fam,nb_theta=theta,verbose=True,penalty=penalty,ctl=ctl)
    endtime = str(datetime.now())
    res['starttime']=starttime
    res['endtime']=endtime

    if fam=='nb':
        res['nb_theta']=res['glmpca_family'].nb_theta

    _ = res.pop('glmpca_family')

    if fam=='poi':
        path = 'glmpca_results/%sglmpca-py_%s_penalty%u%s.pickle' % (dataset_label,fam,penalty,ctl_str)
    elif optimize_nb_theta:
        path = 'glmpca_results/%sglmpca-py_%s_penalty%u%s.pickle' % (dataset_label,fam,penalty,ctl_str)
    else:
        path = 'glmpca_results/%sglmpca-py_%s_fixedTheta%u_penalty%u%s.pickle' % (dataset_label,fam,theta,penalty,ctl_str)

    print('Saving at', path)
    with open(path,'wb') as f:
        pickle.dump(res,f)
    
    
### Code extended from https://github.com/theislab/scanpy/pull/1715 for Deviance residuals
    
## imports for scanpy env
import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata import AnnData

from scanpy import logging as logg

from scanpy._settings import settings, Verbosity
from scanpy._utils import sanitize_anndata, check_nonnegative_integers, view_to_actual
from scanpy.get import _get_obs_rep, _set_obs_rep
from scanpy._compat import Literal
from scanpy.preprocessing._utils import _get_mean_var
from scanpy.preprocessing._distributed import materialize_as_ndarray
from scanpy.preprocessing._simple import filter_genes
    
    
def highly_variable_residuals(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: int = 1000,
    batch_key: Optional[str] = None,
    theta: float = 100,
    clip: Optional[float] = None,
    chunksize: int = 100,
    check_values: bool = True,
    subset: bool = False,
    inplace: bool = True,
    residual_type: str = 'pearson',
    debug=False,
) -> Optional[pd.DataFrame]:
    """\
    See `highly_variable_genes`.

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pd.DataFrame`)
    or updates `.var` with the following fields:

    highly_variable
        boolean indicator of highly-variable genes.
    means
        means per gene.
    variances
        variances per gene.
    residual_variances
        Pearson residual variance per gene. Averaged in the case of multiple
        batches.
    highly_variable_rank
        Rank of the gene according to residual variance, median rank in the
        case of multiple batches.
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are
        detected as HVG.
    highly_variable_intersection : bool
        If batch_key is given, this denotes the genes that are highly variable
        in all batches.
    """

    
    view_to_actual(adata)
    X = _get_obs_rep(adata, layer=layer)
    computed_on = layer if layer else 'adata.X'
    

    # Check for raw counts
    if check_values and (check_nonnegative_integers(X) == False):
        warnings.warn(
            "`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
            UserWarning,
        )

    if batch_key is None:
        batch_info = np.zeros(adata.shape[0], dtype=int)
    else:
        batch_info = adata.obs[batch_key].values
    n_batches = len(np.unique(batch_info))

    # Get pearson residuals for each batch separately
    residual_gene_vars = []
    for batch in np.unique(batch_info):
        
        adata_subset = adata[batch_info == batch]
        

        # Filter out zero genes
        with settings.verbosity.override(Verbosity.error):
            nonzero_genes = filter_genes(adata_subset, min_cells=1, inplace=False)[0]
        adata_subset = adata_subset[:, nonzero_genes]
        

        if layer is not None:
            X_batch = adata_subset.layers[layer]
        else:
            X_batch = adata_subset.X

        # Prepare clipping
        if clip is None:
            n = X_batch.shape[0]
            clip = np.sqrt(n)
        if clip < 0:
            raise ValueError("Pearson residuals require `clip>=0` or `clip=None`.")

        if sp_sparse.issparse(X_batch):
            sums_genes = np.sum(X_batch, axis=0)
            sums_cells = np.sum(X_batch, axis=1)
            sum_total = np.sum(sums_genes).squeeze()
        else:
            sums_genes = np.sum(X_batch, axis=0, keepdims=True)
            sums_cells = np.sum(X_batch, axis=1, keepdims=True)
            sum_total = np.sum(sums_genes)

        # Compute pearson residuals in chunks
        residual_gene_var = np.empty((X_batch.shape[1]))
        for start in np.arange(0, X_batch.shape[1], chunksize):
            stop = start + chunksize

            X_dense = X_batch[:, start:stop].toarray()
            mu = np.array(sums_cells @ sums_genes[:, start:stop] / sum_total) 
            if residual_type == 'pearson':                               
                residuals = (X_dense - mu) / np.sqrt(mu + mu ** 2 / theta)
                residuals = np.clip(residuals, a_min=-clip, a_max=clip)
            elif residual_type == 'deviance':
                residuals = deviance_residuals(X_dense,theta,mu)
            residual_gene_var[start:stop] = np.var(residuals, axis=0)


        # Add 0 values for genes that were filtered out
        zero_gene_var = np.zeros(np.sum(~nonzero_genes))
        residual_gene_var = np.concatenate((residual_gene_var, zero_gene_var))
        # Order as before filtering
        idxs = np.concatenate((np.where(nonzero_genes)[0], np.where(~nonzero_genes)[0]))
        residual_gene_var = residual_gene_var[np.argsort(idxs)]
        residual_gene_vars.append(residual_gene_var.reshape(1, -1))



    residual_gene_vars = np.concatenate(residual_gene_vars, axis=0)

    # Get cutoffs and define hvgs per batch
    residual_gene_vars_sorted = np.sort(residual_gene_vars, axis=1)
    cutoffs_per_batch = residual_gene_vars_sorted[:, -n_top_genes]
    highly_variable_per_batch = np.greater_equal(
        residual_gene_vars.T, cutoffs_per_batch
    ).T

    # Merge hvgs across batches
    highly_variable_nbatches = np.sum(highly_variable_per_batch, axis=0)
    highly_variable_intersection = highly_variable_nbatches == n_batches

    # Get rank per gene within each batch
    # argsort twice gives ranks, small rank means most variable
    ranks_residual_var = np.argsort(np.argsort(-residual_gene_vars, axis=1), axis=1)
    ranks_residual_var = ranks_residual_var.astype(np.float32)
    ranks_residual_var[ranks_residual_var >= n_top_genes] = np.nan
    ranks_masked_array = np.ma.masked_invalid(ranks_residual_var)
    # Median rank across batches,
    # ignoring batches in which gene was not selected
    medianrank_residual_var = np.ma.median(ranks_masked_array, axis=0).filled(np.nan)
    
    means, variances = materialize_as_ndarray(_get_mean_var(X))

    
    df = pd.DataFrame.from_dict(
        dict(
            means=means,
            variances=variances,
            residual_variances=np.mean(residual_gene_vars, axis=0),
            highly_variable_rank=medianrank_residual_var,
            highly_variable_nbatches=highly_variable_nbatches,
            highly_variable_intersection=highly_variable_intersection,
        )
    )
    df = df.set_index(adata.var_names)

    # Sort genes by how often they selected as hvg within each batch and
    # break ties with median rank of residual variance across batches
    df.sort_values(
        ['highly_variable_nbatches', 'highly_variable_rank'],
        ascending=[False, True],
        na_position='last',
        inplace=True,
    )
    df['highly_variable'] = False
    df.highly_variable.iloc[:n_top_genes] = True
    # TODO: following line raises a pandas warning
    # (also for flavor = seurat and cellranger..)
    df = df.loc[adata.var_names]

    if inplace or subset:
        adata.uns['hvg'] = {'flavor': 'pearson_residuals', 'computed_on': computed_on}
        logg.hint(
            'added\n'
            '    \'highly_variable\', boolean vector (adata.var)\n'
            '    \'highly_variable_rank\', float vector (adata.var)\n'
            '    \'highly_variable_nbatches\', int vector (adata.var)\n'
            '    \'highly_variable_intersection\', boolean vector (adata.var)\n'
            '    \'means\', float vector (adata.var)\n'
            '    \'variances\', float vector (adata.var)\n'
            '    \'residual_variances\', float vector (adata.var)'
        )
        adata.var['highly_variable'] = df['highly_variable'].values
        adata.var['highly_variable_rank'] = df['highly_variable_rank'].values
        adata.var['means'] = df['means'].values
        adata.var['variances'] = df['variances'].values
        adata.var['residual_variances'] = df['residual_variances'].values.astype(
            'float64', copy=False
        )
        if batch_key is not None:
            adata.var['highly_variable_nbatches'] = df[
                'highly_variable_nbatches'
            ].values
            adata.var['highly_variable_intersection'] = df[
                'highly_variable_intersection'
            ].values
        if subset:
            adata._inplace_subset_var(df['highly_variable'].values)
    else:
        if batch_key is None:
            df = df.drop(
                ['highly_variable_nbatches', 'highly_variable_intersection'], axis=1
            )
        return df
    
    
    
    
