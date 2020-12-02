import numpy as np
import seaborn as sns
from datetime import datetime
import pickle
from scipy import sparse


import sys
sys.path.append('/tmp/FIt-SNE-master/')
from fast_tsne import fast_tsne
from sklearn.decomposition import PCA

#GLMPCA python
sys.path.append('../../libs/glmpca-py/')
from glmpca import glmpca



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

def add_largedot_legend(ax,loc):
    lgnd = ax.legend(loc=loc,frameon=True)
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

def kobak_tsne(data,name='',n_PCs=50,init=None,do_pca=True,seed=42):    
    
    if do_pca:
        X,PCAinit  = PCA_sklearn(data,n_PCs=n_PCs,seed=seed)
        print('pca is done')
    else:
        X=data
    
    if init is None:
        init = PCAinit
    
    perplexities = [30, int(X.shape[0]/100)]
    learning_rate=X.shape[0]/12

    tsne = dict(perplexity_list=perplexities, initialization=init, learning_rate=learning_rate, seed=seed)
    tsne['coords'] = fast_tsne(X, perplexity_list=perplexities, initialization=init, learning_rate=learning_rate,seed=seed)
    
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
    
def sqrt_full(counts):
    '''The same as sqrt_lazy(), but add sqrt(normalized+1)'''
    normalized = normalize_and_scale(counts,scale_mode='median')
    return np.sqrt(normalized) + np.sqrt(normalized+1)

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

def pearson_residuals(counts, theta):
    '''Computes analytical residuals for NB model with a fixed theta, clipping outlier residuals to sqrt(N)'''
    counts_sum0 = np.sum(counts, axis=0, keepdims=True)
    counts_sum1 = np.sum(counts, axis=1, keepdims=True)
    counts_sum  = np.sum(counts)

    #get residuals
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    n = counts.shape[0]
    z[z >  np.sqrt(n)] =  np.sqrt(n)
    z[z < -np.sqrt(n)] = -np.sqrt(n)
    
    return z

def plot_correlations(ax,correlations,labels,title,vmin=0.5,vmax=1):

    data = np.tril(correlations,k=-1)
    data[data==0]=np.nan

    min_val, max_val = 0, 4
    ind_array = np.arange(min_val, max_val , 1.0)
    x, y = np.meshgrid(ind_array, ind_array)

    ax.imshow(data,cmap='Reds',vmin=vmin,vmax=vmax)

    for val, x_val, y_val in zip(data.flatten(),x.flatten(), y.flatten()):
        c = '%.2f' % (val) if val>0 else ''
        ax.text(x_val, y_val, c, va='center', ha='center')

    ax.set_xlim(min_val-0.5, max_val)
    ax.set_ylim(min_val-0.5, max_val)
    ax.set_xticks(np.arange(max_val))
    ax.set_xticklabels(labels,rotation=25)
    ax.set_yticks(np.arange(max_val))
    ax.set_yticklabels(labels)

    ax.set_title(title)
    sns.despine()

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
    
    
    
    
    
    