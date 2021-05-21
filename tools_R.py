import pickle
from scipy import sparse
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri
import scanpy as sc
import anndata
import pandas as pd
from datetime import datetime



def glmpca_R(data,fam,tol=10**-4,label='',testrun=False, batch_size=1000,seed=202, min_iter=30, max_iter=1000, learning_rate=0.1):    
  

    numpy2ri.activate()
    pandas2ri.activate()

    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1) 
    if not rpackages.isinstalled('glmpca'):
        utils.install_packages("glmpca")
    glmnet = rpackages.importr('glmpca')
    glmpca = robjects.r['glmpca']

    if testrun:
        print('WARNING: Testrun GLMPCA. using only small subset of the data. ')
        data = data[:200,:100]
        adata = anndata.AnnData(data)
        sc.pp.filter_genes(adata,min_counts=1)
        data = adata.X
        label = label + '_testrun_'
        batch_size=100

    basepath = 'glmpca_results/glmpca%s_%s_tol%s_iter%u-%u_lr%s' % (label,fam,tol,min_iter,max_iter,learning_rate)
    print('will save at %s'%basepath)

    print('using input data with shape:',data.shape)

    #create ctl list R object for controlling glmpca
    list_r = robjects.r('list')
    ctl = list_r(tol=tol,
                 batch_size=batch_size,
                 minIter=min_iter,
                 maxIter=max_iter,
                 lr=learning_rate,
                 verbose=True)

    #make use of sparsity
    Y_sparse = data.T.tocoo() ## transpose needed because glmpca wants rows as features
    r_Matrix = rpackages.importr("Matrix")
    Y_r = r_Matrix.sparseMatrix(i=robjects.IntVector(Y_sparse.row + 1),
                                j=robjects.IntVector(Y_sparse.col + 1),
                                x=robjects.IntVector(Y_sparse.data),
                                dims=robjects.IntVector(Y_sparse.shape))

    ### Run GLM PCA in R
    robjects.r("set.seed(%u)"%seed)
    starttime = str(datetime.now())
    res = glmpca(Y_r,50,fam=fam,minibatch='stochastic',ctl=ctl)    
    endtime = str(datetime.now())


    ### Convert to python friendly shape and save
    res_dict = {name:data for name,data in res.items()}
    res_dict['seed']=seed
    res_dict['ctl'] = {name:data for name,data in res_dict['ctl'].items()}
    gf = {name:data for name,data in res_dict['glmpca_family'].items()}

    res_dict_py = {}
    res_dict_py['factors'] = pd.DataFrame(res_dict['factors']).values
    res_dict_py['loadings'] = pd.DataFrame(res_dict['loadings']).values
    res_dict_py['X'] = pd.DataFrame(res_dict['X']).values
    res_dict_py['coefX'] = pd.DataFrame(res_dict['coefX']).values
    res_dict_py['dev'] = res_dict['dev']
    res_dict_py['dev_smooth'] = res_dict['dev_smooth']
    res_dict['ctl'].pop('verbose')
    res_dict_py['ctl'] = res_dict['ctl']
    res_dict_py['offsets'] = res_dict['offsets']
    res_dict_py['optimizer'] = str(res_dict['optimizer']).split('"')[1]
    res_dict_py['minibatch'] = str(res_dict['minibatch']).split('"')[1]
    res_dict_py['seed'] = res_dict['seed']
    res_dict_py['nb_theta'] = float(gf['nb_theta'])
    res_dict_py['starttime']=starttime
    res_dict_py['endtime']=endtime

    with open(basepath + '_results.pickle' ,"wb") as f:
        pickle.dump(res_dict_py,f)