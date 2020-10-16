import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r

MASS = importr('MASS')
stats = importr('stats')


def extract_errors(errors_raw):
    '''unpacks error structure from fitting procedure into array'''
    
    errors = []
    for error in errors_raw:
        if error is None:
            errors.append('no error')
        else:
            errors.append(str(error['error']))
            
    
    return np.array(errors)

def extract_warnings(warning_lists_raw):
    '''unpacks warning structure from fitting procedure and sorts by type'''
    #find out which types of warnings we have in the data
    warning_msgs_unstructured = []
    warning_msgs_structured = []
    for gene_id,warning_list in enumerate(warning_lists_raw):

        if len(warning_list) > 0:
            warnings_thisgene = [str(w) for w in warning_list]
            warning_msgs_unstructured.extend(warnings_thisgene)
            warning_msgs_structured.append(warnings_thisgene)
        else:
            warning_msgs_unstructured.append('"no warning"')
            warning_msgs_structured.append(['"no warning"'])


    warning_types_raw = np.unique(warning_msgs_unstructured)
    
    #now, check for each warning type for which genes it occured
    warning_type_present_idx = []
    for warning_type in warning_types_raw:
        warning_type_present_idx.append(np.array([warning_type in warnings_thisgene for warnings_thisgene in warning_msgs_structured]))

    #clean up warning names
    warning_types_clean = [wt.split('"')[1] for wt in warning_types_raw]
    
    
    return dict(warning_types_clean=np.array(warning_types_clean), warning_types_raw=np.array(warning_types_raw), warning_type_present_idx=np.array(warning_type_present_idx))

def r_define_warning_handler():
    '''function that handles R warnings:
    evalutes 'expr' and returns a list of resulting values
    and warnings that occurred during execution'''
    
    
    r("""withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, list(w))
        invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
    }""")
    
def r_extract_results():
    '''extracts R warnings and makes them python-readable'''
    
    #extract theta values and warnings
    r("theta <- res$value")
    r("warningMessages <- list()")
    r("for (i in 1:length(res$warnings)) {warningMessages <- c(warningMessages,res$warnings[[i]]$message)}")
    
    #convert to python-only datatypes
    theta = list(r("theta"))[0]
    warning_msg_R = list(r("warningMessages"))  
    
    return theta, warning_msg_R


def r_glmPoisson_MLtheta(x,y,thetaML_max_iter=10):
    '''Fits a NB regression in R as in Hafemeister&Satija 2019.
       (See Methods -> Speed considerations -> 1st bullet point)
       First, fits Poisson GLM with intercept (beta0) and slope(beta1).
       Then, ML estimate of overdispersion theta is obtained.'''
    
    robjects.numpy2ri.activate()
    r.assign("y", y)
    r.assign("x", x)
    robjects.numpy2ri.deactivate()
    r("x <- as.matrix(x)")
    r("y <- as.matrix(y)")
    r("fit <- glm(y~x,family='poisson')")
    
    r_define_warning_handler()
    
    r("res <- withWarnings(as.numeric(x = theta.ml(y = y, mu = fit$fitted, limit = %u, trace = TRUE)))" % (thetaML_max_iter))

    theta, warning_msg_R = r_extract_results()
    
    beta0, beta1 = list(r("fit$coefficients"))
    
    return beta0,beta1,theta, warning_msg_R

def r_MLtheta(y,mu,max_iter=10):
    '''Computes ML estimate for overdispersion theta in R, given predicted and observed counts.'''    

    robjects.numpy2ri.activate()
    r.assign("y", y)
    r.assign("mu", mu)
    robjects.numpy2ri.deactivate()
    
    r_define_warning_handler()    
    
    #convert mu and y to correct format and compute
    r("mu <- as.matrix(mu)")
    r("y <- as.matrix(y)")    
    r("res = withWarnings(as.numeric(x = theta.ml(y = y, mu = mu, limit = %u)))" % (max_iter))
        
    return r_extract_results()




