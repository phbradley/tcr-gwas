from scipy.stats import linregress, mannwhitneyu, pearsonr

def linear_reg_phil(x, y):
    regression = linregress(x, y)
    pvalue = regression.pvalue
    return(pvalue)
