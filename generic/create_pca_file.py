"""
Script to create a pca file given an expression file.
Example:
    > python create_pca_file.py 6121.raw.tsv 'count' # if the file contains raw count data
    > python create_pca_file.py 6121.raw.tsv 'array' # if the file contains microarray data

This will create files named 6121.pca.tsv and 6121.pca_attributes.tsv in the current directory.
The expression files must be in tab-separated format, readable by pandas.read_csv(filename, sep='\t', index_col=0).
The columns of expression file will be sample ids, and PCA will be on these sample ids.

Packages required for this script to run: pandas, sklearn:
    > conda install pandas, scikit-learn    
"""
import os, sys, pandas, numpy
from sklearn.decomposition import PCA

def createPCAFile():
    filename = sys.argv[1]  # eg. 6121.raw.tsv
    df = pandas.read_csv(sys.argv[1], sep="\t", index_col=0)
    datatype = sys.argv[2]  # if data is count data sys.argv[2] must be specified as 'count'. If microarray, must be specified as 'array'
    if datatype == 'count':
        df = (df * 1000000).div(df.sum(axis = 0), axis = 1)    # counts per million
        df = numpy.log2(df+1)
    elif datatype == 'array':
        if df.max().max()>100: # log this first
            df = numpy.log2(df)         
    # Auto select number of pca components using this. 
    # From docs (https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html):
    #   If 0 < n_components < 1 and svd_solver == 'full', select the number of components such that the
    #   amount of variance that needs to be explained is greater than the percentage specified by n_components.
    pca = PCA(n_components=0.99, svd_solver='full')
    coords = pca.fit_transform(df.values.T)
    pandas.DataFrame(coords, index=df.columns).to_csv(filename.replace("raw","pca"), sep="\t")
    pandas.DataFrame({"explained_variance_":pca.explained_variance_,
                      "explained_variance_ratio_":pca.explained_variance_ratio_,
                      "singular_values_":pca.singular_values_}).to_csv(filename.replace("raw", "pca_attributes"), sep="\t")

if __name__=="__main__":
    createPCAFile()
