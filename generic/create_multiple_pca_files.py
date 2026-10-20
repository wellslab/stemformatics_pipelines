"""
Script to create multiple pca files given a directory with many expression files.
Example:
    > python create_multiple_pca_files.py pca_files_todo > create_multiple_pca_files.log

This will create files named like 6121.pca.tsv in the specified directory.
The expression files must be in tab-separated format, readable by pandas.read_csv(filename, sep='\t', index_col=0).
The columns of expression file will be sample ids, and PCA will be on these sample ids.
Note the use of log file - if any file fails, it will write to log and continue with others.

Packages required for this script to run: pandas, sklearn:
    > conda install pandas, scikit-learn    
"""
import os, sys, pandas, numpy
from sklearn.decomposition import PCA

def createPCAFiles():
    dirname = sys.argv[1]  # directory name 
    os.chdir(dirname)
    for filename in os.listdir(dirname):
        if 'raw' in filename:
            try:
                df = pandas.read_csv(filename, sep="\t", index_col=0)
                if df.max().max()>100: # log this first
                    df = numpy.log2(df+1)

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
                        
            except:
                print("Failed for dataset", filename)

if __name__=="__main__":
    createPCAFiles()
