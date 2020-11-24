import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.cluster import Birch
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn import metrics
from config import config
import time

def clustering(n_clusters,read_path,label_path):

    n_clusters = n_clusters
    raw_data=pd.read_csv(read_path,header=None)
    y_true = pd.read_csv(label_path,header=None)
    raw_data = pd.DataFrame(raw_data.T,index=raw_data.columns,columns=raw_data.index)
    raw_data = raw_data.fillna(0)
    y_true = pd.DataFrame(y_true.T,index=y_true.columns,columns=y_true.index)
    y_true = np.array(y_true)
    y_true = y_true.flatten()
    iteration = 5
    
    model_kmeans = KMeans(n_clusters=n_clusters,random_state=100,max_iter=iteration)

    Agglomerativeclustering = AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
                        connectivity=None, linkage='ward', memory=None, n_clusters=n_clusters)
    
    y_pre = model_kmeans.fit_predict(raw_data)
    y_Agg_pre = Agglomerativeclustering.fit_predict(raw_data)
    r1 = pd.Series(model_kmeans.labels_).value_counts()
    r2 = pd.DataFrame(model_kmeans.cluster_centers_)
    r = pd.concat([r2, r1], axis = 1)
    r.columns = list(raw_data.columns) + [u'Number of class']
    r = pd.concat([raw_data, pd.Series(model_kmeans.labels_, index = raw_data.index)], axis = 1)
    r.columns = list(raw_data.columns) + [u'Cluster category']

    n_samples,n_features=raw_data.shape
    KMeans_ARI = metrics.adjusted_rand_score(y_true,y_pre)
    KMeans_NMI = metrics.adjusted_mutual_info_score(y_true,y_pre,average_method='arithmetic')
    KMeans_HOM = metrics.homogeneity_score(y_true,y_pre)
    KMeans_COM = metrics.completeness_score(y_true,y_pre)

    Agg_ARI = metrics.adjusted_rand_score(y_true,y_Agg_pre)
    Agg_NMI = metrics.adjusted_mutual_info_score(y_true,y_Agg_pre,average_method='arithmetic')
    Agg_HOM = metrics.homogeneity_score(y_true,y_Agg_pre)
    Agg_COM = metrics.completeness_score(y_true,y_Agg_pre)


#     print("KMeans_ARI=%.4f" % KMeans_ARI)
#     print("KMeans_NMI=%.4f" % KMeans_NMI)
#     print("KMeans_HOM=%.4f" % KMeans_HOM)
#     print("KMeans_COM=%.4f" % KMeans_COM)

    print("Agg_ARI=%.4f" % Agg_ARI)
    print("Agg_NMI=%.4f" % Agg_NMI)
    print("Agg_HOM=%.4f" % Agg_HOM)
    print("Agg_COM=%.4f" % Agg_COM)
    print('+'*10)


