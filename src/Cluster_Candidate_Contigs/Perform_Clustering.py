import numpy as np
from os import mkdir
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

op_dir = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Clustering_Results_Distance/'
try:
	mkdir(op_dir)
except Exception:
	pass

distance_levels = np.arange(0.4,1.0,0.01)
filepath = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Distance_Matrix_K_11.npy'
Distance_Matrix = np.load(filepath)

Distance_Matrix = np.maximum(Distance_Matrix, Distance_Matrix.transpose())
cond_vec = squareform(Distance_Matrix)

single_linkage = linkage(cond_vec, method='single')
complete_linkage = linkage(cond_vec, method='complete')
avg_linkage = linkage(cond_vec, method='average')

ward_linkage = linkage(cond_vec, method = 'ward')
weighted_linkage = linkage(cond_vec, method = 'weighted')
centroid_linkage = linkage(cond_vec, method = 'centroid')

for t in distance_levels:
	print (t)
	cluster_id_sing = fcluster(single_linkage, criterion = 'distance',  t = t)
	cluster_id_comp = fcluster(complete_linkage, criterion = 'distance',  t = t)
	cluster_id_avg = fcluster(avg_linkage, criterion = 'distance',  t = t)
	
	cluster_id_ward = fcluster(ward_linkage, criterion = 'distance', t = t)	
	cluster_id_centroid = fcluster(centroid_linkage, criterion = 'distance', t = t)
	cluster_id_weighted = fcluster(weighted_linkage, criterion = 'distance', t = t)

	np.save(op_dir+'Distance_Single_'+str(round(t, 4)), cluster_id_sing)
	np.save(op_dir+'Distance_Complete_'+str(round(t, 4)), cluster_id_comp)
	np.save(op_dir+'Distance_Average_'+str(round(t, 4)), cluster_id_avg)
	
	np.save(op_dir+'Distance_Ward_'+str(t), cluster_id_ward)
	np.save(op_dir+'Distance_Weighted_'+str(t), cluster_id_weighted)
	np.save(op_dir+'Distance_Centroid_'+str(t), cluster_id_centroid)	
