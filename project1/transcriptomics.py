import numpy as np
import csv



genes = []
data = []

with open("data/GSE150910_gene-level_count_file.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        
        genes.append(row[0])
        data.append(row[1:])


# first row is sample labels
samples = data[0]
data = data[1:]
genes = genes[1:]


# everything is read in as strings, so we need to cast the data to floats
data = np.array(data).astype(float)

# also cast samples and genes to np array
genes = np.array(genes)
samples = np.array(samples)



# shape checks
print(genes.shape)
print(samples.shape)
print(data.shape)


labels = []

for i in range(len(samples)):
    tmp = samples[i].split("_")
    labels.append(tmp[0])
    

labels = np.array(labels)


# set(labels)

data = data[:,labels != "chp"]
samples = samples[labels != "chp"]
labels = labels[labels != "chp"]


print(data.shape)
print(samples.shape)
print(labels.shape)




# rough normalization


for j in range(data.shape[1]):
    column_sum = sum(data[:,j])
    
    data[:,j] = data[:,j] / column_sum







# differential expression analysis

from scipy.stats import ks_2samp


p_values = []
log2_FCs = []


for i in range(data.shape[0]):
    
    
    control = data[i,labels == "control"]
    ipf = data[i,labels == "ipf"]
    
    # get p-value
    ks_statistic, p_value = ks_2samp(control, ipf)
    p_values.append(p_value)
    
    control_mean = np.mean(control)
    if control_mean == 0:
        control_mean = 0.00000001
        # to avoid division by 0 error

        
    ipf_mean = np.mean(ipf)
    if ipf_mean == 0:
        ipf_mean = 0.00000001
        # to avoid division by 0 error
        
        
    # get log2_FC
    FC = ipf_mean / control_mean
    log2_FC = np.log2(FC)
    log2_FCs.append(log2_FC)







print(len(p_values))

p_values = np.array(p_values)
log2_FCs = np.array(log2_FCs)


p_values_bonf = p_values * len(genes)


# print(log2_FCs)


# a = np.array([True, True, False])
# b = np.array([False, True, False])

# print(a & b)


sum((p_values_bonf <= 0.01) & (np.abs(log2_FCs) > 2))


to_keep = (p_values_bonf <= 0.05) & (np.abs(log2_FCs) >= 2)



# print(log2_FCs[to_keep])

print(sum(to_keep))






# saving significant DEGs


print(len(to_keep))
print(len(genes))
print(len(p_values_bonf))
print(len(log2_FCs))
print(data.shape[0])




sig_genes = genes[to_keep]
sig_log2_FCs = log2_FCs[to_keep]
sig_p_values_bonf = p_values_bonf[to_keep]
sig_data = data[to_keep, :]


print(len(sig_genes))

print(sig_data.shape)




# print out DEG lists

for i in range(len(sig_genes)):
    if sig_log2_FCs[i] > 0:
        print(sig_genes[i])





sig_data_log2fc = sig_data.copy()

for i in range(len(sig_genes)):
    
    control_mean =  np.mean(sig_data[i,labels == "control"])
    
    sig_data_log2fc[i,:] = np.log2((sig_data_log2fc[i,:]+0.000001) / (control_mean +0.000001))
    
    
    
import seaborn as sns
import matplotlib.pyplot as plt


sns_plot = sns.clustermap(sig_data_log2fc, xticklabels=samples, yticklabels= sig_genes, cmap="coolwarm", vmin=-2, vmax=2)

# TODO - can this be done with easier syntax?
sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize=2)
sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize=2)


sns_plot.savefig("plots/heatmap.pdf")

# plt.show()




print("testing")