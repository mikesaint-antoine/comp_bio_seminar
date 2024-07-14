import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

snp_info = []
data = []
with open("data/example_data.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=',')

    for row in reader:

        snp_info.append(row[0:3])
        data.append(row[3:])

snp_header = snp_info[0]
snp_info = snp_info[1:]

samples = data[0]
data = data[1:]


snp_info = np.array(snp_info)
data = np.array(data).astype(int)

snp_header = np.array(snp_header)
samples = np.array(samples)



labels = []

for i in range(len(samples)):
    tmp = samples[i].split("_")
    labels.append(tmp[0])
    
labels = np.array(labels)




maf_threshold = 0.05
total_alleles_per_snp = len(samples)*2


to_keep = []

for i in range(data.shape[0]):
    
    ref_allele_count = 2 * np.sum(data[i,:] == 0) + np.sum(data[i,:] == 1)
    var_allele_count = 2 * np.sum(data[i,:] == 2) + np.sum(data[i,:] == 1)

    maf = min(ref_allele_count, var_allele_count) / total_alleles_per_snp

    if maf >= maf_threshold:
        to_keep.append(True)
    else:
        to_keep.append(False)


to_keep = np.array(to_keep)

filtered_snp_info = snp_info[to_keep,:]
filtered_data = data[to_keep,:]


encoded_labels = np.where(labels == "zombie", 1, 0)


snp_genotypes = filtered_data[0, :]


p_values = []

for i in range(filtered_data.shape[0]):
  
    snp_genotypes = filtered_data[i, :]

    contingency_table = np.zeros((2, 3))

    for j in range(len(snp_genotypes)):

        genotype = snp_genotypes[j]
        status = encoded_labels[j]

        contingency_table[status, genotype] += 1

    try:
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
        
    except:
        
        # adding 1 to every element to deal with zero-column problem
        contingency_table = contingency_table + 1
        chi2, p_value, _, _ = chi2_contingency(contingency_table)

    p_values.append(p_value)

p_values = np.array(p_values)

p_values_bonf = p_values * filtered_data.shape[0]



filtered_snp_chromosomes = []
filtered_snp_positions = []


for i in range(filtered_snp_info.shape[0]):
    
    snp_id = filtered_snp_info[i,0]
    
    chromosome, position = snp_id.split(":")
    
    chromosome = chromosome.replace('chr', '')
    
    chromosome = int(chromosome)
    position = int(position)
    
    
    filtered_snp_chromosomes.append(chromosome)
    filtered_snp_positions.append(position)
    
    
filtered_snp_chromosomes = np.array(filtered_snp_chromosomes)
filtered_snp_positions = np.array(filtered_snp_positions)


neg_log_p_values = -np.log10(p_values)



snp_data = np.column_stack((filtered_snp_chromosomes, filtered_snp_positions, neg_log_p_values))

snp_data = snp_data[np.lexsort((snp_data[:, 1], snp_data[:, 0]))]

cumulative_positions = np.zeros_like(filtered_snp_positions)
chromosome_offsets = {}
midpoints = []

current_offset = 0
for chromosome in np.unique(filtered_snp_chromosomes):
    chrom_mask = filtered_snp_chromosomes == chromosome
    chromosome_offsets[chromosome] = current_offset
    cumulative_positions[chrom_mask] = filtered_snp_positions[chrom_mask] + current_offset
    midpoints.append((current_offset + np.max(cumulative_positions[chrom_mask])) / 2)
    current_offset += np.max(filtered_snp_positions[chrom_mask]) + 1  # Add a buffer between chromosomes

plt.figure(figsize=(12, 6))

colors = ['#1f77b4', '#ff7f0e']

dot_size = 10  # Adjust this value to change the size of the dots
current_chromosome = 1
for chromosome in np.unique(snp_data[:, 0]):
    chrom_mask = snp_data[:, 0] == chromosome
    plt.scatter(cumulative_positions[chrom_mask], snp_data[chrom_mask, 2], c=colors[current_chromosome % 2], s=dot_size, label=f'Chromosome {chromosome}')
    current_chromosome += 1

plt.xticks(midpoints, np.unique(filtered_snp_chromosomes))

plt.xlabel('Chromosome')
plt.ylabel('-log10(p-value)')
plt.title('Manhattan Plot of GWAS Results')

plt.show()