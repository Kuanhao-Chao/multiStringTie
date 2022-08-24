# Import libraries
import matplotlib.pyplot as plt
import numpy as np

with open('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio.txt', 'r') as f: 
    lines = f.read().splitlines() 
    arr = np.array(lines)
    # print(arr)
    print(len(arr))
    data = arr.astype(float)
    data = data[~np.isnan(data)]
    print(data)
    print(len(data))

# Creating dataset
print(data)

print(">> Total number of node: ", )
print(len(data))

# 50%
print(">> Ratio > 1.4: ", )
print(np.sum(data > 1.4))
print(">> Ratio < 0.6: ", )
print(np.sum(data < 0.6))

fig = plt.figure(figsize =(10, 7))
 
# Creating plot
# Do not show outlier
# plt.boxplot(data, showfliers=False)
# plt.savefig('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio_no_outlier.png')

# Show outlier
plt.boxplot(data)
plt.savefig('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio_outlier.png')

# show plot
