# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import sys

with open(sys.argv[1], 'r') as f: 
    lines = f.read().splitlines() 
    lines_split = [[float(line.split()[0]), float(line.split()[1])] for line in lines]

    arr = np.array(lines_split)

filter_arr = arr[:,0] < 10000
newarr = arr[filter_arr]
# newarr = arr

#     print(len(arr))
#     data = arr.astype(float)
#     data = data[~np.isnan(data)]
#     print(data)
#     print(len(data))

# # Creating dataset
# print(data)

# print(">> Total number of node: ", )
# print(len(data))

# # 50%
# print(">> Ratio > 1.4: ", )
# print(np.sum(data > 1.4))
# print(">> Ratio < 0.6: ", )
# print(np.sum(data < 0.6))

fig = plt.figure(figsize =(10, 10))
 
# Creating plot
# Do not show outlier
# plt.boxplot(data, showfliers=False)
# plt.savefig('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio_no_outlier.png')


plt.title("StringTie vs multiStringTie Negative coverage (before dividing the read length)")
plt.xlabel('StringTie coverage')
plt.ylabel('multiStringTie coverage')

plt.scatter(newarr[:, 0], newarr[:, 1])
plt.show()
# plt.savefig('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/coverage_sub_chr22.png')
# plt.savefig('/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/coverage_all_chr22_filter.png')

# show plot
