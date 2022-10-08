# dat_values = [train_X[i] for i in range(len(train_X))]
# dat_values=[]
# for i in range(len(train_X)):
#     dat_tmp = train_X[i]
#     for j in range(len(dat_tmp)):
#         dat_values.extend(dat_tmp[j])
# plt.figure(figsize = (15,8))
# sns.set(font_scale = 2,style = 'whitegrid')
# ax = sns.histplot(dat_values,bins = 50,kde = True,color = 'grey')
# sns.despine()
# ax.set_title('Distribution of PSSM values',weight = 'bold')
# ax.set_ylabel('Counts')
# ax.set_xlabel('Values')
# ax.axvline(np.mean(dat_values),color = 'red',linestyle = '--')
# ax.text(np.mean(dat_values),800,'Median = '+str(np.mean(dat_values)))
# plt.savefig('/Users/suhancho/Figures/pssm.value.distribution.pdf')
# plt.close()
