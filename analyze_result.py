import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
#filename = 'search_result.tsv'
df = pd.read_csv(filename,  header = 0, sep = '\t')

nm_df = df[ df['type'] == 'NM' ]

g22_series = nm_df['G2'] >= 2
g31_series = (nm_df['G3'] >= 1) & (nm_df['G1'] >= 1)
both_series = g22_series & g31_series
g2nd_series = nm_df['G1'] >= 2
len_g22 = len(nm_df[g22_series] )
len_g31 = len(nm_df[g31_series] )
len_both = len(nm_df[both_series])
len_g2nd = len(nm_df[g2nd_series])

print("All Entry: {}".format(len(df)))
print("Entry starting with NM_: {}".format( len(nm_df) ))
print("Entry of NM_, and G2 >= 2: {}".format(len(nm_df[g22_series])) )
print("Entry of NM_, and G3 >= 1 and G1 >= 1: {}".format(len(nm_df[g31_series])) )
print("Entry of NM_, and have Both motifs: {}".format(len(nm_df[g22_series & g31_series])) )
print("Entry of NM_, and have either motifs: {}".format(len(nm_df[g22_series | g31_series])) )
print("Entry of NM_, and have no motifs: {}".format(len(nm_df) - len(nm_df[g22_series | g31_series])) )
print("Entry of NM_, and G1 >= 2: {}".format(len_g2nd) )

plt.title('Human mRNAs (NM_*): {}'.format(len(nm_df)))
n_g22_series = len(nm_df[g22_series])
n_g31_series = len(nm_df[g31_series])
n_both_series = len(nm_df[g22_series & g31_series]) 
n_either_series = len(nm_df[g22_series | g31_series])
venn2(
    subsets = ([n_g22_series - n_both_series, n_g31_series - n_both_series, n_both_series] ), 
    set_labels = ('2+2', '3+1')  )
plt.show()
