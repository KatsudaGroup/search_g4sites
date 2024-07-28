import pandas as pd
import sys

if len(sys.argv) > 2:
    raise

#filename = sys.argv[1]
filename = 'result.tsv'
nm_df = pd.read_csv(filename,  header = 0, sep = '\t')

g22_series = nm_df['G2'] >= 2
g31_series = (nm_df['G3'] >= 1) & (nm_df['G1'] >= 1)
both_series = g22_series & g31_series
g2nd_series = nm_df['G1'] >= 2
len_g22 = len(nm_df[g22_series] )
len_g31 = len(nm_df[g31_series] )
len_both = len(nm_df[both_series])
len_g2nd = len(nm_df[g2nd_series])

g22_within_200_series = (nm_df['G2_min_interval'].notna()   ) & ( nm_df["G2_min_interval"] <= 200 )
g31_within_200_series = (nm_df['G3_G1_min_interval'].notna()) & ( nm_df["G3_G1_min_interval"] <= 200 )

print("Entry starting with NM_: {}".format( len(nm_df) ))
print("Entry of NM_, and G2 >= 2: {}".format(len(nm_df[g22_series])) )
print("Entry of NM_, and G3 >= 1 and G1 >= 1: {}".format(len(nm_df[g31_series])) )
print("Entry of NM_, and have Both motifs: {}".format(len(nm_df[g22_series & g31_series])) )
print("Entry of NM_, and have either motifs: {}".format(len(nm_df[g22_series | g31_series])) )
print("Entry of NM_, and have no motifs: {}".format(len(nm_df) - len(nm_df[g22_series | g31_series])) )
print("Entry of NM_',and have two G2 within 200: {}".format(len(nm_df[g22_within_200_series]) ))
print("Entry of NM_',and have two G3 G1 within 200: {}".format(len(nm_df[g31_within_200_series]) ))
print("Entry of NM_, and G1 >= 2: {}".format(len_g2nd) )
#print(nm_df[g22_within_200_series])
#print(nm_df[g22_within_200_series == False])

n_g22_series = len(nm_df[g22_series])
n_g31_series = len(nm_df[g31_series])
n_both_series = len(nm_df[g22_series & g31_series]) 
n_either_series = len(nm_df[g22_series | g31_series])
