import pandas as pd
df = pd.read_csv('search_result.tsv',  header = 0, sep = '\t')

nm_df = df[ df['type'] == 'NM' ]

print("All Entry: {}".format(len(df)))
print("Entry starting with NM_: {}".format( len(nm_df) ))
print("Entry of NM_, and G2 >= 2: {}".format( len(nm_df[nm_df['G2']>=2]) ) )
print("Entry of NM_, and G3 >= 1 and G1 >= 1: {}".format(len( nm_df[ (nm_df['G3']>=1) & (nm_df['G1']>=1) ] )) )
print( nm_df[ (nm_df['G3']>=1) & (nm_df['G1']>=1) ] ) 
