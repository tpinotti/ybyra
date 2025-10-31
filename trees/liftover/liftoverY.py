##
##	liftover ybyra-like files from hg38 to hg37
##	uses CrossMap as a subprocess: https://github.com/liguowang/CrossMap
##
##	TPinotti oct/25
##


import sys
import pandas as pd
import subprocess

bed_file, tsv_file, chain_file = sys.argv[1:]

bed_df = pd.read_csv(bed_file, sep='\t', comment='#', header=None, names=['chr','start','end'])

# remove chr
bed_df['chr'] = bed_df['chr'].str.replace('^chr','',regex=True)

# make bed with snpId
tsv_df = pd.read_csv(tsv_file, sep='\t')
tsv_merge = tsv_df[['snpId','position']].copy()
tsv_merge = tsv_merge.rename(columns={'position':'end'})
merged_bed = pd.merge(bed_df, tsv_merge, on='end', how='inner')[['chr','start','end','snpId']]

# save bed to disk or CrossMap won't work
nochr_bed = f"nochr.{bed_file}"
merged_bed.to_csv(nochr_bed, sep='\t', header=False, index=False)

mapped_bed = f"hg37.{bed_file}"
unmapped_bed = f"liftoverfail.hg38.{bed_file}"

# run CrossMap
subprocess.run(["CrossMap","bed",chain_file,nochr_bed,mapped_bed,"--unmap-file",unmapped_bed],check=True)

mapped_df = pd.read_csv(mapped_bed, sep='\t', header=None, names=['chr','start','end','snpId'])
lifted_tsv_df = pd.merge(mapped_df, tsv_df, on='snpId', how='left')
lifted_tsv_df['position'] = lifted_tsv_df['start'] + 1
lifted_tsv_df.drop(columns=['start','end'], inplace=True)
lifted_tsv_df.to_csv(f"hg37.{tsv_file}", sep='\t', index=False)

unmapped_df = pd.read_csv(unmapped_bed, sep='\t', header=None, names=['chr','start','end','snpId'])
fail_tsv_df = pd.merge(unmapped_df, tsv_df, on='snpId', how='left')
fail_tsv_df.to_csv(f"liftoverfail.{tsv_file}", sep='\t', index=False)

print(f"{len(fail_tsv_df)} positions failed liftover hg38 > hg37")

