###
###	ybyra 3.0
###
###

import pandas as pd
import argparse

#
#	Majority filter (>70% non-clonal reads)
#

def depth_filter(row):
    x, y, w, z = map(int, row['dp4'].split(','))
    total = x + y + w + z
    return (x + y) / total > 0.7 or (w + z) / total > 0.7

#
#	Damage filter
#	Check if damage depending on lib type and if majority without potentially damaged reads
#

def assess_damage(row, lib_type):
    x, y, w, z = map(int, row['dp4'].split(','))
    
    if row['mutation'] == 'C>T':
        if lib_type in ['ss', 'both']:
            return 'yes'
        elif lib_type == 'ds':
            if y + z == 0:
                return 'yes'
            elif (y / (y + z) > 0.7 or z / (y + z) > 0.7):
                return 'no'
            else:
                return 'yes'
    
    elif row['mutation'] == 'G>A':
        if lib_type == 'ss':
            return 'no'
        elif lib_type in ['ds', 'both']:
            if x + w == 0:
                return 'yes'
            elif (x / (x + w) > 0.7 or w / (x + w) > 0.7):
                return 'no'
            else:
                return 'yes'
    
    return 'no'  # default – plus all other mutations

#
#	Fiat lux
#

def main(lib_type, out_prefix, alleles_file, snpinfo_file):
    
    snpinfo = pd.read_csv(snpinfo_file, sep='\t')
    alleles = pd.read_csv(alleles_file, sep='\t', names=['position', 'ref', 'alt', 'coverage', 'qual', 'dp4', 'geno'])

    derived_calls = pd.merge(alleles, snpinfo, left_on=['position', 'geno'], right_on=['position', 'der'])[['position', 'snpId', 'mutation', 'id', 'parent', 'coverage', 'dp4']]
    derived_calls['state'] = 'derived'
    flt_derived_calls = derived_calls[derived_calls.apply(depth_filter, axis=1)]
    nopass_derived_calls = derived_calls[~derived_calls.index.isin(flt_derived_calls.index)]
    
    ancestral_calls = pd.merge(alleles, snpinfo, left_on=['position', 'geno'], right_on=['position', 'anc'])[['position', 'snpId', 'mutation', 'id', 'parent', 'coverage', 'dp4']]
    ancestral_calls['state'] = 'ancestral'
    flt_ancestral_calls = ancestral_calls[ancestral_calls.apply(depth_filter, axis=1)]
    nopass_ancestral_calls = ancestral_calls[~ancestral_calls.index.isin(flt_ancestral_calls.index)]

    mergecolumns = ['position', 'snpId', 'mutation', 'id', 'parent', 'coverage', 'dp4', 'state']
    hits = pd.concat([flt_derived_calls[mergecolumns], flt_ancestral_calls[mergecolumns]])
    nopass = pd.concat([nopass_derived_calls[mergecolumns], nopass_ancestral_calls[mergecolumns]])

    nopass.to_csv(f'{out_prefix}.nopass', sep='\t', index=False)

	# damage filter
    hits['damage'] = hits.apply(assess_damage, axis=1, lib_type=lib_type)
    hits.to_csv(f'{out_prefix}.calls', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Match to tree and filter Y-chromosome genotypes")
    parser.add_argument('--lib', choices=['ds', 'ss', 'both'], required=True, help="Library build: 'ds' for double-stranded, 'ss' for single-stranded, 'both' if both were merged in bam")
    parser.add_argument('--out', required=True, help="Output file prefix")
    parser.add_argument('--alleles', required=True, help="Path to the alleles file")
    parser.add_argument('--snpinfo', required=True, help="Path to the SNP info file")
    args = parser.parse_args()

    main(args.lib, args.out, args.alleles, args.snpinfo)
    
    


