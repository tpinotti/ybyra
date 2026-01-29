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
#	Check if damage depending on lib type and damage model
#

def assess_damage(row, lib_type, dmg_model):
    x, y, w, z = map(int, row['dp4'].split(','))
    mut = row['mutation']
    state = row['state']

    if lib_type in ['ds', 'both']:
        right = ['C>T', 'G>A']
        left = ['T>C', 'A>G']
    else:	# 'ss'
        right = ['C>T']
        left = ['T>C']

    #
    #	naive mode – simply remove all C>T (or G>A)
    #

    if dmg_model == 'naive':
        if mut in right:
            return 'yes'
        return 'no'

    #
    #	unidirectional mode – remove all derived C>T (or G>A); keep all ancestral
    #

    elif dmg_model == 'uni':
        if mut in right and state == 'derived':
            return 'yes'
        return 'no'

    #
    #	bidirectional mode – remove all derived C>T (or G>A); remove all ancestral T>C (or A>G)
    #

    elif dmg_model == 'bi':
        if mut in right and state == 'derived':
            return 'yes'
        if mut in left and state == 'ancestral':
            return 'yes'
        return 'no'

    elif dmg_model != 'none':
        raise Exception("Invalid damage model provided via --dmg")

    return 'no'

#
#	Fiat lux
#

def main(lib_type, out_prefix, alleles_file, snpinfo_file, dmg_model):

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

	# damage assesment
    if not hits.empty:
        hits['damage'] = hits.apply(assess_damage, axis=1, lib_type=lib_type, dmg_model=dmg_model)
    else:
        hits['damage'] = []

    hits.to_csv(f'{out_prefix}.calls', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Match to tree and filter Y-chromosome genotypes")
    parser.add_argument('--lib', choices=['ds', 'ss', 'both'], required=True, help="Library build: 'ds' for double-stranded, 'ss' for single-stranded, 'both' if both were merged in bam")
    parser.add_argument('--dmg', choices=['none','naive','uni','bi'], required=True, help="Damage model to assess damage: 'naive' to remove all C>T (and G>A depending on --lib), 'uni' for unidirectional damage model (only derived), 'bi' for bidirectional damage (both derived and ancestral)")
    parser.add_argument('--out', required=True, help="Output file prefix")
    parser.add_argument('--alleles', required=True, help="Path to the alleles file")
    parser.add_argument('--snpinfo', required=True, help="Path to the SNP info file")
    args = parser.parse_args()

    main(args.lib, args.out, args.alleles, args.snpinfo, args.dmg)
