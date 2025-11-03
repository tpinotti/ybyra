import pandas as pd
import sys

if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help'):
    print("Usage: python ynvert.py <aggregate.yplace>")
    sys.exit()

input = sys.argv[1]
output = f"inv.{input}"

df = pd.read_csv(input, sep='\t')
df['tree_path'] = df['tree_path'].apply(lambda x: '<'.join(x.split('<')[::-1]))
df.to_csv(output, sep='\t', index=False)
