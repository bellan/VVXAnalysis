#!/usr/bin/env python3

from json import load

with open('data/systematics.json') as fin:
    syst_values = load(fin)

variable = 'mWZG'  # 'mZZG'
sample0 = 'WZGTo3LNuG'  # 'ZZGTo4LG'

print('#'*5, variable, '#'*5)

systs = syst_values[sample0][variable].keys()
print('{:16s}'.format('PROCESS'), end='\t')
for syst in systs:
    print('{:>11.11s}'.format(syst), end='\t')
print()

for sample in syst_values:
    out_line = '{:16.16s} '.format(sample)
    for syst in systs:
        shifts = syst_values[sample][variable][syst]
        up = shifts['up']
        dn = shifts['dn']
        if( abs(up-1) > 1.5 or abs(dn-1) > 1.5 ):
            out_line += '\t{:12.12s}'.format('largeShift')
        else:
            out_line += '\t{:4.4f}'.format(1 + abs(up - dn) / 2)
            #out_line += '{:2.2f} / {:2.2f}  '.format(up, dn)
    print(out_line)
