#!/usr/bin/env python
import pandas as pd
import os 
import sys
aa_tbl_jned = pd.read_table(sys.argv[2])
aa_tbl_jned.index = aa_tbl_jned['#IID']
cur_col = aa_tbl_jned[sys.argv[1]]
cur_col.to_csv(sys.argv[3] + sys.argv[1],index_label='#IID',sep='\t',na_rep='NA')

