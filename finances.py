import os, yaml
from glob import glob
import numpy as np
import pandas as pd 

pd.set_option('display.max_rows', 10)
pd.set_option('expand_frame_repr', False)
DATADIR = '/data/finances'
METADATADIR = '/source/sandbox'


class MonthStatement(object):
    def __init__(self, year, month):
        self.categories = yaml.load(open(os.path.join(METADATADIR, 'categories.yml')))

        self.data = os.path.join(DATADIR, 'statement_%s%02d.csv' %(year, month))
        self.df = pd.read_csv(self.data, index_col=0)
        # self.categorize()

    def categorize(self):
        # removing ingnorables
        self.df = self.df[self.df.category != 'ignore']
        self.df = self.df.fillna(-999)

        self.df['category2'] = -999

        # atributing pre-determined categories
        cat1, cat2 = [], []
        saveidx = []

        for idx, cat in enumerate(self.df.category):
            if cat != -999:
                cat1.append(cat.split(',')[0])
                cat2.append(cat.split(',')[1])
                saveidx.append(idx)
            else:
                cat1.append(-999)
                cat2.append(-999)

        self.df.category = cat1
        self.df.category2 = cat2


        c1, c2 = [], []
        # discovering categories based on beneficiaries
        # allindexes, indexes = [], []
        for idx, ben in enumerate(self.df.beneficiary):
            if idx not in saveidx:
                # allindexes.append(idx)
                for key1, dict1 in self.categories.items():
                    for key2, value2 in dict1.items():
                        if ben in value2:
                            # indexes.append(idx)
                            # print idx, " Beneficiary %s is contained in cat2 = %s" %(ben, key2)
                            c1.append(key1)
                            c2.append(key2)
            else:
                c1.append(self.df.category[idx])
                c2.append(self.df.category2[idx])

        self.df.category = c1
        self.df.category2 = c2

    def get_total_expenses(self):
        df = self.df
        total_debit = df.amount[df.amount < 0].sum() * -1 
        savings     = df.amount[df.category == 'savings'].sum() * -1
        return total_debit - savings














