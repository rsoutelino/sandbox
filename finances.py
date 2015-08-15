import os, sys,  yaml
from glob import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging 
import datetime as dt

# TO-DO list
#     - put wildcards for things that enumerate: YouTube, Tank
#     - hook directly to google drive instead of keeping exporting

pd.set_option('display.max_rows', 10)
pd.set_option('expand_frame_repr', False)
mpl.rcParams['figure.facecolor'] = 'w'
# mpl.style.use('ggplot')
DATADIR = '/data/finances'
METADATADIR = '/source/sandbox'



def parse1(obstime): 
    return dt.datetime.strptime(obstime, "%d/%m/%Y")


class MonthStatement(object):
    def __init__(self, year, month):
        self.categories = yaml.load(open(os.path.join(METADATADIR, 'categories.yml')))
        self.data = os.path.join(DATADIR, 'statement_%s%02d.csv' %(year, month))
        self.df = pd.read_csv(self.data, index_col=0, date_parser=parse1)
        self.categorize()
        self.month = "%s%02d" %(year, month)

    def categorize(self):
        # removing ingnorables
        try:
            self.df = self.df[self.df.category != 'ignore']
            self.df = self.df.fillna(-999)
        except TypeError:
            self.df.category = -999

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
        notFoundList = False
        for idx, ben in enumerate(self.df.beneficiary):
            foundBen = False
            if idx not in saveidx:
                # allindexes.append(idx)
                for key1, dict1 in self.categories.items():
                    for key2, value2 in dict1.items():
                        if ben in value2:
                            foundBen = True
                            # indexes.append(idx)
                            # print idx, " Beneficiary %s is contained in cat2 = %s" %(ben, key2)
                            c1.append(key1)
                            c2.append(key2)
            else:
                foundBen = True
                c1.append(self.df.category[idx])
                c2.append(self.df.category2[idx])

            if not foundBen:
                notFoundList = True
                logging.error('No category found. Categorize "%s" and retry.' %ben)


        if notFoundList:
            sys.exit('Categorization FAILED. Categorize the items above and retry.')

        self.df.category = c1
        self.df.category2 = c2

    def get_total_expenses(self):
        df = self.df
        total_debit = df.amount[df.amount < 0].sum() * -1 
        savings     = df.amount[df.category == 'savings'].sum() * -1
        return total_debit - savings

    def get_by_category(self, category):
        df = self.df
        result = df.amount[df.category == category].sum()
        return np.abs(result)


class Period(object):
    """docstring for Period"""
    def __init__(self, start, end):
        """
        start   ::  YYYYMM   [Ex: 201503]
        end     ::  YYYYMM   [Ex: 201507] 
        """
        self.start = dt.datetime.strptime(start, '%Y%m')
        self.end   = dt.datetime.strptime(end, '%Y%m')
        self.daterange = pd.date_range(self.start, self.end, freq='M')
        self.daterange = [date - dt.timedelta(days=15) for date in self.daterange]

        monthList = []
        for month in self.daterange:
            df = MonthStatement(month.year, month.month)
            monthList.append(df)

        self.monthList = monthList

    def get_summary(self):
        total_expenses = [month.get_total_expenses() for month in self.monthList]
        summary = pd.DataFrame(index=self.daterange, data=dict(total_exp=total_expenses))

        for cat in self.monthList[0].categories:
            summary[cat] = [month.get_by_category(cat) for month in self.monthList]

        summary.pop('uncategorized')
        self.summary = summary

    def plot_summary(self):
        mpl_style_dicts = ['#e34a33', 
                           '#2ca25f', 
                           '#8856a7', 
                           '#43a2ca', 
                           '#c994c7', 
                           '#756bb1', 
                           '#d95f0e', 
                           '#636363', 
                           '#31a354', 
                           '#dd1c77'] 

        self.get_summary()
        self.summary.plot(subplots=True, sharex=True, layout=(5, 2), 
                          kind='area', stacked=False, style=mpl_style_dicts,
                          figsize=(17,12), grid=False)
        plt.tight_layout()












