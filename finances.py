import os, sys,  yaml, urllib2
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
#     - create functionality to plot savings time series

pd.set_option('display.max_rows', 30)
pd.set_option('expand_frame_repr', False)
mpl.rcParams['figure.facecolor'] = 'w'
# mpl.style.use('ggplot')
DATADIR = '/data/finances'
METADATADIR = '/source/sandbox'



def parse1(obstime): 
    return dt.datetime.strptime(obstime, "%d/%m/%Y")

def get_currency(currency_in, currency_out):
    request = 'http://finance.yahoo.com/d/quotes.csv?e=.csv&f=sl1d1t1&s=%s%s=X' %(currency_in, currency_out)
    req = urllib2.Request(request)
    response = urllib2.urlopen(req)
    factor = float( response.read().split(',')[1] )
    return factor

def convert(value, factor):
    result = value * factor
    return result

def convert_online(value, currency_in, currency_out):
    factor = get_currency(currency_in, currency_out)
    result = value * factor
    return result


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
        total_debit = np.abs( df.amount[df.amount < 0].sum() ) 
        savings     = np.abs( df.amount[df.category == 'savings'].sum() )
        return total_debit - savings

    def get_basic_expenses(self):
        df = self.df
        basic = df.amount[df.category == 'groceries'].sum() +\
                df.amount[df.category == 'restaurant'].sum() +\
                df.amount[df.category == 'car'].sum() +\
                df.amount[df.category == 'living'].sum()
        return np.abs(basic)         

    def get_soft_expenses(self):
        df = self.df
        soft = df.amount[df.category == 'shopping'].sum() +\
                df.amount[df.category == 'wellness'].sum() +\
                df.amount[df.category == 'entertainement'].sum()
        return np.abs(soft) 

    def get_by_category(self, category):
        df = self.df
        result = df.amount[df.category == category].sum()
        return np.abs(result)

    def get_by_category2(self, category2):
        df = self.df
        result = df.amount[df.category2 == category2].sum()
        return np.abs(result)


class Period(object):
    """docstring for Period"""
    def __init__(self, start, end):
        """
        start   ::  YYYYMM   [Ex: 201503]
        end     ::  YYYYMM   [Ex: 201507] 
        """
        self.goals = yaml.load(open(os.path.join(METADATADIR, 'goals.yml')))
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
        basic_expenses = [month.get_basic_expenses() for month in self.monthList]
        soft_expenses  = [month.get_soft_expenses()  for month in self.monthList]

        summary = pd.DataFrame(index=self.daterange, 
                               data=dict(total_expenses=total_expenses, 
                                         basic_expenses=basic_expenses,
                                         soft_expenses=soft_expenses))

        for cat in self.monthList[0].categories:
            summary[cat] = [month.get_by_category(cat) for month in self.monthList]

        rescue = [month.get_by_category2('rescue') for month in self.monthList]
        summary.savings = summary.savings - rescue
        self.summary = summary

    def get_by_subcategory(self, cat):
        data = [month.get_by_category2(cat) for month in self.monthList]
        df = pd.DataFrame(index=self.daterange, data={cat: data})
        return df
        

    def plot_summary(self, scale=False):
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
        dfplot = self.summary.copy()
        dfplot.pop('uncategorized')
        dfplot.pop('savings')
        dfplot.pop('total_expenses')

        axes = dfplot.plot(subplots=True, sharex=True, layout=(5, 2), legend=False,
                           kind='bar', stacked=False, style=mpl_style_dicts,
                           figsize=(17,12), grid=True, width=0.8, alpha=0.5)
        ax = plt.gca()
        ticks = [dt.strftime('%b/%Y') for dt in dfplot.index.to_pydatetime()]
        ax.set_xticklabels(ticks, rotation=30)

        a, b = axes.shape
        for i in range(a):
            for j in range(b):
                ax = axes[i,j]
                title = ax.get_title()
                dfplot['value'] = self.goals[title]
                ax.plot(dfplot.index, dfplot.value) # DON'T KNOW WHY THIS IS NOT WORKING



                # ax.plot([dfplot.index[0], dfplot.index[-1]], [self.goals[title], self.goals[title]], 
                        # 'k', linewidth=3, alpha=0.2)
                if scale:
                    if ax.get_title() != 'income':
                        ax.set_ylim([0, 3500])
            
        plt.tight_layout()
        plt.show()

    def plot_budget(self):
        self.get_summary()
        self.budget = self.summary.income - (self.summary.total_expenses)
        self.budget.plot(kind='bar', color='#756bb1', stacked=False, 
                         alpha=0.5, width=0.8, figsize=(17,6))
        ax = plt.gca()
        ticks = [dt.strftime('%b/%Y') for dt in self.budget.index.to_pydatetime()]
        ax.set_xticklabels(ticks, rotation=30)

    def plot_by_subcategory(self, cat):
        df = self.get_by_subcategory(cat)
        self.budget = self.summary.income - (self.summary.total_expenses)
        df.plot(kind='bar', color='#756bb1', stacked=False, 
                alpha=0.5, width=0.8, figsize=(17,6))
        ax = plt.gca()
        ticks = [dt.strftime('%b/%Y') for dt in self.budget.index.to_pydatetime()]
        ax.set_xticklabels(ticks, rotation=30)


class Projection(object):
    """docstring for Projection"""
    def __init__(self, start=None, years=5):
        self.goals = yaml.load(open(os.path.join(METADATADIR, 'goals.yml')))
        self.sources = self.goals['sources']
        self.sinks   = self.goals['sinks']

        today = dt.date.today()
        start = start or dt.datetime(today.year, today.month, 1)
        end   = dt.datetime(today.year + years, today.month, 1)
        self.daterange = pd.date_range(start, end, freq='M')
        self.daterange = [date - dt.timedelta(days=15) for date in self.daterange]
        self.exchange_factor = get_currency('BRL', 'NZD') 

    def get_source_sink_from_goals(self, yearly_trip=False, dad_deduction=False, cef_deduction=False):
        source, sink = [], []
        for month in self.daterange:
            sources, sinks = 0, 0
            for key, src in self.sources.items():
                # if key == 'juliana_wec':
                    # import pdb; pdb.set_trace()
                # stablishing the value based on the dates this source is valid
                if src.has_key('starts') and src.has_key('expires'):
                    starts = dt.datetime.strptime(str(src['starts']), '%Y%m%d')
                    ends   = dt.datetime.strptime(str(src['expires']), '%Y%m%d')
                    if month < starts or month > ends:
                        value = 0
                    else:
                        value = src['value']
                elif src.has_key('starts'):
                    starts = dt.datetime.strptime(str(src['starts']), '%Y%m%d')
                    if month < starts:
                        value = 0
                    else:
                        value = src['value']
                elif src.has_key('expires'):
                    ends = dt.datetime.strptime(str(src['expires']), '%Y%m%d')
                    if month > ends:
                        value = 0
                    else:
                        value = src['value']
                else:
                    value = src['value']

                if src['currency'] != 'NZD':
                    sources += convert(value, self.exchange_factor)
                else:
                    sources += value
                
                del value

            for key, snk in self.sinks.items():
                # stablishing the value based on the dates this sink is valid
                if snk.has_key('starts') and snk.has_key('expires'):
                    starts = dt.datetime.strptime(str(snk['starts']), '%Y%m%d')
                    ends   = dt.datetime.strptime(str(snk['expires']), '%Y%m%d')
                    if month < starts or month > ends:
                        value = 0
                    else:
                        value = snk['value']
                elif snk.has_key('starts'):
                    starts = dt.datetime.strptime(str(snk['starts']), '%Y%m%d')
                    if month < starts:
                        value = 0
                    else:
                        value = snk['value']
                elif snk.has_key('expires'):
                    ends = dt.datetime.strptime(str(snk['expires']), '%Y%m%d')
                    if month > ends:
                        value = 0
                    else:
                        value = snk['value']
                else:
                    value = snk['value']

                if snk['currency'] != 'NZD':
                    sinks += convert(value, self.exchange_factor)
                else:
                    sinks += value
                
                del value

            if yearly_trip:
                if month.month == 12:
                    sinks += self.goals['extra']['yearly_trip']

            if dad_deduction:
                if month.month == 3 or month.month == 9:
                    sinks += self.goals['extra']['dad_deduction'] / 2 
            
            if cef_deduction:
                if month.month == 1:
                    sinks += self.goals['extra']['cef_deduction']


            source.append(sources)
            sink.append(sinks)

        return np.array(source), np.array(sink)


    def from_goals(self, yearly_trip=False, dad_deduction=False, 
                         cef_deduction=False):
        mpl_style_dicts = ['#99d8c9', 
                           '#fa9fb5', 
                           '#3182bd'] 
        source, sink = self.get_source_sink_from_goals(yearly_trip=yearly_trip,
                                                       dad_deduction=dad_deduction,
                                                       cef_deduction=cef_deduction)
        self.df = pd.DataFrame(index=self.daterange, data={'source': source})
        self.df['sink']  = sink
        self.df['gains'] = source - sink 

        plt.figure(figsize=(15,8))
        ax1 = plt.subplot(211)
        ax1.set_title('Projection from stablished goals', fontweight='bold')
        self.df.plot(ax=ax1, kind='area', stacked=False, style=mpl_style_dicts)
        
        self.df['goods'] = np.cumsum(source - sink) 
        self.df.goods += self.goals['extra']['initial_balance']

        ax2 = plt.subplot(212)
        self.df.goods.plot(ax=ax2, kind='area', stacked=False, 
                           color=mpl_style_dicts[-1], legend=True)
        print self.df.tail()
        plt.show()

        
        






                

    def from_averages(self):
        pass

        
class Savings(object):
    """docstring for Savings"""
    def __init__(self):
        data1  = os.path.join(DATADIR, 'online_saver.csv')
        data2  = os.path.join(DATADIR, 'online_bonus_saver.csv')
        df1    = pd.read_csv(data1, index_col=0, date_parser=parse1)
        df2    = pd.read_csv(data2, index_col=0, date_parser=parse1)
        goods1 = np.cumsum(df1.Amount) 
        goods2 = np.cumsum(df2.Amount)

        self.df = pd.DataFrame(index=df1.index, data={'online_saver': goods1})
        self.df = self.df.resample('M')
        self.df['bonus_saver'] = goods2.resample('M')
        
        



        









