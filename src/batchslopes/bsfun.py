import numpy as np
import pandas as pd
import os
import datetime

############################
#
############################
def MinMaxRange(x):
    '''
    Find the indices between minimum and maximum OD.
    
    Arguments:
        x:  numpy vector, OD measurements
        
    Return:
            numpy vector, indices of x
    '''   
    return np.arange(np.argmin(x), np.argmax(x)+1)
############################
#
############################
def MakeBins(x,partitions=2):
    '''
    The input vector is binned into equal sized partitions.

    Arguments:
        x:          numpy vector, OD measurements
        partitions: integer, number of bins
    Returns:
        bins:       numpy array, indices of start (column: 0) and end (column: 1) of bins
                    0, if fewer than three measurements would occupy the bins
        
    '''
    mb = int(len(x)/partitions)
    if mb > 0:
        mall = np.arange(0,len(x)+1,mb)
        bins = np.vstack([mall[:-1],mall[1:]]).T
    else:
        bins = 0
        
    return bins
############################
#
############################
def SlopeCalc(t, x, bins):
    '''
    Calculation of linear regression on ln(OD) and time.
    
    Arguments:
        t,x:        numpy vector, OD and time measurements
        bins:       numpy array, indices of start (column: 0) and end (column: 1) of bins
        
    Returns:
        (m,c,r2):   m: float, slope; c: float, y-intersection; r2: float, correlation coefficient
        myRange:    numpy vector, updated indices of x used for linear regression
    '''
    myRange = np.arange(bins[0],bins[1])
    MaxRange = MinMaxRange(x[myRange])
    myRange = myRange[MaxRange]
    if len(myRange) < 4:
#         print('too few data points in partition.')
        return False
    else:            
        Time = np.vstack([t[myRange], np.ones(len(myRange))]).T
        m,c = np.linalg.lstsq(Time, np.log(x[myRange]), rcond=None)[0]
        r2 = np.corrcoef(t[myRange],np.log(x[myRange]))[0][1]
        return ((m, c, r2), myRange)
############################
#
############################    
def DetectR2MaxSingle(t, x, partitions):
    '''
    Iterative partitioning of the data range to find the bin with the highest correlation coefficient for exponential growth.
    
    Arguments:
        t,x:        numpy vector, OD and time measurements
        partitions: integer, number of bins
    '''
    # Only checking region between min and max OD values
    MaxRange = MinMaxRange(x)
    ttest, xtest = t[MaxRange], x[MaxRange]
    
    # Initiating loop parameters
    counter = 0
    R2_ref = 0
    R2_tst = .0000001
    Slope_ref = 0
    Slope_tst = .000001
    # Generating new partitions in the data until the correlation coefficient of the new sub-bins is lower than the combined parent bin
    while R2_tst > R2_ref or R2_tst > .95: 
        # overwriting reference with current solution
        Slope_ref = Slope_tst
        R2_ref = R2_tst    
        # partitioning new range
        bins = MakeBins(xtest,partitions)
        # Return False if too few data points remain for sensible binning (<4) 
        if np.sum(bins) == 0:
#             print('bins is false')
            return False
        # if there are too many partitions in the beginning we reduce it further down into the loop
        if partitions > 3: 
            partitions = np.ceil(partitions/2)
            
        # calculating the regression details
        result = np.array([SlopeCalc(ttest, xtest, mybin) for mybin in bins])
        # finding subsequences with too little sample points and removing them from further analysis
        Posidx = [myres!=False for myres in result]
        result = result[Posidx]
        # if enough data points could be analysed
        if result.any():
            # extracting details for each partition
            RegrList = list()
            RangeList = list()
            [[RegrList.append(indi[0]),RangeList.append(indi[1])] for indi in result]
            RegrList = np.array(RegrList)
            RangeList = np.array(RangeList)
            slope_MaxIdx = np.argmax(RegrList[:,0])
            Slope_tst = RegrList[slope_MaxIdx,0]
            R2_tst = RegrList[slope_MaxIdx,2]   
            Ycorr = RegrList[slope_MaxIdx,1]
            # selecting the range with the previous highest slope
            myRange = RangeList[slope_MaxIdx] # np.arange(bins[slope_MaxIdx, 0], bins[slope_MaxIdx,1])
            xtest = xtest[myRange]
            ttest = ttest[myRange]
            # If we would decide to break while loop based on the number of iterations
            counter +=1 
            
            # Storing best solution
            if R2_tst > R2_ref:
                myResult = {'R2':R2_tst, 'Slope': Slope_tst, 'ycorrect':Ycorr, 'time':ttest, 'OD':xtest}
        else:
            if 'myResult' in locals():
                return myResult
            else:
                myResult = False
                return myResult
        
    return myResult
############################
#
############################
def CorrectedOD(GV,expo1,expo2):
    '''
        Correction of the OD
    '''
    od = (GV ** expo1) * (np.exp(expo2))
    return od
############################
#
############################
def StatsCalc(df, repnum):
    '''
        mean and standard deviation for replicates
        taken from df and repnum
    '''
    r = int(len(df.columns)/repnum)
    mean = pd.DataFrame()
    stdev = pd.DataFrame()
    for idxr in range(r):
        idxr += 1
        idxnew = idxr * repnum 
        mean[str(idxr)] = (df.iloc[:,(idxnew-repnum):(idxnew)].mean(1))
        stdev[str(idxr)] = (df.iloc[:,(idxnew-repnum):(idxnew)].std(1))
    
    return mean, stdev
############################
#
############################
def Batchfilenames(dirname):
    '''
    only relevant csv files!
    '''
    files = os.listdir(dirname)    
    csv_files = list(filter(lambda f: f.endswith('.csv'), files))
    myFiles = list()
    for idx in range(len(csv_files)):
        myFiles.append(os.path.join('data', csv_files[idx]))
    return(myFiles)
############################
#
############################
def Fileparams(Source):
    '''
    if synthetic -> all parameter output (skiprows, etc)
    else -> take all csv's
    creating a dictionary
    '''
    if Source == 0:
    # 0 means synthetic, 1 means real
        myparams= {'skiprows': 0, 'decimal' : '.', 'TimeUnit' : 1}
    else:
        myparams= {'skiprows': 23, 'decimal' : ',', 'TimeUnit' : 60}
    return(myparams)
############################
#
############################
def Foldertest(Name):
    extension = os.path.splitext(Name)[1]
    if extension == '':
        return 1
    else:
        return 0
############################
#
############################
def createExcelFiles(tk, means, stdev, myCols, mu_list, r2_list, repnum, df, OD, myFiles):

    today = datetime.date.today()
    
    myr2table = pd.DataFrame()
    myr2table['Wells'] = myCols
    myr2table['Mu Value'] = mu_list
    myr2table['R2 Value'] = r2_list

    myr2table2 = pd.DataFrame()
#     myr2table2['Wells'] = myCols
    myr2table2['Means of Mu'] = list(np.mean(np.array(mu_list).reshape(-1, repnum), axis = 1))
    myr2table2['StandardDeviation of Mu'] = list(np.std(np.array(mu_list).reshape(-1, repnum), axis = 1))
    myr2table2['Means of R2'] = list(np.mean(np.array(r2_list).reshape(-1, repnum), axis = 1))
    myr2table2['StandardDeviation of R2'] = list(np.std(np.array(r2_list).reshape(-1, repnum), axis = 1))
    myr2table2.index = np.arange(1, len(myr2table2)+1)
# df[['A1', 'A2','A3','A4', 'D5']].plot()

# int((len(df.columns))/repnum)
    if isinstance(myFiles, str):
        with pd.ExcelWriter(today.strftime("%y%m%d")+'_'+myFiles[:-4]+'_results.xlsx') as writer:
            df.to_excel(writer, sheet_name='1 raw data')
            OD.to_excel(writer, sheet_name='2 correctedOD')
            means.to_excel(writer, sheet_name='3 means')
            stdev.to_excel(writer, sheet_name='4 stdev')
            myr2table.to_excel(writer, sheet_name='5 Growthrates', index = False)
            myr2table2.to_excel(writer, sheet_name='6 means&stdev of Growthrates', index = False)
    else:
        with pd.ExcelWriter(today.strftime("%y%m%d")+'_'+myFiles[tk][5:-4]+'_results.xlsx') as writer:
            df.to_excel(writer, sheet_name='1 raw data')
            OD.to_excel(writer, sheet_name='2 correctedOD')
            means.to_excel(writer, sheet_name='3 means')
            stdev.to_excel(writer, sheet_name='4 stdev')
            myr2table.to_excel(writer, sheet_name='5 Growthrates', index = False)
            myr2table2.to_excel(writer, sheet_name='6 means&stdev of Growthrates', index = False)