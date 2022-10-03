import matplotlib as mpl
mpl.use('TkAgg')
import numpy as np
import pandas as pd
from ModifiedFedergruenAndZheng2 import ModifiedFedergruenAndZheng
from DoubleRQ_ChuaFengSengaViswanathan_Continuous import SimulateDoubleRQ
from F_SimulationOptimization2 import SimulationOptimization


#function that gets optimal double RQ
def OptimalDoubleRQ(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r):
    return ModifiedFedergruenAndZheng(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r)

#function that gets near optimal sQsS
def NearOptimalsQsS(r_r, r_p, Q_r, Q_p, Lambda_d, Lambda_r, K_r, K_p,c_p, c_r, h_p, h_r, p, L_p, L_r):
    #from POMS Chicago PPT: s_m = r_r*, Q_m = Q_r*, s_r = r_p*, S_r = r_p* + Q_p*
    DeltaNeighborhood = 3
    TimeToSteadyState = 1000.0
    TimeUnits = 10000.0
    sQsS_SimulationDataFrame, sQsS_AverageDataFrame, sQsS_MinimumDataFrame = SimulationOptimization(r_r, r_p, Q_r, Q_p, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, DeltaNeighborhood, TimeToSteadyState, TimeUnits)
    return sQsS_SimulationDataFrame, sQsS_AverageDataFrame, sQsS_MinimumDataFrame

#function that connects to simulation of Double RQ
def SimulateRQ(r_r, r_p, Q_r, Q_p, Lambda_d, Lambda_r, K_r, K_p,c_p, c_r, h_p, h_r, p, L_p, L_r, seed):
    TimeToSteadyState = 1000.0
    TimeUnits = 10000.0
    return SimulateDoubleRQ(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r,seed, TimeToSteadyState, TimeUnits)

def WriteToMainDataFrame(Outputs, DataFrame, BaselinesQsSAverage, s_m, Q_m, s_r, S_r, sQsSAverage,Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p, L_p, L_r,RQGreaterThansQsSMinDifference,RQGreaterThansQsSMaxDifference,RQLessThanOrEqualTosQsSMinDifference, RQLessThanOrEqualTosQsSMaxDifference,PercentLessThan):
        DataFrame["sQsS Average Using Double RQ"] = BaselinesQsSAverage
        DataFrame["sQsS Average"] = sQsSAverage
        DataFrame["s_m"] = s_m
        DataFrame["Q_m"] = Q_m
        DataFrame["s_r"] = s_r
        DataFrame["S_r"] = S_r
        DataFrame["(sQsS - Average)/sQsS"] = (sQsSAverage - DataFrame["Average Result"])/sQsSAverage
        DataFrame["(BaselinesQsS - Average)/BaselinesQsS"] = (BaselinesQsSAverage - DataFrame["Average Result"]) / BaselinesQsSAverage
        DataFrame["(sQsS - BaselinesQsS)/sQsS"] = (sQsSAverage - BaselinesQsSAverage) / sQsSAverage
        DataFrame["Percent of RQ <= sQsS"] = PercentLessThan
        DataFrame["Lambda_d"] = Lambda_d
        DataFrame["Lambda_r"] = Lambda_r
        DataFrame["K_r"] = K_r
        DataFrame["K_p"] = K_p
        DataFrame["c_p"] = c_p
        DataFrame["c_r"] = c_r
        DataFrame["h_p"] = h_p
        DataFrame["h_r"] = h_r
        DataFrame["p"] = p
        DataFrame["L_p"] = L_p
        DataFrame["L_r"] = L_r

        ForReturning = Outputs.append(DataFrame)
        return ForReturning

def Main(FilePath, FileCounter):
    rowcounter = 1
    Outputs_sQsS = pd.DataFrame()
    Outputs_DoubleRQ = pd.DataFrame()
    Outputs = pd.DataFrame()

    DoubleRQSimList = []
    SeedListRQ = []
    SimulationParameters = pd.read_csv(FilePath)

    #iterate through each of the simulation parameters
    for index,row in SimulationParameters.iterrows():

        #Clearing The List
        del DoubleRQSimList[:]
        del SeedListRQ[:]

        ID = row["ID"]
        Lambda_d = row["Lambda_d"]
        Lambda_r = row["Lambda_r"]
        K_r = row["K_r"]
        K_p = row["K_p"]
        c_p = row["c_p"]
        c_r = row["c_r"]
        h_p = row["h_p"]
        h_r = row["h_r"]
        p = row["p"]
        L_p = row["L_p"]
        L_r = row["L_r"]
        print ID
        #obtain optimal Double RQ parameters and Simulation Results
        r_r, r_p, Q_r, Q_p, Z = OptimalDoubleRQ(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r)

        randomnumbers = pd.read_csv("RandomNumbers.csv")
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            DoubleRQSimList.append(SimulateRQ(r_r, r_p, Q_r, Q_p, Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p, L_p, L_r, seed))
            SeedListRQ.append(seed)


        DoubleRQSimAverage = sum(DoubleRQSimList)/float(len(DoubleRQSimList))
        DoubleRQ_ResultDataFrame = pd.DataFrame()
        DoubleRQ_ResultDataFrame = DoubleRQ_ResultDataFrame.append({"ID": rowcounter, "r_r": r_r, "r_p": r_p, "Q_r": Q_r, "Q_p": Q_p, "Average Result": DoubleRQSimAverage, "Theoretical Double RQ": Z}, ignore_index=True)

        # Obtain near optimal sQsS parameters and Simulation Results
        sQsS_SimulationDataFrame, sQsS_AverageDataFrame, sQsS_MinimumDataFrame = NearOptimalsQsS(r_r, r_p, Q_r, Q_p, Lambda_d, Lambda_r, K_r, K_p,c_p, c_r, h_p, h_r, p, L_p, L_r)
        sQsSAverage = sQsS_MinimumDataFrame["Average Result"].item()
        s_m = sQsS_MinimumDataFrame["s_m"].item()
        Q_m = sQsS_MinimumDataFrame["Q_m"].item()
        s_r = sQsS_MinimumDataFrame["s_r"].item()
        S_r = sQsS_MinimumDataFrame["S_r"].item()

        #Compare Simulation Results of Double RQ and sQsS
        ComparisonRQsQsSDataFrame = pd.DataFrame({"Double RQ Results": DoubleRQSimList})
        SeedDataFrameRQ = pd.DataFrame({"SeedRQ": SeedListRQ})
        ComparisonRQsQsSDataFrame["SeedRQ"] = pd.Series(SeedDataFrameRQ["SeedRQ"])

        PlaceholderDataFrame = sQsS_SimulationDataFrame[(sQsS_SimulationDataFrame["s_m"] == s_m) & (sQsS_SimulationDataFrame["Q_m"] == Q_m) & (sQsS_SimulationDataFrame["s_r"] == s_r) & (sQsS_SimulationDataFrame["S_r"] == S_r)]
        PlaceholderDataFrame = PlaceholderDataFrame.reset_index(drop=True)

        s_mBaseLine = r_p
        Q_mBaseLine = Q_p
        s_rBaseLine = r_r
        S_rBaseLine = r_r + Q_r

        BaselinesQsSAverageDataFrame = sQsS_AverageDataFrame[(sQsS_AverageDataFrame["s_m"] == s_mBaseLine) & (sQsS_AverageDataFrame["Q_m"] == Q_mBaseLine) & (sQsS_AverageDataFrame["s_r"] == s_rBaseLine) & (sQsS_AverageDataFrame["S_r"] == S_rBaseLine)]
        BaseLinesQsSAverage = BaselinesQsSAverageDataFrame["Average Result"].item()

        ComparisonRQsQsSDataFrame["Simulation Result"] = pd.Series(PlaceholderDataFrame["Simulation Result"])
        ComparisonRQsQsSDataFrame["Seed"] = pd.Series(PlaceholderDataFrame["Seed"])



        #Get Minimum and Maximum Differences from Simulation Results
        ComparisonRQsQsSDataFrame["RQ - sQsS"] = ComparisonRQsQsSDataFrame["Double RQ Results"] - ComparisonRQsQsSDataFrame["Simulation Result"]
        # Get how many have RQ <= sQsS
        ComparisonRQsQsSDataFrame["RQ Less Than or Equal To?"] = np.where((ComparisonRQsQsSDataFrame["Double RQ Results"] <= ComparisonRQsQsSDataFrame["Simulation Result"]), True, False)

        RQGreaterThansQsS = ComparisonRQsQsSDataFrame.loc[ComparisonRQsQsSDataFrame["RQ Less Than or Equal To?"] == False]
        RQLessThanOrEqualTosQsS = ComparisonRQsQsSDataFrame.loc[ComparisonRQsQsSDataFrame["RQ Less Than or Equal To?"] == True]

        RQGreaterThansQsSMinDifference = np.nan
        RQGreaterThansQsSMaxDifference = np.nan

        #sQss > RQ Results
        if len(RQGreaterThansQsS) > 0:
            RQGreaterThansQsSMinimumDifferenceDataFrame = RQGreaterThansQsS.ix[RQGreaterThansQsS["RQ - sQsS"].idxmin()]
            RQGreaterThansQsSMaximumDifferenceDataFrame = RQGreaterThansQsS.ix[RQGreaterThansQsS["RQ - sQsS"].idxmax()]
            RQGreaterThansQsSMinDifference = RQGreaterThansQsSMinimumDifferenceDataFrame["RQ - sQsS"].item()
            RQGreaterThansQsSMaxDifference = RQGreaterThansQsSMaximumDifferenceDataFrame["RQ - sQsS"].item()

        RQLessThanOrEqualTosQsSMinDifference = np.nan
        RQLessThanOrEqualTosQsSMaxDifference = np.nan

        #sQsS <= RQ Results
        if len(RQLessThanOrEqualTosQsS) > 0:
            RQLessThanOrEqualTosQsSMinimumDifferenceDataFrame = RQLessThanOrEqualTosQsS.ix[RQLessThanOrEqualTosQsS["RQ - sQsS"].idxmin()]
            RQLessThanOrEqualTosQsSMaximumDifferenceDataFrame = RQLessThanOrEqualTosQsS.ix[RQLessThanOrEqualTosQsS["RQ - sQsS"].idxmax()]
            RQLessThanOrEqualTosQsSMinDifference = RQLessThanOrEqualTosQsSMinimumDifferenceDataFrame["RQ - sQsS"].item()
            RQLessThanOrEqualTosQsSMaxDifference = RQLessThanOrEqualTosQsSMaximumDifferenceDataFrame["RQ - sQsS"].item()


        #Get Percent of Simulation RQ <= sQsS
        NumRows = ComparisonRQsQsSDataFrame.shape[0]
        NumTrue = len(ComparisonRQsQsSDataFrame[ComparisonRQsQsSDataFrame["RQ Less Than or Equal To?"] == True])
        PercentLessThan = NumTrue/float(NumRows)

        #Combine into DataFrames
        Outputs = WriteToMainDataFrame(Outputs, DoubleRQ_ResultDataFrame, BaseLinesQsSAverage, s_m, Q_m, s_r, S_r, sQsSAverage, Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p, L_p, L_r, RQGreaterThansQsSMinDifference, RQGreaterThansQsSMaxDifference,RQLessThanOrEqualTosQsSMinDifference, RQLessThanOrEqualTosQsSMaxDifference, PercentLessThan)

        #next set of parameters
        rowcounter += 1

    # Output into CSV Files
    Outputs.to_csv("Simulation_Output_" + str(FileCounter) + ".csv", encoding='utf-8', index=False)