import matplotlib as mpl
mpl.use('TkAgg')
import pandas as pd
from vanDerLaan_sQsS_Simulation_Continuous import SimulatesQsS

def SimWriteToDataFrame(DataFrame, s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r,seed, rowcounter, TimeToSteadyState, TimeUnits):
    ForReturning = DataFrame.append({"ID": rowcounter, "s_m": s_m, "Q_m": Q_m, "s_r": s_r, "S_r": S_r, "Simulation Result":SimulatesQsS(s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, TimeToSteadyState, TimeUnits), "Seed": seed}, ignore_index=True)
    return ForReturning

def AverageWriteToDataFrame(AverageDataFrame, FilteredDataFrame, s_m, Q_m, s_r, S_r, rowcounter):
    #print s_m, Q_m, s_r, S_r
    ForReturning = AverageDataFrame.append({"ID": rowcounter, "s_m": s_m, "Q_m": Q_m, "s_r": s_r, "S_r": S_r, "Average Result": FilteredDataFrame["Simulation Result"].mean()}, ignore_index=True)
    return ForReturning

def SimulationOptimization(r_r, r_p, Q_r, Q_p, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, DeltaNeighborhood, TimeToSteadyState, TimeUnits):
    randomnumbers = pd.read_csv("RandomNumbers.csv")
    s_m = r_p
    Q_m = Q_p
    s_r = r_r
    S_r = r_r + Q_r
    print s_m, Q_m, s_r, S_r

    #change s_m
    counter = 0
    rowcounter = 1

    SimulationDataFrame_s_m = pd.DataFrame()

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            #if s_m - counter >= 0:
            SimulationDataFrame_s_m = SimWriteToDataFrame(SimulationDataFrame_s_m, s_m - counter, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
            rowcounter += 1

            if s_m + counter <= s_r:
                SimulationDataFrame_s_m = SimWriteToDataFrame(SimulationDataFrame_s_m, s_m + counter, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

        counter += 1


    # Get the Average Cost of Simulated SamplePaths
    counter = 0
    rowcounter = 1

    AverageDataFrame_s_m = pd.DataFrame()

    while counter <= DeltaNeighborhood:
        #if s_m - counter >= 0:
        s_mMinusCounter = SimulationDataFrame_s_m[(SimulationDataFrame_s_m["s_m"] == s_m - counter) & (SimulationDataFrame_s_m["Q_m"] == Q_m) & (SimulationDataFrame_s_m["s_r"] == s_r) & (SimulationDataFrame_s_m["S_r"] == S_r)]
        AverageDataFrame_s_m = AverageWriteToDataFrame(AverageDataFrame_s_m, s_mMinusCounter, s_m - counter, Q_m, s_r, S_r, rowcounter)
        rowcounter += 1

        if s_m + counter <= s_r:
            s_mPlusCounter = SimulationDataFrame_s_m[(SimulationDataFrame_s_m["s_m"] == s_m + counter) & (SimulationDataFrame_s_m["Q_m"] == Q_m) & (SimulationDataFrame_s_m["s_r"] == s_r) & (SimulationDataFrame_s_m["S_r"] == S_r)]
            AverageDataFrame_s_m = AverageWriteToDataFrame(AverageDataFrame_s_m, s_mPlusCounter, s_m + counter, Q_m, s_r, S_r, rowcounter)
            rowcounter += 1

        counter += 1

    MinimumDataFrame_s_m = AverageDataFrame_s_m.ix[AverageDataFrame_s_m["Average Result"].idxmin()]



######################################################################################################################################################################################################
    counter = 0
    rowcounter = 1

    SimulationDataFrame_Q_m = pd.DataFrame()

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if Q_m - counter >= 1:
                SimulationDataFrame_Q_m = SimWriteToDataFrame(SimulationDataFrame_Q_m, s_m, Q_m - counter, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

            #Q_m + counter
            SimulationDataFrame_Q_m = SimWriteToDataFrame(SimulationDataFrame_Q_m, s_m, Q_m + counter, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
            rowcounter += 1

        counter += 1

    counter = 0
    rowcounter = 1

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if Q_m - counter >= 1:
                SimulationDataFrame_Q_m = SimWriteToDataFrame(SimulationDataFrame_Q_m, MinimumDataFrame_s_m["s_m"].item(), Q_m - counter, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1


            SimulationDataFrame_Q_m = SimWriteToDataFrame(SimulationDataFrame_Q_m, MinimumDataFrame_s_m["s_m"].item(), Q_m + counter, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
            rowcounter += 1

        counter += 1

    #Get Average Q_m
    counter = 0
    rowcounter = 1

    AverageDataFrame_Q_m = pd.DataFrame()

    while counter <= DeltaNeighborhood:

        if Q_m - counter >= 1:
            Q_mMinusCounter = SimulationDataFrame_Q_m[(SimulationDataFrame_Q_m["s_m"] == s_m) & (SimulationDataFrame_Q_m["Q_m"] == Q_m - counter) & (SimulationDataFrame_Q_m["s_r"] == s_r) & (SimulationDataFrame_Q_m["S_r"] == S_r)]
            AverageDataFrame_Q_m = AverageWriteToDataFrame(AverageDataFrame_Q_m, Q_mMinusCounter, s_m, Q_m - counter, s_r, S_r, rowcounter)
            rowcounter += 1

            Q_mMinusCounter = SimulationDataFrame_Q_m[(SimulationDataFrame_Q_m["s_m"] == MinimumDataFrame_s_m["s_m"].item()) & (SimulationDataFrame_Q_m["Q_m"] == Q_m - counter) & (SimulationDataFrame_Q_m["s_r"] == s_r) & (SimulationDataFrame_Q_m["S_r"] == S_r)]
            AverageDataFrame_Q_m = AverageWriteToDataFrame(AverageDataFrame_Q_m, Q_mMinusCounter, MinimumDataFrame_s_m["s_m"].item(), Q_m - counter, s_r, S_r, rowcounter)
            rowcounter += 1


        Q_mPlusCounter = SimulationDataFrame_Q_m[(SimulationDataFrame_Q_m["s_m"] == s_m) & (SimulationDataFrame_Q_m["Q_m"] == Q_m + counter) & (SimulationDataFrame_Q_m["s_r"] == s_r) & (SimulationDataFrame_Q_m["S_r"] == S_r)]
        AverageDataFrame_Q_m = AverageWriteToDataFrame(AverageDataFrame_Q_m, Q_mPlusCounter, s_m, Q_m + counter, s_r, S_r, rowcounter)
        rowcounter += 1

        Q_mPlusCounter = SimulationDataFrame_Q_m[(SimulationDataFrame_Q_m["s_m"] == MinimumDataFrame_s_m["s_m"].item()) & (SimulationDataFrame_Q_m["Q_m"] == Q_m + counter) & (SimulationDataFrame_Q_m["s_r"] == s_r) & (SimulationDataFrame_Q_m["S_r"] == S_r)]
        AverageDataFrame_Q_m = AverageWriteToDataFrame(AverageDataFrame_Q_m, Q_mPlusCounter, MinimumDataFrame_s_m["s_m"].item(), Q_m + counter, s_r, S_r, rowcounter)
        rowcounter += 1

        counter += 1

    MinimumDataFrame_Q_m = AverageDataFrame_Q_m.ix[AverageDataFrame_Q_m["Average Result"].idxmin()]

####################################################################################################################################################################################
    # Change s_r (original s_m and Q_m)
    counter = 0
    rowcounter = 1

    SimulationDataFrame_s_r = pd.DataFrame()

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if s_r - counter >= s_m:
                SimulationDataFrame_s_r = SimWriteToDataFrame(SimulationDataFrame_s_r, s_m, Q_m, s_r - counter, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

            #s_r + counter
            if s_r + counter < S_r:
                SimulationDataFrame_s_r = SimWriteToDataFrame(SimulationDataFrame_s_r, s_m, Q_m, s_r + counter, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

        counter += 1

    # Change s_r (minimum combined s_m and Q_m)
    counter = 0
    rowcounter = 1

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if s_r - counter >= MinimumDataFrame_Q_m["s_m"].item():
                SimulationDataFrame_s_r = SimWriteToDataFrame(SimulationDataFrame_s_r, MinimumDataFrame_Q_m["s_m"].item(), MinimumDataFrame_Q_m["Q_m"].item(), s_r - counter, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

            #s_r + counter
            if s_r + counter < S_r:
                SimulationDataFrame_s_r = SimWriteToDataFrame(SimulationDataFrame_s_r, MinimumDataFrame_Q_m["s_m"].item(), MinimumDataFrame_Q_m["Q_m"].item(), s_r + counter, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

        counter += 1

    # Get Average s_r
    counter = 0
    rowcounter = 1

    AverageDataFrame_s_r = pd.DataFrame()

    while counter <= DeltaNeighborhood:

        if s_r - counter >= MinimumDataFrame_Q_m["s_m"].item():
            s_rMinusCounter = SimulationDataFrame_s_r[(SimulationDataFrame_s_r["s_m"] == s_m) & (SimulationDataFrame_s_r["Q_m"] == Q_m) & (SimulationDataFrame_s_r["s_r"] == s_r - counter) & (SimulationDataFrame_s_r["S_r"] == S_r)]
            AverageDataFrame_s_r = AverageWriteToDataFrame(AverageDataFrame_s_r, s_rMinusCounter, s_m, Q_m, s_r - counter, S_r, rowcounter)
            rowcounter += 1

            s_rMinusCounter = SimulationDataFrame_s_r[(SimulationDataFrame_s_r["s_m"] == MinimumDataFrame_Q_m["s_m"].item()) & (SimulationDataFrame_s_r["Q_m"] == MinimumDataFrame_Q_m["Q_m"].item()) & (SimulationDataFrame_s_r["s_r"] == s_r - counter) & (SimulationDataFrame_s_r["S_r"] == S_r)]
            AverageDataFrame_s_r = AverageWriteToDataFrame(AverageDataFrame_s_r, s_rMinusCounter, MinimumDataFrame_Q_m["s_m"].item(), MinimumDataFrame_Q_m["Q_m"].item(), s_r - counter, S_r, rowcounter)
            rowcounter += 1




        s_rPlusCounter = SimulationDataFrame_s_r[(SimulationDataFrame_s_r["s_m"] == s_m) & (SimulationDataFrame_s_r["Q_m"] == Q_m) & (SimulationDataFrame_s_r["s_r"] == s_r + counter) & (SimulationDataFrame_s_r["S_r"] == S_r)]
        AverageDataFrame_s_r = AverageWriteToDataFrame(AverageDataFrame_s_r, s_rPlusCounter, s_m, Q_m, s_r + counter, S_r, rowcounter)
        rowcounter += 1

        s_rPlusCounter = SimulationDataFrame_s_r[(SimulationDataFrame_s_r["s_m"] == MinimumDataFrame_Q_m["s_m"].item()) & (SimulationDataFrame_s_r["Q_m"] == MinimumDataFrame_Q_m["Q_m"].item()) & (SimulationDataFrame_s_r["s_r"] == s_r + counter) & (SimulationDataFrame_s_r["S_r"] == S_r)]
        AverageDataFrame_s_r = AverageWriteToDataFrame(AverageDataFrame_s_r, s_rPlusCounter, MinimumDataFrame_Q_m["s_m"].item(), MinimumDataFrame_Q_m["Q_m"].item(), s_r + counter, S_r, rowcounter)
        rowcounter += 1

        counter += 1

    MinimumDataFrame_s_r = AverageDataFrame_s_r.ix[AverageDataFrame_s_r["Average Result"].idxmin()]

####################################################################################################################################################################################
    # Change S_r (original s_m, Q_m, s_r)
    counter = 0
    rowcounter = 1

    SimulationDataFrame_S_r = pd.DataFrame()

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if S_r - counter > s_r:
                SimulationDataFrame_S_r = SimWriteToDataFrame(SimulationDataFrame_S_r, s_m, Q_m, s_r, S_r - counter, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

            #S_r + counter

            SimulationDataFrame_S_r = SimWriteToDataFrame(SimulationDataFrame_S_r, s_m, Q_m, s_r, S_r + counter, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
            rowcounter += 1

        counter += 1

    # Change S_r (minimum combined s_m, Q_m, S_r)
    counter = 0
    rowcounter = 1

    while counter <= DeltaNeighborhood:
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            if S_r - counter > MinimumDataFrame_s_r["s_r"].item():
                SimulationDataFrame_S_r = SimWriteToDataFrame(SimulationDataFrame_S_r, MinimumDataFrame_s_r["s_m"].item(), MinimumDataFrame_s_r["Q_m"].item(), MinimumDataFrame_s_r["s_r"].item(), S_r - counter, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
                rowcounter += 1

            #S_r + counter
            SimulationDataFrame_S_r = SimWriteToDataFrame(SimulationDataFrame_S_r, MinimumDataFrame_s_r["s_m"].item(), MinimumDataFrame_s_r["Q_m"].item(), MinimumDataFrame_s_r["s_r"].item(), S_r + counter, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate,L_m, L_r, seed, rowcounter, 100,500)
            rowcounter += 1

        counter += 1

    # Get Average S_r
    counter = 0
    rowcounter = 1

    AverageDataFrame_S_r = pd.DataFrame()

    while counter <= DeltaNeighborhood:

        if S_r - counter > MinimumDataFrame_s_r["s_r"].item():
            S_rMinusCounter = SimulationDataFrame_S_r[(SimulationDataFrame_S_r["s_m"] == s_m) & (SimulationDataFrame_S_r["Q_m"] == Q_m) & (SimulationDataFrame_S_r["s_r"] == s_r) & (SimulationDataFrame_S_r["S_r"] == S_r - counter)]
            AverageDataFrame_S_r = AverageWriteToDataFrame(AverageDataFrame_S_r, S_rMinusCounter, s_m, Q_m, s_r, S_r - counter, rowcounter)
            rowcounter += 1

            S_rMinusCounter = SimulationDataFrame_S_r[(SimulationDataFrame_S_r["s_m"] == MinimumDataFrame_s_r["s_m"].item()) & (SimulationDataFrame_S_r["Q_m"] == MinimumDataFrame_s_r["Q_m"].item()) & (SimulationDataFrame_S_r["s_r"] == MinimumDataFrame_s_r["s_r"].item()) & (SimulationDataFrame_S_r["S_r"] == S_r - counter)]
            AverageDataFrame_S_r = AverageWriteToDataFrame(AverageDataFrame_S_r, S_rMinusCounter, MinimumDataFrame_s_r["s_m"].item(), MinimumDataFrame_s_r["Q_m"].item(), MinimumDataFrame_s_r["s_r"].item(), S_r - counter, rowcounter)
            rowcounter += 1




        S_rPlusCounter = SimulationDataFrame_S_r[(SimulationDataFrame_S_r["s_m"] == s_m) & (SimulationDataFrame_S_r["Q_m"] == Q_m) & (SimulationDataFrame_S_r["s_r"] == s_r) & (SimulationDataFrame_S_r["S_r"] == S_r+ counter)]
        AverageDataFrame_S_r = AverageWriteToDataFrame(AverageDataFrame_S_r, S_rPlusCounter, s_m, Q_m, s_r, S_r + counter, rowcounter)
        rowcounter += 1

        S_rPlusCounter = SimulationDataFrame_S_r[(SimulationDataFrame_S_r["s_m"] == MinimumDataFrame_s_r["s_m"].item()) & (SimulationDataFrame_S_r["Q_m"] == MinimumDataFrame_s_r["Q_m"].item()) & (SimulationDataFrame_S_r["s_r"] == MinimumDataFrame_s_r["s_r"].item()) & (SimulationDataFrame_S_r["S_r"] == S_r + counter)]
        AverageDataFrame_S_r = AverageWriteToDataFrame(AverageDataFrame_S_r, S_rPlusCounter, MinimumDataFrame_s_r["s_m"].item(), MinimumDataFrame_s_r["Q_m"].item(), MinimumDataFrame_s_r["s_r"].item(), S_r + counter, rowcounter)
        rowcounter += 1

        counter += 1

    MinimumDataFrame_S_r = AverageDataFrame_S_r.ix[AverageDataFrame_S_r["Average Result"].idxmin()]

    ####################################################################################################################################################################################

    SimulationDataFrame = pd.DataFrame()

    SimulationDataFrame = SimulationDataFrame.append(SimulationDataFrame_s_m)
    SimulationDataFrame = SimulationDataFrame.append(SimulationDataFrame_Q_m)
    SimulationDataFrame = SimulationDataFrame.append(SimulationDataFrame_s_r)
    SimulationDataFrame = SimulationDataFrame.append(SimulationDataFrame_S_r)

    #rearrange dataframe columns
    cols = SimulationDataFrame.columns.tolist()
    cols[1], cols[6] = cols[6], cols[1]
    cols[2], cols[6] = cols[6], cols[2]
    cols[5], cols[3] = cols[3], cols[5]
    cols[6], cols[4] = cols[4], cols[6]
    cols[5], cols[6] = cols[6], cols[5]
    SimulationDataFrame = SimulationDataFrame[cols]

    AverageDataFrame = pd.DataFrame()

    AverageDataFrame = AverageDataFrame.append(AverageDataFrame_s_m)
    AverageDataFrame = AverageDataFrame.append(AverageDataFrame_Q_m)
    AverageDataFrame = AverageDataFrame.append(AverageDataFrame_s_r)
    AverageDataFrame = AverageDataFrame.append(AverageDataFrame_S_r)

    AveCols = AverageDataFrame.columns.tolist()
    AveCols[0], AveCols[1] = AveCols[1], AveCols[0]
    AveCols[1], AveCols[4] =AveCols[4], AveCols[1]
    AveCols[3], AveCols[4] = AveCols[4], AveCols[3]
    AveCols[3], AveCols[5] =AveCols[5], AveCols[3]
    AverageDataFrame = AverageDataFrame[AveCols]

    MinimumDataFrameC = pd.DataFrame()
    MinimumDataFrameC = MinimumDataFrameC.append(MinimumDataFrame_s_m,ignore_index=True)
    MinimumDataFrameC = MinimumDataFrameC.append(MinimumDataFrame_Q_m,ignore_index=True)
    MinimumDataFrameC = MinimumDataFrameC.append(MinimumDataFrame_s_r,ignore_index=True)
    MinimumDataFrameC = MinimumDataFrameC.append(MinimumDataFrame_S_r,ignore_index=True)


    MinimumDataFrame = MinimumDataFrameC.ix[MinimumDataFrameC["Average Result"].idxmin()]
    #print MinimumDataFrame

    #get full simulation of original and the minimum:
    FinalSimulationsOrig = pd.DataFrame()
    FinalSimulationsMin = pd.DataFrame()
    AverageSimOrig = pd.DataFrame()
    AverageSimMin = pd.DataFrame()

    rowcounter = 1
    for index, row in randomnumbers.iterrows():
        seed = row["RandomNumbers"]
        FinalSimulationsOrig = SimWriteToDataFrame(FinalSimulationsOrig, s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate,return_rate, L_m, L_r, seed, rowcounter, TimeToSteadyState, TimeUnits)
        rowcounter += 1
    AverageSimOrig = AverageWriteToDataFrame(AverageSimOrig, FinalSimulationsOrig, s_m, Q_m, s_r, S_r, 1)

    rowcounter = 1
    for index, row in randomnumbers.iterrows():
        seed = row["RandomNumbers"]
        FinalSimulationsMin = SimWriteToDataFrame(FinalSimulationsMin, MinimumDataFrame["s_m"].item(),  MinimumDataFrame["Q_m"].item(),  MinimumDataFrame["s_r"].item(),  MinimumDataFrame["S_r"].item(), K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate,return_rate, L_m, L_r, seed, rowcounter, TimeToSteadyState, TimeUnits)
        rowcounter += 1

    AverageSimMin = AverageWriteToDataFrame(AverageSimMin, FinalSimulationsMin, MinimumDataFrame["s_m"].item(),  MinimumDataFrame["Q_m"].item(),  MinimumDataFrame["s_r"].item(),  MinimumDataFrame["S_r"].item(), 1)






    return SimulationDataFrame, AverageSimOrig, AverageSimMin

