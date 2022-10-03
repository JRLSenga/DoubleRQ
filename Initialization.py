import pandas as pd
from ModifiedFedergruenandZheng2 import ModifiedFedergruenAndZheng
from DoubleRQ import SimulateDoubleRQ
from VanDerLaanandTeunter_Parameter_Approximation import VdlAndTeunter_Estimate
from vanDerLaan_sQsS_Simulation_Continuous import SimulatesQsS

def Main():
    rowcounter = 1
    Outputs_sQsS = pd.DataFrame()
    Outputs_DoubleRQ = pd.DataFrame()
    Outputs = pd.DataFrame()
    Overall_Consolated_Results = []


    DoubleRQSimList = []
    PullPolicySimList_DoubleRQ_Params = []
    PullPolicySimList_VDLTeunter_Params = []

    Headers = []

    Headers.append("ID")
    Headers.append("lambda_d")
    Headers.append("lambda_r")
    Headers.append("h_p")
    Headers.append("h_r")
    Headers.append("p")
    Headers.append("K_p")
    Headers.append("K_r")
    Headers.append("L")
    Headers.append("Double RQ")
    Headers.append("Z")
    Headers.append("Pull with Double RQ")
    Headers.append("Pull with Teunter")
    Headers.append("2rQ - s_m")
    Headers.append("2rQ - Q_m")
    Headers.append("2rQ - s_r")
    Headers.append("2rQ - S_r")
    Headers.append("VT - s_m")
    Headers.append("VT - Q_m")
    Headers.append("VT - s_r")
    Headers.append("VT - S_r")
    Headers.append("2rq vs Z")
    Headers.append("2rq vs pull(2rqparams)")
    Headers.append("VdlTeunter vs pull(2rqparams)")
    Headers.append("Type")


    SimulationParameters = pd.read_csv("SampleParameters.csv") #change input filename here

    #iterate through each of the simulation parameters
    for index,row in SimulationParameters.iterrows():

        #Clearing The List
        del DoubleRQSimList[:]
        del PullPolicySimList_DoubleRQ_Params[:]
        del PullPolicySimList_VDLTeunter_Params[:]

        Results = []


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
        TimeToSteadyState = row["TimeToSteadyState"]
        TimeUnits = row["TimeUnits"]

        h_r_ratio = row["h_r ratio"]
        p_ratio = row["p ratio"]
        K_r_ratio = row["K_r ratio"]
        Lambda_r_ratio = row["lambda_r ratio"]

        print(ID)
        #obtain optimal Double RQ parameters and Simulation Results
        r_r, r_p, Q_r, Q_p, Z = ModifiedFedergruenAndZheng(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r)

        s_mBaseLine = r_p
        Q_mBaseLine = Q_p
        s_rBaseLine = r_r
        S_rBaseLine = r_r + Q_r
        print(s_mBaseLine, Q_mBaseLine, s_rBaseLine, S_rBaseLine)
        s_mVDLTeunter, Q_mVDLTeunter, s_rVDLTeunter, S_rVDLTeunter, Type = VdlAndTeunter_Estimate(h_p, h_r, p, Lambda_d,Lambda_r, L_p, K_p, K_r)
        print(s_mVDLTeunter, Q_mVDLTeunter, s_rVDLTeunter, S_rVDLTeunter, Type)
        # print(s_mBaseLine, Q_mBaseLine, s_mBaseLine, S_rBaseLine)
        # print(s_mVDLTeunter, Q_mVDLTeunter, s_rVDLTeunter,S_rVDLTeunter)



        randomnumbers = pd.read_csv("RandomNumbers.csv")#"/Users/juansenga/Dropbox/Projects/Double R,Q/DoubleRQ_Python_Files/RandomNumbers.csv")
        for index, row in randomnumbers.iterrows():
            seed = row["RandomNumbers"]

            DoubleRQSimList.append(SimulateDoubleRQ(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits))
            PullPolicySimList_DoubleRQ_Params.append(SimulatesQsS(s_mBaseLine, Q_mBaseLine, s_rBaseLine, S_rBaseLine, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed,TimeToSteadyState, TimeUnits))
            PullPolicySimList_VDLTeunter_Params.append(SimulatesQsS(s_mVDLTeunter, Q_mVDLTeunter, s_rVDLTeunter, S_rVDLTeunter, K_p, K_r, c_p, c_r, h_p, h_r, p,Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits))


        #print DoubleRQSimList
        DoubleRQSim_Average = sum(DoubleRQSimList)/float(len(DoubleRQSimList))
        PullPolicyDoubleRQ_Average = sum(PullPolicySimList_DoubleRQ_Params)/float(len(PullPolicySimList_DoubleRQ_Params))
        PullPolicyVDLTeunter_Average = sum(PullPolicySimList_VDLTeunter_Params) / float(len(PullPolicySimList_VDLTeunter_Params))


        Results.append(ID)
        Results.append(Lambda_d)
        Results.append(Lambda_r_ratio)
        Results.append(h_p)
        Results.append(h_r_ratio)
        Results.append(p_ratio)
        Results.append(K_p)
        Results.append(K_r_ratio)
        Results.append(L_p)
        Results.append(DoubleRQSim_Average)
        Results.append(Z)
        Results.append(PullPolicyDoubleRQ_Average)
        Results.append(PullPolicyVDLTeunter_Average)
        Results.append(s_mBaseLine)
        Results.append(Q_mBaseLine)
        Results.append(s_rBaseLine)
        Results.append(S_rBaseLine)
        Results.append(s_mVDLTeunter)
        Results.append(Q_mVDLTeunter)
        Results.append(s_rVDLTeunter)
        Results.append(S_rVDLTeunter)
        Results.append(100*abs(DoubleRQSim_Average-Z)/Z)
        Results.append(100*(DoubleRQSim_Average - PullPolicyDoubleRQ_Average)/PullPolicyDoubleRQ_Average)
        Results.append(100*(PullPolicyVDLTeunter_Average - PullPolicyDoubleRQ_Average)/PullPolicyDoubleRQ_Average)
        Results.append(Type)

        Overall_Consolated_Results.append(Results)
        print("Double RQ vs Z = ", round(100*abs(DoubleRQSim_Average-Z)/Z,2), "DoubleRQ vs Pull-DoubleRQ Param = ", round(100*(DoubleRQSim_Average - PullPolicyDoubleRQ_Average)/PullPolicyDoubleRQ_Average,2), "Pull-Double RQ vs Teunter = ", round(100*(PullPolicyVDLTeunter_Average - PullPolicyDoubleRQ_Average)/PullPolicyDoubleRQ_Average,2))
        Data = pd.DataFrame(Overall_Consolated_Results, columns = Headers)


        Data.to_csv("OutputFile_Sample.csv", encoding = 'utf-8', index = False) # change output file name here


Main()