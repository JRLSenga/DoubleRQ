

import math




def CDFLeadTimeDemand(y, Lambda, L):
    DeltaG = 0.0
    j = 0.0
    evalue = math.exp(-Lambda * L)
    while j <= y - 1:
        # DeltaG = 0.0
        i = 0.0
        while i <= j:
            if i == 0.0:
                DeltaG = DeltaG + (evalue * ((Lambda * L) ** i)) / math.factorial(i)
            else:
                logdenominator = 0.0
                k = 1.0
                ##log factorial computation
                while k <= i:
                    logdenominator = logdenominator + math.log(k)
                    k += 1.0
                lognumerator = math.log(evalue) + (math.log(Lambda * L) * i)
                logdivision = lognumerator - logdenominator
                converted = math.exp(logdivision)
                DeltaG = DeltaG + converted
            i += 1.0
        j += 1.0
    return DeltaG



def VdlAndTeunter_Estimate(h_s, h_r, b, Lambda_d,Lambda_r, L, K_m, K_r):

    Q_m = math.sqrt((2*K_m*(Lambda_d))/((Lambda_r/(Lambda_d*h_r))+(1-Lambda_r/Lambda_d)*h_s))
    Q_r = math.sqrt((2*K_r*Lambda_r)/(h_r*Lambda_r/Lambda_d + h_s))

    threshold_s_m = 1 - (h_s*Q_m)/(b*Lambda_d)
    threshold_s_r = 1 - (h_s*Q_r)/(b*Lambda_d)

    s_m = 0.0
    s_r = 0.0
    while CDFLeadTimeDemand(s_m, Lambda_d, L) <= threshold_s_m:
        s_m += 1

    while CDFLeadTimeDemand(s_r, Lambda_d, L) <= threshold_s_r:
        s_r +=1

    if math.ceil(s_m) <= math.ceil(s_r) and math.ceil(s_r) <= math.ceil(s_m) + math.ceil(Q_m):
        Type = "General"
    else:
        Type = "Simple"
        same_threshold = 1- h_s/((b*((Lambda_d - Lambda_r)/Q_m)+ Lambda_r/Q_r))
        s_m = 0.0
        s_r = 0.0
        while CDFLeadTimeDemand(s_m, Lambda_d, L) <= same_threshold:
            s_m += 1
            s_r += 1



    return math.ceil(s_m), math.ceil(Q_m), math.ceil(s_r), math.ceil(s_r + Q_r), Type
