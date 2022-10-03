
import math

def G(y, Lambda,c_p, c_r, h_p, h_r, p, L):
    Gy = 0.0
    DeltaG = 0.0
    j = 0.0
    evalue = math.exp(-Lambda*L)
    while j <= y-1:
        #DeltaG = 0.0
        i = 0.0
        while i <= j:
            # print (i, y-1)
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
    Gy = (h_p+p)*DeltaG + p*(Lambda*L - y)
    return Gy

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



def GPlus1(PreviousG, y, Lambda, c_p, c_r, h_p, h_r, p, L):
    evalue = math.exp(-Lambda*L)
    i = 0.0
    DeltaG = 0.0
    lognumerator = 0.0
    logdenominator = 0.0
    converted = 0.0
    logdivision = 0.0
    while i <= y - 1:
        #print i
        if i == 0.0:
            DeltaG = DeltaG + (evalue * ((Lambda * L) ** i)) / math.factorial(i)
        else:
            logdenominator = 0.0
            j = 1
            ##log factorial computation
            while j <= i:
                logdenominator = logdenominator + math.log(j)
                j += 1
            lognumerator = math.log(evalue) + (math.log(Lambda * L) * i)
            logdivision = lognumerator - logdenominator
            converted = math.exp(logdivision)
            DeltaG = DeltaG + converted
        i += 1
    DeltaG = (h_p+p)*DeltaG - p
    DeltaG = DeltaG + PreviousG
    return DeltaG

def E_IL_r(Q_r, Q_p, Lambda_d, Lambda_r):
    Approx = (Lambda_r)/(Lambda_d - Lambda_r)*Q_p + (Q_r - 1.0)/2.0
    return Approx

def RQ(Lambda_d, Lambda_r, K, c_p, c_r, h_p, h_r, p, L_p, L_r, Type, gamma):
    #find minimum of G(y)
    min_G = -100
    min_Gplus1 = min_G + 1
    while G(min_G, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) > G(min_Gplus1, Lambda_d, c_p, c_r, h_p, h_r, p,L_p):
        min_G = min_Gplus1
        min_Gplus1 = min_G + 1

    if Type == "remfg":
        Q = 1.0
        r = min_G #r+1
        Z = K*Lambda_d + G(r, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) + gamma*2
        Storage = Z
        GStorage = G(r + Q+1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)

        while Z > min(G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p), GStorage) + h_r*(Lambda_d / (2*Lambda_r))*(2*Q + 1.0):#h_r * ((Rho / (1.0 - 2.0 * Rho) + 1.0 / 2) * (2.0 * (Q) + 1.0) - 0.5):#h_r * E_IL_r(Q+1, Lambda_d, Lambda_r):

            if G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) < GStorage:
                Storage = (Storage - gamma*Q)*Q + G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
                r = r - 1
            else:
                Storage = (Storage - gamma*Q)*Q + GStorage
                GStorage = GPlus1(GStorage, r + Q + 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
            Q = Q + 1
            Z = Storage / Q + gamma*Q
            Storage = Storage/Q + gamma*Q

            print(Z, G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) + h_r*(Lambda_d / (2*Lambda_r))*(Q + 1.0), GStorage + h_r*(Lambda_d / (2*Lambda_r))*(2*Q + 1.0), r - 1, Q)

        return (r - 1, Q, Z)


    if Type == "mfg":
        Q = 1.0
        r = min_G
        Z = K*Lambda_d + G(r, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) + gamma*2
        Storage = Z
        GStorage = G(r + Q+1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)

        while Z > min(G(r-1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p), GStorage) + h_r*((Lambda_r*Lambda_d)/((Lambda_d-Lambda_r)**2))*(2*Q+1):

            #print Z, G(r-1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p), GStorage, r, Q

            if G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) < GStorage:
                # if r - 1 < 0:
                #     Storage += G(r,Lambda_d, c_p, c_r, h_p, h_r, p,L_p) - GPlus1(0, r-1,Lambda_d,c_p, c_r, h_p, h_r, p, L_p)
                # else:
                Storage =  (Storage - gamma*Q)*Q + G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
                r -= 1
            else:
                Storage = (Storage - gamma*Q)*Q + GStorage
                GStorage = GPlus1(GStorage, r+Q+1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)

            Q += 1
            Z = Storage / float(Q)+ gamma*Q
            Storage = Storage / Q + gamma * Q
            print (Z, G(r - 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p) +  h_r*((Lambda_r*Lambda_d)/((Lambda_d - Lambda_r)**2))*(Q+1), GStorage + h_r*((Lambda_r*Lambda_d)/((Lambda_d - Lambda_r)**2))*(2*Q+1), r - 1, Q)

        return (r-1, Q, Z)


def SatisfyConstraint(r_r, r_p,Lambda_d, Lambda_r, K_r, K_p,c_p, c_r, h_p, h_r, p,L_p, L_r, gamma_r, gamma_p):

    #r_r < r_p
    Z = 10000000.0
    Z_temp = 0.0
    r_r_final = 0
    r_p_final = 0
    for x_r in range(r_r, r_p + 1):
        #print x, r_r, r_p
        Q_p = 1.0
        Q_r = 1.0
        G_r_Storage = G(x_r + 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
        G_p_Storage = G(x_r + 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
        Z_r_temp = K_r*Lambda_d  + G_r_Storage
        Z_p_temp = K_p*Lambda_d  + G_p_Storage
        Storage_r = Z_r_temp
        Storage_p = Z_p_temp


        #computing for cost - need to factor in the linear term
        while Z_r_temp > G_r_Storage + h_r*(Lambda_d / (2*Lambda_r))*(Q_r + 1.0):
            Storage_r = Storage_r + G_r_Storage #+ gamma_r
            G_r_Storage = GPlus1(G_r_Storage, x_r + Q_r + 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
            Q_r += 1
            Z_r_temp = Storage_r / float(Q_r)

        while Z_p_temp > G_p_Storage + h_r*((Lambda_r*Lambda_d)/((Lambda_d - Lambda_r)**2))*(Q_p+1):
            Storage_p = Storage_p + G_p_Storage #+ gamma_p
            G_p_Storage = GPlus1(G_p_Storage, x_r + Q_p + 1, Lambda_d, c_p, c_r, h_p, h_r, p, L_p)
            Q_p += 1
            Z_p_temp = Storage_p / float(Q_p)

        Z_temp += (((Lambda_d - Lambda_r)/ float(Lambda_d)) * Z_p_temp) + ((Lambda_r /float(Lambda_d)) * Z_r_temp + h_r * E_IL_r(Q_r, Q_p, Lambda_d, Lambda_r))

        #print Z_temp, Z
        if Z_temp <= Z:
            Z = Z_temp
            r_r_final = x_r
            r_p_final = x_r

    return r_r_final, Q_r, r_p_final, Q_p, Z

def ModifiedFedergruenAndZheng(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r):

    #check if the K_r and K_p should be different, and if that's the case should the L_p be used in p and L_r be used in r when computing for Gy.
    r_r, Q_r, Z_remfg = RQ(Lambda_d, Lambda_r, K_r, c_p, c_r, h_p, h_r, p,L_p, L_r, "remfg", h_r*Lambda_d/(2*Lambda_r) )
    print("Z_remfg ", Z_remfg)
    r_p, Q_p, Z_mfg = RQ(Lambda_d, Lambda_r, K_p, c_p, c_r, h_p, h_r, p,L_p, L_r, "mfg",h_r*(Lambda_r*Lambda_d/(Lambda_d - Lambda_r)**2))
    print("Z_mfg ", Z_mfg)
    # Z_remfg = Z_remfg - h_r * E_IL_r(Q_r, Lambda_d, Lambda_r)
    # print "Z_remfg after -h_r ", Z_remfg
    Z = (((Lambda_d - Lambda_r)/ float(Lambda_d)) * Z_mfg) + ((Lambda_r / float(Lambda_d)) * Z_remfg) + h_r * E_IL_r(Q_r,Q_p, Lambda_d, Lambda_r)
    print(Z)

    #Combining
    if r_r >= r_p:
        return r_r, r_p, Q_r, Q_p, Z
    else:
        print("entered Modified: ", r_r, Q_r, r_p, Q_p)
        r_r, Q_r, r_p, Q_p, Z = SatisfyConstraint(r_r, r_p, Lambda_d, Lambda_r, K_r, K_p,c_p, c_r, h_p, h_r, p,L_p, L_r,  h_r*Lambda_d/(2*Lambda_r), h_r*(Lambda_r*Lambda_d/(Lambda_d - Lambda_r)**2))
        return r_r, r_p, Q_r, Q_p, Z
