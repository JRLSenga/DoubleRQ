
import numpy as np
from ModifiedFedergruenandZheng2 import ModifiedFedergruenAndZheng

class Simulation:
    def __init__(self, r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, b, demand_rate, return_rate, L_p, L_r, TimeToSteadyState): #initializer
        self.serviceable_inventory = r_p + Q_p
        self.inventory_position = r_p + Q_p
        self.returns_inventory = 0.0

        #initialization of clock, and interarrival times
        self.clock = 0.0
        self.demand_rate = demand_rate
        self.return_rate = return_rate
        self.t_demand = self.generate_interarrival_demand()
        self.t_returns = self.generate_interarrival_returns()
        self.t_mfgLeadTime = [float("inf")]
        self.t_remfgLeadTime = [float("inf")]

        #cost initialization
        self.cost_total = 0.0
        self.cost_orders_mfg = 0.0
        self.cost_orders_remfg = 0.0
        self.cost_holding_serviceable = 0.0
        self.cost_holding_returns = 0.0
        self.cost_backorder = 0.0

        #parameter setups - taken from the inputs of the simulation
        self.r_p = r_p
        self.Q_p = Q_p
        self.r_r = r_r
        self.Q_r = Q_r
        self.c_p = c_p
        self.c_r = c_r
        self.K_p = K_p
        self.K_r = K_r
        self.h_p = h_p
        self.h_r = h_r
        self.b = b

        self.L_p = L_p
        self.L_r = L_r

        self.TimeToSteadyState = TimeToSteadyState

        # print self.cost_total
    #advances time in the simulation
    def advance_time(self):
        # print "t_event: ", self.t_demand, self.t_returns, self.t_mfgLeadTime[0], self.t_remfgLeadTime[0]
        t_event = min(self.t_demand, self.t_returns, self.t_remfgLeadTime[0], self.t_mfgLeadTime[0])


        if self.clock >= self.TimeToSteadyState:
            self.cost_holding_serviceable += (self.h_p * max(self.serviceable_inventory,0))*(t_event - self.clock)
            self.cost_holding_returns += (self.h_r * self.returns_inventory)*(t_event - self.clock)
            self.cost_backorder -= (self.b * min(self.serviceable_inventory, 0))*(t_event - self.clock)


            self.cost_total += ((self.h_p * max(self.serviceable_inventory,0)) + (self.h_r * self.returns_inventory) - (self.b * min(self.serviceable_inventory, 0)))*(t_event - self.clock)
            # print self.t_demand, self.t_returns, self.t_remfgLeadTime, self.t_mfgLeadTime
            # print "costs: ", self.cost_total, self.cost_holding_serviceable, self.cost_holding_returns, self.cost_backorder, self.cost_orders_remfg, self.cost_orders_mfg, "serviceable: ", self.serviceable_inventory, "returns: ", self.returns_inventory, "position: ", self.inventory_position, t_event - self.clock, t_event, self.clock, self.cost_total / (
            #             self.clock + 1)
            # print " "
            # self.lol = raw_input()
            # if self.lol == "Y":
            #     pass

        self.clock = t_event

        if t_event == self.t_demand:
            self.handle_demand_event()
        elif t_event == self.t_returns:
            self.handle_returns_event()
        elif t_event == self.t_remfgLeadTime[0]:
            self.handle_remfg_event()
        elif t_event == self.t_mfgLeadTime[0]:
            self.handle_mfg_event()

    #Demand Arrivals
    def handle_demand_event(self):
        self.serviceable_inventory -= 1.0
        self.inventory_position -= 1.0

        #Order if at cutoff
        if self.inventory_position == self.r_r:
            if self.returns_inventory >= self.Q_r:
                # print "remfg"
                self.inventory_position += self.Q_r
                self.returns_inventory -= self.Q_r
                if self.t_remfgLeadTime[0] == float("inf"):
                    self.t_remfgLeadTime.pop(0)
                    self.t_remfgLeadTime.append(self.clock + self.L_r)
                else:
                    self.t_remfgLeadTime.append(self.clock + self.L_r)

                if self.clock >= self.TimeToSteadyState:
                    self.cost_orders_remfg += self.K_r + self.c_r * self.Q_r
                    self.cost_total += self.K_r + self.c_r * self.Q_r

        if self.inventory_position == self.r_p:
            # print "mfg"
            self.inventory_position += self.Q_p
            if self.t_mfgLeadTime[0] == float("inf"):
                self.t_mfgLeadTime.pop(0)
                self.t_mfgLeadTime.append(self.clock + self.L_p)
            else:
                self.t_mfgLeadTime.append(self.clock + self.L_p)

            #update mfg ordering cost
            if self.clock >= self.TimeToSteadyState:
                self.cost_orders_mfg += self.K_p + self.c_p * self.Q_p
                self.cost_total += self.K_p + self.c_p * self.Q_p

        #update next demand arrival
        self.t_demand = self.clock + self.generate_interarrival_demand()

    #Return Arrivals
    def handle_returns_event(self):
        self.returns_inventory += 1

        #update next return
        self.t_returns = self.clock + self.generate_interarrival_returns()

    #Handle
    def handle_mfg_event(self):
        self.serviceable_inventory += self.Q_p
        if len(self.t_mfgLeadTime) == 1:
            self.t_mfgLeadTime.pop(0)
            self.t_mfgLeadTime.append(float("inf"))
        else:
            self.t_mfgLeadTime.pop(0)

    def handle_remfg_event(self):
        self.serviceable_inventory += self.Q_r
        if len(self.t_remfgLeadTime) == 1:
            self.t_remfgLeadTime.pop(0)
            self.t_remfgLeadTime.append(float("inf"))
        else:
            self.t_remfgLeadTime.pop(0)

    def generate_interarrival_demand(self):
        return np.random.exponential(1.0/self.demand_rate)#uses scale, not the lambda. so if Lambda = 5, we need to put 1/5

    def generate_interarrival_returns(self):
        return np.random.exponential(1.0/self.return_rate)




#def __init__(self, r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, b, demand_rate, return_rate, L_p, L_r): #initializer
#note that  the assumption is that r_r >= r_p

def SimulateDoubleRQ(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits):

    np.random.seed(seed)



    Trial = Simulation(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, TimeToSteadyState)
    IP = []
    SI = []
    RI = []
    Time = []
    counter = 0

    while Trial.clock <= TimeUnits:
        Trial.advance_time()
        #for graphing
        # counter += 1
        # Time.append(counter)
        # IP.append(Trial.inventory_position)
        # SI.append(Trial.serviceable_inventory)
        # RI.append(Trial.returns_inventory)


        #print Trial.t_mfgLeadTime
        #print Trial.t_remfgLeadTime
        # print Trial.inventory_position
        # print Trial.serviceable_inventory
        # print Trial.returns_inventory
        # print " "

    #remember to remove the holding cost pls
    # print("YOO")
    return float(Trial.cost_total)/(float(Trial.clock) - float(TimeToSteadyState))#, float(Trial.cost_holding_returns)/(float(Trial.clock) - float(TimeToSteadyState)), float(Trial.cost_orders_remfg)/(float(Trial.clock) - float(TimeToSteadyState))

# Lambda_d = 1.0
# Lambda_r = 0.1
# L_p = 4
# L_r = L_p
# h_r = 0.5
# h_p = 1.0
# c_p = 0
# c_r = 0
# K_r = 10.0
# K_p = 100.0
# p = 10.0
#
# print(SimulateDoubleRQ(2.0, 15.0, 5.0, 15.0, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, 17, 1000, 15000))
    #this would save the image of the graph
    # plt.plot(Time, IP, label = "IP")
    # plt.plot(Time, SI, label = "SI")
    # plt.plot(Time, RI, label = "RI")
    # plt.legend(loc = "upper left")
    # plt.savefig("trial.png")
        #print s.cost_total
        #print s.cost_holding_serviceable
        #print s.cost_holding_returns
        #print s.cost_backorder
        #print " "
# c_r = 0.0
# c_p = 0.0
# K_r = 350.0
# K_p = 500.0
# h_p = 50.0
# h_r = 20.0
# p = 60.0
# Lambda_d = 0.5
# Lambda_r = 0.4
# L_p = 2.0
# TimeToSteadyState = 100.0
# TimeUnits = 1000.0
# L_r = 2.0
# seed = 1150

# c_r = 0.0
# c_p = 0.0
# K_r = 30.0
# K_p = 100.0
# h_p = 50.0
# h_r = 0.5
# p = 60.0
# Lambda_d = 1.0
# Lambda_r = 0.75
# L_p = 5.0
# TimeToSteadyState = 0.0
# TimeUnits = 100.0
# L_r = 5.0
# seed = 1150

# h_p = 1.0
# h_r = 0.4
# p = 1.2
#
# #Lambda_r < Lambda_d
# Lambda_d = 0.5
# Lambda_r = 0.25
#
# L_p = 2.0
# L_r = 2.0
# c_p = 0.0
# c_r = 0.0
#
# #K_r < K_p
# K_r = 10.0
# K_p = 100.0
# TimeToSteadyState = 1000.0
# TimeUnits = 10000.0
# seed = 115045
# #
# r_r, r_p, Q_r, Q_p, Z = ModifiedFedergruenAndZheng(Lambda_d, Lambda_r, K_r, K_p, c_p, c_r, h_p, h_r, p, L_p, L_r)
# print r_r, Q_r, r_p, Q_p, Z
# # S, HoldingCostReturns, Remanufacturing = SimulateDoubleRQ(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits)
# # print HoldingCostReturns, Remanufacturing, S, Z, (S-Z)/Z
# # print SimulateDoubleRQ(-2, 6, 0, 1, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits)
# print SimulateDoubleRQ(r_p, Q_p, r_r, Q_r, K_p, K_r, c_p, c_r, h_p, h_r, p, Lambda_d, Lambda_r, L_p, L_r, seed, TimeToSteadyState, TimeUnits)