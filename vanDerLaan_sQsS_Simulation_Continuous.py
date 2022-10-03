
import numpy as np

class Simulation:
    def __init__(self, s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, TimeToSteadyState): #initializer
        self.serviceable_inventory = s_m + Q_m
        self.inventory_position = s_m + Q_m
        self.returns_inventory = 0.0

        #initialization of clock, and interarrival times
        self.clock = 0.0
        self.demand_rate = demand_rate
        self.return_rate = return_rate
        self.t_demand = self.generate_interarrival_demand()
        self.t_returns = self.generate_interarrival_returns()
        self.t_mfgLeadTime = [float("inf")]
        self.t_remfgLeadTime = [float("inf")]
        self.remfg_amount = []

        #cost initialization
        self.cost_total = 0.0
        self.cost_orders_mfg = 0.0
        self.cost_orders_remfg = 0.0
        self.cost_holding_serviceable = 0.0
        self.cost_holding_returns = 0.0
        self.cost_backorder = 0.0

        #parameter setups - taken from the inputs of the simulation
        self.s_m = s_m
        self.Q_m = Q_m
        self.s_r = s_r
        self.S_r = S_r
        self.c_m = c_m
        self.c_r = c_r
        self.K_m = K_m
        self.K_r = K_r
        self.h_m = h_m
        self.h_r = h_r
        self.b = b

        self.L_m = L_m
        self.L_r = L_r

        self.TimeToSteadyState = TimeToSteadyState


    #advances time in the simulation
    def advance_time(self):
        t_event = min(self.t_demand, self.t_returns, self.t_remfgLeadTime[0], self.t_mfgLeadTime[0])
        if self.clock >= self.TimeToSteadyState:
            self.cost_holding_serviceable += (self.h_m * max(self.serviceable_inventory,0))*(t_event - self.clock)
            self.cost_holding_returns += (self.h_r * self.returns_inventory)*(t_event - self.clock)
            self.cost_backorder -= (self.b * min(self.serviceable_inventory, 0))*(t_event - self.clock)

            self.cost_total += ((self.h_m * max(self.serviceable_inventory,0)) + (self.h_r * self.returns_inventory) - (self.b * min(self.serviceable_inventory, 0)))*(t_event - self.clock)

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
        self.serviceable_inventory -= 1
        self.inventory_position -= 1

        #Order if at cutoff
        # if self.s_m <> self.s_r:

        if self.inventory_position <= self.s_r and self.returns_inventory >= self.S_r - self.inventory_position:
            self.order_up_to_mechanism()

        elif self.inventory_position == self.s_m and self.returns_inventory < self.S_r - self.inventory_position:
            self.inventory_position += self.Q_m
            if self.t_mfgLeadTime[0] == float("inf"):
                self.t_mfgLeadTime.pop(0)
                self.t_mfgLeadTime.append(self.clock + self.L_m)
            else:
                self.t_mfgLeadTime.append(self.clock + self.L_m)


            #update mfg ordering cost
            if self.clock >= self.TimeToSteadyState:
                self.cost_orders_mfg += self.K_m + self.c_m * self.Q_m
                self.cost_total += self.K_m + self.c_m * self.Q_m


        #update next demand arrival
        self.t_demand = self.clock + self.generate_interarrival_demand()

    #Return Arrivals
    def handle_returns_event(self):
        self.returns_inventory += 1

        #Order if there is enough returns
        if self.inventory_position <= self.s_r:
            self.order_up_to_mechanism()

        #update next return
        self.t_returns = self.clock + self.generate_interarrival_returns()

    #Handle
    def handle_mfg_event(self):
        #print "mfg"
        self.serviceable_inventory += self.Q_m
        if len(self.t_mfgLeadTime) == 1:
            self.t_mfgLeadTime.pop(0)
            self.t_mfgLeadTime.append(float("inf"))
        else:
            self.t_mfgLeadTime.pop(0)


    def handle_remfg_event(self):
        self.serviceable_inventory += self.remfg_amount[0]
        self.remfg_amount.pop(0)
        if len(self.t_remfgLeadTime) == 1:
            self.t_remfgLeadTime.pop(0)
            self.t_remfgLeadTime.append(float("inf"))
        else:
            self.t_remfgLeadTime.pop(0)

    def order_up_to_mechanism(self):
        if self.inventory_position <= self.s_r:
            if self.returns_inventory >= self.S_r - self.inventory_position:
                self.returns_inventory = self.returns_inventory - (self.S_r - self.inventory_position)
                self.remfg_amount.append(self.S_r - self.inventory_position)
                self.inventory_position = self.S_r
                if self.t_remfgLeadTime[0] == float("inf"):
                    self.t_remfgLeadTime.pop(0)
                    self.t_remfgLeadTime.append(self.clock + self.L_r)
                else:
                    self.t_remfgLeadTime.append(self.clock + self.L_r)

                # update remfg ordering cost
                if self.clock >= self.TimeToSteadyState:
                    self.cost_orders_remfg += self.K_r + self.c_r * (self.S_r - self.inventory_position)
                    self.cost_total += self.K_r + self.c_r * (self.S_r - self.inventory_position)

    def generate_interarrival_demand(self):
        return np.random.exponential(1.0 / self.demand_rate)#uses scale, not the lambda. so if Lambda = 5, we need to put 1/5

    def generate_interarrival_returns(self):
        return np.random.exponential(1.0 / self.return_rate)



#def __init__(self, s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r): #initializer
#note that  the assumption is that s_r >= s_m
def SimulatesQsS(s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, seed, TimeToSteadyState, TimeUnits):
    np.random.seed(seed)
    Trial = Simulation(s_m, Q_m, s_r, S_r, K_m, K_r, c_m, c_r, h_m, h_r, b, demand_rate, return_rate, L_m, L_r, TimeToSteadyState)
    IP = []
    SI = []
    RI = []
    Time = []
    counter = 0

    while Trial.clock <= TimeUnits:
        Trial.advance_time()
        #for graphing
        counter += 1
        Time.append(counter)
        IP.append(Trial.inventory_position)
        SI.append(Trial.serviceable_inventory)
        RI.append(Trial.returns_inventory)

    return Trial.cost_total / (Trial.clock - Trial.TimeToSteadyState)

