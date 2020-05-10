import numpy as np
import scipy.stats as ss

S0 = float(input("S0: "))
r = float(input("r: "))
q = float(input("q: "))
sigma = float(input("sigma: "))
T = float(input("T: "))
K1 = float(input("K1: "))
K2 = float(input("K2: "))
K3 = float(input("K3: "))
K4 = float(input("K4: "))


def option_price(S0, r, q, sigma, T, K1, K2, K3, K4):  # compute the option price through derived closed-form formula

    dA1 = (np.log(S0 / K1) + (r - q - (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dA2 = (np.log(S0 / K2) + (r - q - (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dA3 = (np.log(S0 / K1) + (r - q + (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dA4 = (np.log(S0 / K2) + (r - q + (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dC1 = (np.log(S0 / K3) + (r - q - (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dC2 = (np.log(S0 / K4) + (r - q - (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dC3 = (np.log(S0 / K3) + (r - q + (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    dC4 = (np.log(S0 / K4) + (r - q + (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))

    payoff1 = S0 * np.exp((r - q) * T) * (ss.norm.cdf(dA3) - ss.norm.cdf(dA4)) - \
            K1 * (ss.norm.cdf(dA1) - ss.norm.cdf(dA2))

    payoff2 = (K2 - K1) * (ss.norm.cdf(dA2) - ss.norm.cdf(dC1))

    payoff3 = ((K1 - K2) / (K3 - K4)) * (K4 * (ss.norm.cdf(dC1) - ss.norm.cdf(dC2)) - \
            S0 * np.exp((r - q) * T) * (ss.norm.cdf(dC3) - ss.norm.cdf(dC4)))

    call_price = np.exp(- r * T) * (payoff1 + payoff2 + payoff3)

    return call_price


output = option_price(S0, r, q, sigma, T, K1, K2, K3, K4)
print("Through Martingale pricing method, the price would be: " + str(output))

# Monte Carlo Simulation


def Monte_Carlo(S0, r, q, sigma, T, K1, K2, K3, K4):
    ST_list = []
    option_list = []
    for random in range(10000):
        lnST = np.log(S0) + (r - q - (0.5 * sigma * sigma)) * T + sigma * np.sqrt(T) * ss.norm.ppf(np.random.rand())
        ST = np.exp(lnST)
        ST_list.append(ST)

    for i in range(len(ST_list)):
        if ST_list[i] <= K1:
            option = 0
            option_list.append(option)
        elif K1 < ST_list[i] <= K2:
            option = ST_list[i] - K1
            option_list.append(option)
        elif K2 < ST_list[i] <= K3:
            option = K2 - K1
            option_list.append(option)
        elif K3 < ST_list[i] <= K4:
            option = ((K1 - K2) / (K3 - K4)) * (K4 - ST_list[i])
            option_list.append(option)
        else:
            option = 0
            option_list.append(option)

    summing = 0
    for x in range(len(option_list)):
        summing = summing + option_list[x]

    call = np.exp(-r * T) * (summing / len(option_list))
    return call


output = Monte_Carlo(S0, r, q, sigma, T, K1, K2, K3, K4)
print("Through Monte Carlo Simulation, the call price would be: " + str(output))

# Compute the confidence interval

confidence_list = []
for j in range(20):
    Ans = Monte_Carlo(S0, r, q, sigma, T, K1, K2, K3, K4)
    confidence_list.append(Ans)

sum_for_interval = 0
for a in confidence_list:
    sum_for_interval = sum_for_interval + a

average = sum_for_interval / len(confidence_list)

diff_for_interval = 0
for a in confidence_list:
    diff_for_interval = diff_for_interval + ((a - average) ** 2)

standard_deviation = np.sqrt(diff_for_interval / len(confidence_list))

upper = average + 2 * standard_deviation
lower = average - 2 * standard_deviation
print("95% Confidence Interval would be: " + "[" + str(lower) + " - " + str(upper) + "]")
