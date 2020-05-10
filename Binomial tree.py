import scipy.stats as ss
import numpy as np
import math

S0 = float(input("S0: "))
K = float(input("K: "))
r = float(input("r: "))
q = float(input("q: "))
sigma = float(input("sigma: "))
T = float(input("T: "))
n = int(input("Number of Simulations: "))
num_rep = int(input("Number of repetitions: "))
option_type = input("Option Type(please enter \"c\" or \"p\"): ")
E_or_A = input("Early exercise or not(please enter \"E\" or \"A\"): ")


# Black-Scholes formula(only for european option)


def Black_Scholes(S0, r, q, sigma, T, K, option_type):
    d1 = (np.log(S0 / K) + (r - q + (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r - q - (sigma ** 2) / 2) * T) / (sigma * np.sqrt(T))
    if option_type == "c" or option_type == "C":
        BS_price = S0 * np.exp(-q * T) * ss.norm.cdf(d1) - K * np.exp(-r * T) * ss.norm.cdf(d2)
    if option_type == "p" or option_type == "P":
        BS_price = K * np.exp(-r * T) * (1 - ss.norm.cdf(d2)) - S0 * np.exp(-q * T) * (1 - ss.norm.cdf(d1))
    return BS_price


if E_or_A == "e" or E_or_A == "E":
    BS_output1= Black_Scholes(S0, r, q, sigma, T, K, option_type)
    print("Through Black-Scholes formulas, the option price would be: " + str(BS_output1))


# Monte Carlo Simulation(only for European option)


def Monte_Carlo(S0, r, q, sigma, T, K, n, option_type):
    ST_list = []
    option_list = []
    for random in range(n):
        lnST = np.log(S0) + (r - q - (0.5 * sigma * sigma)) * T + sigma * np.sqrt(T) * ss.norm.ppf(np.random.rand())
        ST = np.exp(lnST)
        ST_list.append(ST)

    if option_type == "c" or option_type == "C":
        for i in range(len(ST_list)):
            if ST_list[i] <= K:
                option = 0
                option_list.append(option)
            else:
                option = ST_list[i] - K
                option_list.append(option)
    elif option_type == "p" or option_type == "P":
        for i in range(len(ST_list)):
            if ST_list[i] >= K:
                option = 0
                option_list.append(option)
            else:
                option = K - ST_list[i]
                option_list.append(option)
    summing = 0
    for x in range(len(option_list)):
        summing = summing + option_list[x]

    MC_price = np.exp(-r * T) * (summing / len(option_list))
    return MC_price


if E_or_A == "e" or E_or_A == "E":
    MC_output = Monte_Carlo(S0, r, q, sigma, T, K, n, option_type)
    print("Through Monte Carlo Simulation, the option price would be: " + str(MC_output))

# Compute the confidence interval

if E_or_A == "e" or E_or_A == "E":
    confidence_list = []
    for j in range(num_rep):
        Ans = Monte_Carlo(S0, r, q, sigma, T, K, n, option_type)
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

# CRR binomial tree model - basic requirement with 2 dimensions

delta_t = T / n
u = np.exp(sigma * np.sqrt(delta_t))
d = np.exp(-sigma * np.sqrt(delta_t))
p = (np.exp((r - q) * (delta_t)) - d) / (u - d)
print(u)

strictprice_list = []
for i in range(n + 1):
    temp_list = []
    for j in range(i + 1):
        Sij = S0 * (u ** (i - j)) * (d ** j)
        temp_list.append(Sij)
    strictprice_list.append(temp_list)

option_list = [[]]
for j in range(n + 1):
    if option_type == "c" or option_type == "C":
        Cnj = max((strictprice_list[n][j] - K), 0)
        option_list[0].append(Cnj)
    elif option_type == "p" or option_type == "P":
        Pnj = max((K - strictprice_list[n][j]), 0)
        option_list[0].append(Pnj)

for i in range(n, 0, -1):
    temp_list = []
    for j in range(i):
        if E_or_A == "e" or E_or_A == "E":
            price = np.exp(-r * (delta_t)) * (p * option_list[n - i][j] + (1 - p) * option_list[n - i][j + 1])
        else:
            temp_price = np.exp(-r * (delta_t)) * (
                        p * option_list[n - i][j] + (1 - p) * option_list[n - i][j + 1])
            if option_type == "c" or option_type == "C":
                price = max(strictprice_list[i - 1][j] - K, temp_price)
            elif option_type == "p" or option_type == "P":
                price = max(K - strictprice_list[i - 1][j], temp_price)
        temp_list.append(price)
    option_list.append(temp_list)

print("Through CRR binomial tree model, the option price would be: " + str(option_list[n][0]))

# CRR binomial tree model - basic requirement with one column vector

single_option_list = []
for j in range(n + 1):
    if option_type == "c" or option_type == "C":
        Cnj = max((strictprice_list[n][j] - K), 0)
        single_option_list.append(Cnj)
    elif option_type == "p" or option_type == "P":
        Pnj = max((K - strictprice_list[n][j]), 0)
        single_option_list.append(Pnj)

for i in range(n, 0, -1):
    for j in range(i):
        if E_or_A == "e" or E_or_A == "E":
            single_option_list[j] = np.exp(-r * (delta_t)) * (p * single_option_list[j] + (1 - p) * single_option_list[j + 1])
        else:
            temp_price = np.exp(-r * (delta_t)) * (p * single_option_list[j] + (1 - p) * single_option_list[j + 1])
            if option_type == "c" or option_type == "C":
                single_option_list[j] = max(strictprice_list[i - 1][j] - K, temp_price)
            elif option_type == "p" or option_type == "P":
                single_option_list[j] = max(K - strictprice_list[i - 1][j], temp_price)

print("Through CRR binomial tree model, the option price would be: " + str(single_option_list[0]))


# Combinatorial Method for binomial tree

if E_or_A == "e" or E_or_A == "E":
    com_price = 0
    for j in range(n + 1):
        if j > 0:
            factorial = 1
            down = 1
            for x in range(1, j + 1):
                factorial = factorial * (n - x + 1)
                down = down * x
            factor = math.log(factorial) - math.log(down)
        else:
            factor = math.log(1)
        com_price = com_price + np.exp(factor + (n - j) * math.log(p) + j * math.log(1 - p)) * option_list[0][j]

    com_price = com_price * np.exp(-r * T)
    print("Through Combinatorial Method, the option price would be: " + str(com_price))
