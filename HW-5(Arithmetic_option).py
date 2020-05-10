import numpy as np
import scipy.stats as ss
from timeit import default_timer as timer
import copy

constantList = input().split()  # in the order of St K r q sigma T_minus_t M n Savg passing_t - input1
St = float(constantList[0])
K = float(constantList[1])
r = float(constantList[2])
q = float(constantList[3])
sigma = float(constantList[4])
T_minus_t = float(constantList[5])
M = int(constantList[6])  # you approximate S with M parts
n = int(constantList[7])
Savg = float(constantList[8])
passing_t = float(constantList[9])
noS = 5000
noR = 20
pass_n = int((passing_t / T_minus_t) * n)

# Binomial tree model for European and American arithmetic calls
dT = T_minus_t / n
u = np.exp(sigma * np.sqrt(dT))
d = np.exp(-sigma * np.sqrt(dT))
p = (np.exp((r - q) * (dT)) - d) / (u - d)

strictprice_list = []
for i in range(n + 1):
    temp_list = []
    for j in range(i + 1):
        Sij = round(St * (u ** (i - j)) * (d ** j), 5)
        temp_list.append(Sij)
    strictprice_list.append(temp_list)

average_list = []
for i in range(n + 1):
    temp_list = []
    for j in range(i + 1):
        #Amax = (St * (1-u **(i-j+1)) / (1-u) + St * u **(i-j) *d*(1-d **j)/(1-d)) / (i+1)
        #Amin = (St * (1-d**(j+1))/(1-d) + St * (d**j) *u * (1-u**(i-j))/(1-u))/(i+1)
        Amax = (Savg * (pass_n + 1) + St * u * (1 - (u ** (i - j))) / (1 - u) + St * (u ** (i - j)) * d * (1 - d ** j) / (1 - d)) / (pass_n + i + 1)
        Amin = (Savg * (pass_n + 1) + St * d * (1 - (d ** j)) / (1 - d) + St * (d ** j) * u * (1 - u ** (i - j)) / (1 - u)) / (pass_n + i + 1)
        temp_list2 = []
        for k in range(M + 1):
            Aijk = (M - k) / M * Amax + k / M * Amin
            temp_list2.append(Aijk)
        temp_list.append(temp_list2)
    average_list.append(temp_list)
# print(average_list)
call_list_E = [[]]
for j in range(n + 1):
    temp_list2 = []
    for k in range(M + 1):
        bt_price = max((average_list[-1][j][k] - K), 0)
        temp_list2.append(bt_price)
    call_list_E[0].append(temp_list2)
# print(call_list_E)
call_list_A = copy.deepcopy(call_list_E)

for i in range(n, 0, -1):
    temp_list = []
    for j in range(i):
        temp_list2 = []
        for k in range(M + 1):
            #Au = ((i+1)*average_list[i - 1][j][k] + strictprice_list[i][j]) / (i + 2)
            #Ad = ((i+1) * average_list[i - 1][j][k] + strictprice_list[i][j + 1]) / (i + 2)
            Au = ((i + pass_n) * average_list[i - 1][j][k] + strictprice_list[i][j]) / (i + pass_n + 1)
            Ad = ((i + pass_n) * average_list[i - 1][j][k] + strictprice_list[i][j + 1]) / (i + pass_n + 1)
            up_no = -1
            down_no = -1
            if j == 0:
                up_no = 0
                down_no = 0
            elif j == i - 1:
                up_no = M
                down_no = M
            else:
                for z in range(M):
                    # if (average_list[i][j][z] - Au) * (average_list[i][j][z + 1] - Au) < 0:
                    #     up_no = z
                    #     break
                    if  average_list[i][j][z] > Au:
                        up_no = up_no + 1
                    else:
                        break     
                for x in range(M):
                    # if (average_list[i][j + 1][x] - Ad) * (average_list[i][j + 1][x + 1] - Ad) < 0:
                    #     down_no = x
                    #     break
                    if  average_list[i][j+1][x] > Ad:
                        down_no = down_no + 1
                    else:
                        break
            if (up_no == 0 and down_no == 0) or (up_no == M and down_no == M):
                Cu = call_list_E[n - i][j][up_no]
                Cd = call_list_E[n - i][j + 1][down_no]
            else:
                Wu = (average_list[i][j][up_no] - Au) / (average_list[i][j][up_no] - average_list[i][j][up_no + 1])
                Cu = Wu * call_list_E[n - i][j][up_no + 1] + (1 - Wu) * call_list_E[n - i][j][up_no]
                Wd = (average_list[i][j + 1][down_no] - Ad) / (average_list[i][j + 1][down_no] - average_list[i][j + 1][down_no + 1])
                Cd = Wd * call_list_E[n - i][j + 1][down_no + 1] + (1 - Wd) * call_list_E[n - i][j + 1][down_no]
            bt_price = np.exp(-r * dT) * (p * Cu + (1 - p) * Cd)
            temp_list2.append(bt_price)
        temp_list.append(temp_list2)
    call_list_E.append(temp_list)
# print(call_list_E[n-2])
print('Binomial tree European arithmetic average call:', call_list_E[n][0][0])

for i in range(n, 0, -1):
    temp_list = []
    for j in range(i):
        temp_list2 = []
        for k in range(M + 1):
            Au = ((i + pass_n) * average_list[i - 1][j][k] + strictprice_list[i][j]) / (i + pass_n + 1)
            Ad = ((i + pass_n) * average_list[i - 1][j][k] + strictprice_list[i][j + 1]) / (i + pass_n + 1)
            up_no = -1
            down_no = -1
            if j == 0:
                up_no = 0
                down_no = 0
            elif j == i - 1:
                up_no = M
                down_no = M
            else:
                for z in range(M):
                    # if (average_list[i][j][z] - Au) * (average_list[i][j][z + 1] - Au) < 0:
                    #     up_no = z
                    #     break
                    if  average_list[i][j][z] > Au:
                        up_no = up_no + 1
                    else:
                        break     
                for x in range(M):
                    # if (average_list[i][j + 1][x] - Ad) * (average_list[i][j + 1][x + 1] - Ad) < 0:
                    #     down_no = x
                    #     break
                    if  average_list[i][j+1][x] > Ad:
                        down_no = down_no + 1
                    else:
                        break
            if up_no == down_no:
                Cu = call_list_A[n - i][j][up_no]
                Cd = call_list_A[n - i][j + 1][down_no]
            else:
                Wu = (average_list[i][j][up_no] - Au) / (average_list[i][j][up_no] - average_list[i][j][up_no + 1])
                Cu = Wu * call_list_A[n - i][j][up_no + 1] + (1 - Wu) * call_list_A[n - i][j][up_no]
                Wd = (average_list[i][j + 1][down_no] - Ad) / (average_list[i][j + 1][down_no] - average_list[i][j + 1][down_no + 1])
                Cd = Wd * call_list_A[n - i][j + 1][down_no + 1] + (1 - Wd) * call_list_A[n - i][j + 1][down_no]
            bt_price = max(np.exp(-r * dT) * (p * Cu + (1 - p) * Cd), average_list[i - 1][j][k] - K)
            temp_list2.append(bt_price)
        temp_list.append(temp_list2)
    call_list_A.append(temp_list)
# print(call_list_A[n-2])
print('Binomial tree American arithmetic average call:', call_list_A[n][0][0])


# Monte Carlo Simulation to price European arithmetic average calls
def Monte_Carlo(St, K, r, q, sigma, T_minus_t, n, Savg, pass_n):

    dT = T_minus_t / n
    option_sum = 0
    for i in range(noS):
        summing = 0
        St_dep = St
        for j in range(n):
            St_dep = np.exp(np.log(St_dep) + (r - q - (0.5 * sigma * sigma)) * dT + sigma * np.sqrt(dT) * np.random.normal())
            summing = summing + St_dep
        call_MC = max((((Savg * (pass_n + 1)) + summing) / (pass_n + n + 1)) - K, 0)  
        option_sum = option_sum + call_MC
    
    MC_price = np.exp(-r * (T_minus_t)) * option_sum / noS

    return MC_price

MC_price = Monte_Carlo(St, K, r, q, sigma, T_minus_t, n, Savg, pass_n)
print('Monte Carlo European arithmetic average call:', str(MC_price))

# Compute the confidence interval
confidence_list = [MC_price]
start = timer()
for j in range(noR - 1):
    Ans = Monte_Carlo(St, K, r, q, sigma, T_minus_t, n, Savg, pass_n)
    confidence_list.append(Ans)

average = np.mean(confidence_list)
standard_deviation = np.std(confidence_list)

upper = average + 2 * standard_deviation
lower = average - 2 * standard_deviation

duration = timer() - start
print(duration, 'seconds')

print("95% Confidence Interval would be: " + "[" + str(lower) + " - " + str(upper) + "]")