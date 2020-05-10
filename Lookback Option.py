import scipy.stats as ss
import numpy as np
import copy

constantList = input().split()  # in the order of St r q sigma t T Smax n -  input1
St = float(constantList[0])
r = float(constantList[1])
q = float(constantList[2])
sigma = float(constantList[3])
t = float(constantList[4])
T = float(constantList[5])
Smax = float(constantList[6])
n = int(constantList[7])
noS = 1000
noR = 10

# Binomial tree model for European and American lookback puts
dT = (T - t) / n
u = np.exp(sigma * np.sqrt(dT))
d = np.exp(-sigma * np.sqrt(dT))
p = (np.exp((r - q) * (dT)) - d) / (u - d)

strictprice_list = []
for i in range(n + 1):
    temp_list = []
    for j in range(i + 1):
        Sij = round(St * (u ** (i - j)) * (d ** j), 4)
        temp_list.append(Sij)
    strictprice_list.append(temp_list)

smax_list = []
for i in range(n + 1):
    temp_list2 = []
    if i == 0:
        temp_list2.append(Smax)
        smax_list.append([temp_list2])
        continue
    else:
        for j in range(i + 1):
            temp_list3 = []
            if j == 0:
                for k in smax_list[i - 1][j]:
                    if k >= strictprice_list[i][j]:
                        if k not in temp_list3:
                            temp_list3.append(k)
                    else:
                        if strictprice_list[i][j] not in temp_list3:
                            temp_list3.append(strictprice_list[i][j])
                temp_list2.append(temp_list3)
            else:
                if j == i:
                    temp_list3 = []
                    for k in smax_list[i -1][j -1]:
                        if k >= strictprice_list[i][j]:
                            if k not in temp_list3:
                                temp_list3.append(k)
                        else:
                            if strictprice_list[i][j] not in temp_list3:
                                temp_list3.append(strictprice_list[i][j])
                    temp_list2.append(temp_list3)
                else:
                    temp_list3 = []
                    for k in smax_list[i - 1][j - 1]:
                        if k >= strictprice_list[i][j]:
                            if k not in temp_list3:
                                temp_list3.append(k)
                        else:
                            if strictprice_list[i][j] not in temp_list3:
                                temp_list3.append(strictprice_list[i][j])
                    for k in smax_list[i - 1][j]:
                        if k >= strictprice_list[i][j]:
                            if k not in temp_list3:
                                temp_list3.append(k)
                        else:
                            if strictprice_list[i][j] not in temp_list3:
                                temp_list3.append(strictprice_list[i][j])
                    temp_list2.append(temp_list3)
        smax_list.append(temp_list2)

put_list_E = [[]]
for j in range(n + 1):
    temp_list2 = []
    for k in smax_list[-1][j]:
        bt_price = max((k - strictprice_list[-1][j]), 0)
        temp_list2.append(bt_price)
    put_list_E[0].append(temp_list2)
# print(put_list_E)
put_list_A = copy.deepcopy(put_list_E)

for i in range(n, 0, -1):
    temp_list2 = []
    for j in range(i):
        temp_list3 = []
        for k in range(len(smax_list[i - 1][j])):
            up_price_no = -1
            down_price_no = -1
            for z in range(len(smax_list[i][j])):
                if smax_list[i - 1][j][k] == smax_list[i][j][z]:
                    up_price_no = z
                    break
            if up_price_no == -1:
                for z in range(len(smax_list[i][j])):
                    if smax_list[i][j][z] == strictprice_list[i][j]:
                        up_price_no = z
                        break
            for z in range(len(smax_list[i][j + 1])):
                if smax_list[i - 1][j][k] == smax_list[i][j + 1][z]:
                    down_price_no = z
                    break   
            up_price = put_list_Ｅ[n - i][j][up_price_no]
            down_price = put_list_Ｅ[n - i][j + 1][down_price_no]
            bt_price = np.exp(-r * (dT)) * (p * up_price + (1 - p) * down_price)
            temp_list3.append(bt_price)
        temp_list2.append(temp_list3)
    put_list_E.append(temp_list2)

print('Binomial tree European lookback put:', put_list_E[n][0][0])

for i in range(n, 0, -1):
    temp_list2 = []
    for j in range(i):
        temp_list3 = []
        for k in range(len(smax_list[i - 1][j])):
            up_price_no = -1
            down_price_no = -1
            for z in range(len(smax_list[i][j])):
                if smax_list[i - 1][j][k] == smax_list[i][j][z]:
                    up_price_no = z
                    break
            if up_price_no == -1:
                for z in range(len(smax_list[i][j])):
                    if smax_list[i][j][z] == strictprice_list[i][j]:
                        up_price_no = z
                        break
            for z in range(len(smax_list[i][j + 1])):
                if smax_list[i - 1][j][k] == smax_list[i][j + 1][z]:
                    down_price_no = z
                    break   
            up_price = put_list_A[n - i][j][up_price_no]
            down_price = put_list_A[n - i][j + 1][down_price_no]
            bt_price = np.exp(-r * (dT)) * (p * up_price + (1 - p) * down_price)
            temp_list3.append(max(bt_price, (smax_list[i - 1][j][k] - strictprice_list[i - 1][j])))
        temp_list2.append(temp_list3)
    put_list_A.append(temp_list2)

print('Binomial tree American lookback put:', put_list_A[n][0][0])

# Implement Cheuk and Vorst to price European and American lookback puts
# u = np.exp(sigma * np.sqrt(dT))
# d = np.exp(-sigma * np.sqrt(dT))
# p = (np.exp((r - q) * (dT)) - d) / (u - d)
# dT = (T - t) / 1000
cheuk = (np.exp((r - q) * dT) * u - 1) / (np.exp((r - q) * dT) * (u - d))

strictprice_list = []
for i in range(n + 1):
    temp_list = []
    for j in range(i + 1):
        Sij = u ** j
        temp_list.append(Sij)
    temp_list.reverse()
    strictprice_list.append(temp_list)
# print(strictprice_list)

put_list_E = []
for i in range(n + 1):
    bt_price = St * max((strictprice_list[-1][i] - 1), 0)
    put_list_E.append(bt_price)
put_list_A = copy.deepcopy(put_list_E)

for i in range(n, 0, -1):
    for j in range(i):
        if j == i - 1:
            put_list_E[j] = np.exp((r - q) * dT) * np.exp(-r * (dT)) * ((1-cheuk) * put_list_E[j] + cheuk * put_list_E[j + 1])
        else:
            put_list_E[j] = np.exp((r - q) * dT) * np.exp(-r * (dT)) * ((1-cheuk) * put_list_E[j] + cheuk * put_list_E[j + 2])

print('Binomial tree European lookback put:', put_list_E[0])

for i in range(n, 0, -1):
    for j in range(i):
        if j == i - 1:
            put_list_A[j] = max((np.exp((r - q) * dT) * np.exp(-r * (dT)) * ((1-cheuk) * put_list_A[j] + cheuk * put_list_A[j + 1])), St * (strictprice_list[i-1][j]-1))     
        else:
            put_list_A[j] = max((np.exp((r - q) * dT) * np.exp(-r * (dT)) * ((1-cheuk) * put_list_A[j] + cheuk * put_list_A[j + 2])), St * (strictprice_list[i-1][j]-1))

print('Binomial tree American lookback put:', put_list_A[0])

# Smax list quick version
# while True:
#     quick_list = []
#     udList = input().split()
#     i_search = int(udList[0]) + int(udList[1])
#     j_search = int(udList[1])

#     for j in range(j_search + 1):
#         max_num = St * u ** (i_search - j_search - j)
#         if max_num >= St:
#             quick_list.append(max_num)

#     if Smax > St:
#         quick_list2 = [Smax]
#         for i in quick_list:
#             if i > Smax:
#                 quick_list2.append(i)
#         print(quick_list2) 
#     else:
#         print(quick_list)
#     print(smax_list[i_search][j_search])


# Monte Carlo Simulation to price European lookback puts
def Monte_Carlo(St, r, q, sigma, t, T, n, Smax):
    dT = (T - t) / n
    ST_list = []
    option_list = []
    for random in range(noS):
        temp = [St]
        for i in range(n):
            lnSt_plus_1 = np.log(temp[i]) + (r - q - (0.5 * sigma * sigma)) * dT + sigma * np.sqrt(dT) * ss.norm.ppf(np.random.rand())
            St_plus_1 = np.exp(lnSt_plus_1)
            temp.append(St_plus_1)
        ST_list.append(temp)
    
    for i in range(noS):
        Smax1 = max(max(ST_list[i]), Smax)
        put_price = max(Smax1 - ST_list[i][-1], 0)
        option_list.append(put_price)
    
    summing = 0
    for i in range(noS):
        summing = summing + option_list[i]
    
    MC_price = np.exp(-r * (T - t)) * (summing / noS)

    return MC_price

MC_price = Monte_Carlo(St, r, q, sigma, t, T, n, Smax)
print('Monte Carlo European lookback put:', str(MC_price))

# Compute the confidence interval
confidence_list = []
for j in range(noR):
    Ans = Monte_Carlo(St, r, q, sigma, t, T, n, Smax)
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
