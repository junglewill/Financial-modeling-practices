import scipy.stats as ss
import numpy as np
from numpy.linalg import inv
import copy

# the Cholesky decomposition method to price rainbow option
constantList = input().split()  # in the order of K r T n  -  input1
K = float(constantList[0])
r = float(constantList[1])
T = float(constantList[2])
n = int(constantList[3])
noS = 10000  # number of simulations
noR = 20  # number of repetitions
S0List = input().split()  # input2
qList = input().split()  # input3
sigmaList = input().split()  # input4
lo_ijList = []
for i in range(n):  # input5 ~ n+5
    temp = input().split()
    lo_ijList.append(temp)

for i in range(n):
    S0List[i] = float(S0List[i])
    qList[i] = float(qList[i])
    sigmaList[i] = float(sigmaList[i])
    for j in range(n):
        lo_ijList[i][j] = float(lo_ijList[i][j])

C = []  # get the covariance matrix
for i in range(n):
    temp2 = []
    for j in range(n):
        cov = sigmaList[i] * sigmaList[j] * lo_ijList[i][j] * T
        temp2.append(cov)
    C.append(temp2)

A = []
for i in range(n):  # apply the algorithm to compute A
    zeroList = []
    for j in range(n):
        zeroList.append(0)
    A.append(zeroList)

for i in range(n):
    for j in range(n):
        if i == 0:
            if j == 0:
                A[i][j] = np.sqrt(C[i][j])
            else:
                A[i][j] = C[i][j] / A[i][i]
        else:
            if i == j:
                tempsum = 0
                for k in range(i):
                    tempsum = tempsum + (A[k][i] ** 2)
                A[i][j] = np.sqrt(C[i][j] - tempsum)
            elif j > i:
                tempsum = 0
                for k in range(i):
                    tempsum = tempsum + (A[k][i] * A[k][j])
                A[i][j] = (C[i][j] - tempsum) / A[i][i]

def Cholesky(K, r, T, n, noS, S0List, qList, sigmaList):

    ST_List = []
    for i in range(noS):
        zeroList = []
        for j in range(n):
            zeroList.append(0)
        ST_List.append(zeroList)

    antitheticList1 = []
    for i in range(noS):
        temp2 = []
        for j in range(n):
            random = ss.norm.ppf(np.random.rand())
            temp2.append(random)
        if i < (noS // 2):
            antitheticList1.append(temp2)
        for k in range(n):
            rSum = 0
            for x in range(n):
                rSum = rSum + (temp2[x] * A[x][k])
            lnST = np.log(S0List[k]) + (r - qList[k] - (sigmaList[k] * sigmaList[k] * 0.5)) * T + rSum
            ST_List[i][k] = np.exp(lnST)

    sum_of_payoff = 0
    for i in range(noS):
        payoff = max(max(ST_List[i]) - K, 0)
        sum_of_payoff = sum_of_payoff + payoff
    Cho_call = np.exp(-r * T) * (sum_of_payoff / noS)

    return Cho_call, antitheticList1


call_price, list_for_mm = Cholesky(K, r, T, n, noS, S0List, qList, sigmaList)
print("Through Cholesky decomposition, the rainbow option price would be: " + str(call_price))


# Apply the antithetic variable approach and moment matching to reduce variation

def Cholesky_reduce(K, r, T, n, noS, S0List, qList, sigmaList, antitheticList):

    for i in range(noS // 2):  # Antithetic variable approach
        temp3 = []
        for j in range(n):
            temp3.append((-1 * antitheticList[i][j]))
        antitheticList.append(temp3)
    # a = np.mat(antitheticList)
    # print(a[4999], a[9999])

    for i in range(n):
        tempsum2 = 0
        for j in range(noS):
            tempsum2 = tempsum2 + (antitheticList[j][i] ** 2)  # through antithetic approach, the mean is 0 now

        tempsum2 = np.sqrt(tempsum2 / (noS - 1))
        for k in range(noS):
            antitheticList[k][i] = antitheticList[k][i] / tempsum2
    saveList = copy.deepcopy(antitheticList)

    ST_List1 = []
    for i in range(noS):
        zeroList1 = []
        for j in range(n):
            zeroList1.append(0)
        ST_List1.append(zeroList1)

    for i in range(noS):
        # print(antitheticList[i])
        for k in range(n):
            rsum2 = 0
            for x in range(n):
                rsum2 = rsum2 + (antitheticList[i][x] * A[x][k])
                # print(antitheticList[i][x], A[x][k])
            lnST = np.log(S0List[k]) + (r - qList[k] - (sigmaList[k] * sigmaList[k] * 0.5)) * T + rsum2
            ST_List1[i][k] = np.exp(lnST)

    sum_of_payoff2 = 0
    for i in range(noS):
        payoff = max(max(ST_List1[i]) - K, 0)
        sum_of_payoff2 = sum_of_payoff2 + payoff

    Cho_call_reduce = np.exp(-r * T) * (sum_of_payoff2 / noS)

    return Cho_call_reduce, saveList


call_price_reduce, saveList = Cholesky_reduce(K, r, T, n, noS, S0List, qList, sigmaList, list_for_mm)
print("When with antithetic variate and moment matching, the rainbow option price would be: " + str(call_price_reduce))

# Implement inverse Cholesky method

inverseList = []
for i in range(n):
    tempList = []
    for j in range(noS):
        tempList.append(saveList[j][i])
    inverseList.append(tempList)
inverseList = np.mat(inverseList)

x = np.vstack(inverseList)
cov = np.cov(x)
cov = cov.tolist()
if isinstance(cov, float):
    cov = [[cov]]

A2 = []
for i in range(n):  # apply the algorithm to compute A2
    zeroList = []
    for j in range(n):
        zeroList.append(0)
    A2.append(zeroList)

for i in range(n):
    for j in range(n):
        if i == 0:
            if j == 0:
                A2[i][j] = np.sqrt(cov[i][j])
            else:
                A2[i][j] = cov[i][j] / A2[i][i]
        else:
            if i == j:
                tempsum5 = 0
                for k in range(i):
                    tempsum5 = tempsum5 + (A2[k][i] ** 2)
                A2[i][j] = np.sqrt(cov[i][j] - tempsum5)
            elif j > i:
                tempsum5 = 0
                for k in range(i):
                    tempsum5 = tempsum5 + (A2[k][i] * A2[k][j])
                A2[i][j] = (cov[i][j] - tempsum5) / A2[i][i]

A2 = np.mat(A2)
inverse_A2 = np.linalg.inv(A2)
inverse_A2 = inverse_A2.tolist()

def inverse_Cholesky(K, r, T, n, noS, S0List, qList, sigmaList, saveList):

    for i in range(noS):
        for j in range(n):
            rsum3 = 0
            for k in range(n):
                rsum3 = rsum3 + saveList[i][k] * inverse_A2[j][k]
            saveList[i][j] = rsum3

    ST_List2 = []
    for i in range(noS):
        zeroList2 = []
        for j in range(n):
            zeroList2.append(0)
        ST_List2.append(zeroList2)

    for i in range(noS):
        for j in range(n):
            rsum3 = 0
            for k in range(n):
                rsum3 = rsum3 + saveList[i][k] * A[k][j]
            lnST = np.log(S0List[j]) + (r - qList[j] - (sigmaList[j] * sigmaList[j] * 0.5)) * T + rsum3
            ST_List2[i][j] = np.exp(lnST)

    sum_of_payoff3 = 0
    for i in range(noS):
        payoff = max(max(ST_List2[i]) - K, 0)
        sum_of_payoff3 = sum_of_payoff3 + payoff

    call_price_inverse = np.exp(-r * T) * (sum_of_payoff3 / noS)

    return call_price_inverse


call_price_inverse = inverse_Cholesky(K, r, T, n, noS, S0List, qList, sigmaList, saveList)
print("When implement inverse Cholesky, the rainbow option price would be: " + str(call_price_inverse))

# Compute the confidence interval
confidence_list = [call_price]
confidence_list_reduce = [call_price_reduce]
confidence_list_inverse = [call_price_inverse]

for j in range(noR - 1):
    Ans, list_for_mm2 = Cholesky(K, r, T, n, noS, S0List, qList, sigmaList)
    confidence_list.append(Ans)
    Ans_reduce, saveList2 = Cholesky_reduce(K, r, T, n, noS, S0List, qList, sigmaList, list_for_mm2)
    confidence_list_reduce.append(Ans_reduce)
    Ans_inverse = inverse_Cholesky(K, r, T, n, noS, S0List, qList, sigmaList, saveList2)
    confidence_list_inverse.append(Ans_inverse)

sum_for_interval = 0
for a in confidence_list:
    sum_for_interval = sum_for_interval + a

average = sum_for_interval / len(confidence_list)

diff_for_interval = 0
for a in confidence_list:
    diff_for_interval = diff_for_interval + ((a - average) ** 2)

standard_deviation = np.sqrt(diff_for_interval / (len(confidence_list) - 1))

upper = average + 2 * standard_deviation
lower = average - 2 * standard_deviation

print("95% Confidence Interval would be: " + "[" + str(lower) + " - " + str(upper) + "]")

sum_for_interval_reduce = 0
for b in confidence_list_reduce:
    sum_for_interval_reduce = sum_for_interval_reduce + b

average_reduce = sum_for_interval_reduce / len(confidence_list_reduce)

diff_for_interval_reduce = 0
for b in confidence_list_reduce:
    diff_for_interval_reduce = diff_for_interval_reduce + ((b - average_reduce) ** 2)

standard_deviation_reduce = np.sqrt(diff_for_interval_reduce / (len(confidence_list_reduce) - 1))

upper_reduce = average_reduce + 2 * standard_deviation_reduce
lower_reduce = average_reduce - 2 * standard_deviation_reduce

print("95% Confidence Interval would be: " + "[" + str(lower_reduce) + " - " + str(upper_reduce) + "]")

sum_for_interval_inverse = 0
for b in confidence_list_inverse:
    sum_for_interval_inverse = sum_for_interval_inverse + b

average_inverse = sum_for_interval_inverse / len(confidence_list_inverse)

diff_for_interval_inverse = 0
for b in confidence_list_inverse:
    diff_for_interval_inverse = diff_for_interval_inverse + ((b - average_inverse) ** 2)

standard_deviation_inverse = np.sqrt(diff_for_interval_inverse / (len(confidence_list_inverse) - 1))

upper_inverse = average_inverse + 2 * standard_deviation_inverse
lower_inverse = average_inverse - 2 * standard_deviation_inverse

print("95% Confidence Interval would be: " + "[" + str(lower_inverse) + " - " + str(upper_inverse) + "]")

