import os
import time
import numpy as np
from utils import read_parameters
import sys


try:
    P = read_parameters(sys.argv[1])
except:
    print("Where parameters path??")
    sys.exit()

directory = P["OutputFolder"]
MAXITER = 3600
i=0
tot_num = P["muBtilde"][2]*P["muQtilde"][2]*P["muStilde"][2]
RunMode = P["RunMode"]
NumberTask = P["Ttilde"][2]
while True and (i<MAXITER):
    L = [] 
    for filename in os.listdir(directory):
        if filename.startswith("TEMP_unordered_inversion_"):
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as f:
                line_count = sum(1 for line in f)
                L.append((int(filename.split("_")[-1][:-4]), int(line_count)))
    os.system('clear')
    L = np.array(L)
    L[:] = L[np.argsort(L[:, 0])]
    m, M = np.argmin(L[:,1]), np.argmax(L[:,1])
    for j, e in enumerate(L):
        percent = e[1]/tot_num * 100
        if int(percent) == 100:
            print("Task number:", str(e[0]), ": Done!")
        else:
            if RunMode == 0:
                print("Temperature", str(e[0])+" over "+str(NumberTask)+", Completed:", np.round(percent, 3), "%")
            elif RunMode == 1:
                if j==m:
                    print("Task number: ", str(e[0])+",", np.round(percent, 3),"% <--- MIN")
                elif j==M:
                    print("Task number: ", str(e[0])+",", np.round(percent, 3),"% <--- MAX")
                else:
                    print("Task number: ", str(e[0])+",", np.round(percent, 3),"%")
    time.sleep(1)
    i+=1
