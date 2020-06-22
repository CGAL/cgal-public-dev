import os
import subprocess
import pandas as pd
import argparse

parser = argparse.ArgumentParser() 
parser.add_argument("-d", "--dataPath", help = "Path to the data.")
parser.add_argument('-r','--rayList', nargs='+', help='Number of rays to be shot.', required=True)
parser.add_argument('-p','--point', nargs='+', help='Point to shoot the rays from.', required=True)

args = parser.parse_args()

programs = ['./AabbTest', './EmbreeTest', './AabbParallel', './EmbreeParallel']
dfList = []

for program in programs:
    tempList = []
    for rays in args.rayList:
        command  = program + " " + args.dataPath + " " + rays + " " + args.point[0] + " " + args.point[1] + " " + args.point[2]
        arglist = command.split(" ")
        process = subprocess.run(arglist, 
                         stdout=subprocess.PIPE, 
                         universal_newlines=True)
        funcIndex = process.stdout.find("Function() time:")+17
        tempList.append(process.stdout[funcIndex:-1])
    dfList.append(tempList)

dfdata = {'No of Rays' : args.rayList}
for i in range (len(programs)):
    dfdata.update({programs[i][2:] : dfList[i]} )

df = pd.DataFrame(data = dfdata)
df.to_csv("OutTable.csv", index = False)

print(df)
