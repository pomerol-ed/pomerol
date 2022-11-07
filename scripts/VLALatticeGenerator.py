#
# This file is part of pomerol, an exact diagonalization library aimed at
# solving condensed matter models of interacting fermions.
#
# Copyright (C) 2016-2022 A. Antipov, I. Krivenko and contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json

U=4.0
t=[0.5,1,0.3,0.4]
n=len(t)+1

print "Creating Lattice file from scratch for VLA calc with 1 correlated and " + str(n-1) + " bath sites"

filename="LatticeTest.json"

lattice=dict()
lattice["sites"]={}
sites=[]
lattice["sites"]["0"]={}
Master=lattice["sites"]["0"]
Master["U"]=U
Master["LocalMu"]=U/2
Master["type"]="s"
Master["hopping"]=[]
for j in range(1,n):
    Master["hopping"].append({})
    Master["hopping"][j-1]["value"]=t[j-1]
    Master["hopping"][j-1]["to"]=str(j)
    Master["hopping"][j-1]["orbital_from"]=0
    Master["hopping"][j-1]["orbital_to"]=0


for i in range(1,n):
    sites.append(str(i))
    lattice["sites"][str(i)]={}
    lattice["sites"][str(i)]["U"]=0.0
    lattice["sites"][str(i)]["LocalMu"]=0.0
    lattice["sites"][str(i)]["type"]="s"
    lattice["sites"][str(i)]["hopping"]=[]
    lattice["sites"][str(i)]["hopping"].append({})
    lattice["sites"][str(i)]["hopping"][0]["value"]=t[i-1]
    lattice["sites"][str(i)]["hopping"][0]["orbital_from"]=0
    lattice["sites"][str(i)]["hopping"][0]["orbital_to"]=0
    lattice["sites"][str(i)]["hopping"][0]["to"]="0"


output_file=open(filename,"w")
output_file.write(json.dumps(lattice, sort_keys=True, indent=4))
output_file.close()
print "Saved to " + filename
