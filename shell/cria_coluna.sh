#!/bin/bash

infile=posicoes_proab1.dat

grep -v %  $infile | grep -v '^\s*$' | awk '{ print $1,"&",$3"grau", $4"min","&",$5"grau", $6"min","&",$7,"&",$8":"$9",",$10"/"$11"/20"$12,"&","CTD"$2"\\\\" }' | sed "s/min/'/g" | sed "s:grau:"$^\circ'$'":g" | sed "s/CTD1/CTD/g" | sed "s/CTD2/XBT/g" | sed "s:CTD3:CTD/XBT:g" 


#linha de exemplo
# 1 & 44$^\circ$ 34.96' & 24$^\circ$ 30.02' & 149 & 25/01/08, 10:51 & CTD/XBT\\
