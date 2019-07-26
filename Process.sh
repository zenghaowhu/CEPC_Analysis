#!/bin/bash
path=`pwd`

cd $path
#pars=("uusd" "ccds" "cuxx" "cc_nots" "dtdt" "utut" "uu_notd")
pars=("qqh_bb" "qqh_cc" "qqh_gg" "ww_h0ccbs" "ww_h0cuxx" "ww_h0uusd" "ww_h0ccds" "zz_h0dtdt" "zz_h0utut" "zz_h0uu_notd" "zz_h0cc_nots")
#Rparas=("1" "2" "2.5")
#Pparas=("-1" "0" "1" "2" "3")

#Rparas=("1" "1.3" "1.5" "2" "2.5")
#Pparas=("-1" "0" "0.5" "1" "1.5" "2" "3")
#Gparas=("-1" "-0.5" "0" "0.5" "1" "2")
#dirs=("div4" "div5" "div6" "div7" "div8" "div10" "div12")



ipar=3
#while [ "$ipar" -lt "7" ]
#do
par=${pars[$ipar]}

#idir=0
#while [ "$idir" -lt "1" ]
#do
#dir=${dirs[$idir]}


#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/ztoqq
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/350GeV/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/4fermions/${par}
jobPath=/cefs/higgs/zengh/Reconstruction/E240.Pww_h.${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/jetclustering/ww/${par}

#numSlcio=`find ${jobPath} -name "*.slcio" |wc -l`
#NumSlcio=`ls ${jobPath}/RecoJet/*.slcio |wc -l`
#NumSlcio=`ls ${jobPath}/MCPSLCIO/*.slcio |wc -l`
NumSlcio=`ls ${jobPath}/*.slcio |wc -l`
numSlcio=$((10#${NumSlcio}+1))
#value=$((10#${a}+32))

#iRpara=0
#while [ "$iRpara" -lt "3" ]
#do
#Rpara=${Rparas[$iRpara]}

#iPpara=0
#while [ "$iPpara" -lt "5" ]
#do
#Ppara=${Pparas[${iPpara}]}

#iGpara=2
#while [ "$iGpara" -lt "3" ]
#do
#Gpara=${Gparas[${iGpara}]}


j=1
while [ "$j" -lt "${numSlcio}" ]
#while [ "$j" -lt "2" ]
do

# if [ "$j" -lt 10 ]
# then i="0000"$j
# elif [ "$j" -ge 10 -a "$j" -lt 100 ]
# then i="000"$j
# elif [ "$j" -ge 100 -a "$j" -lt 1000 ]
# then i="00"$j
# fi
i=`printf "%05d" $j`

export RecoWorkDir=${jobPath}
#OUTPUTDATA=newSigma/afterPair/Out${Rpara}:${Ppara}
#OUTPUTDATA=rootFiles/Out${Rpara}${Ppara}
OUTPUTDATA=rootFiles
echo $RecoWorkDir/$OUTPUTDATA/
mkdir -p $RecoWorkDir/$OUTPUTDATA/
#outlog=$RecoWorkDir/$OUTPUTDATA/recolog_${j}GeV.log
cp -fr /afs/ihep.ac.cn/users/z/zengh/script/Analysis/knight.xml  $RecoWorkDir/$OUTPUTDATA/${par}_${i}.xml
#inputlcio=../Coral/result/Pe2e2h_e3e3.eU.pU.02.${j}.slcio 
#if [ $j -lt 10 ]
#then
#inputlcio=${jobPath}/VLC/RecoJet/Out${Rpara}:${Ppara}:${Gpara}/Reco_${i}.slcio
#inputlcio=${jobPath}/RecoJet/Out${Rpara}-${Ppara}/RECO_${i}.slcio 
#inputlcio=${jobPath}/RecoJet/Out${Rpara}${Ppara}/Reco_${i}.slcio
#inputlcio=${jobPath}/RecoJet/Out21/Reco_${i}.slcio
#inputlcio=${jobPath}/RecoJet/Reco_${i}.slcio
inputlcio=${jobPath}/${par}.e0.p0.${i}_rec.slcio
#inputlcio=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/zPole/MCPSLCIO/${i}.slcio
#inputlcio=${jobPath}/Reco_${i}.slcio
#inputlcio=${jobPath}/240GeV2ferminosdd.slcio
#inputlcio=${jobPath}/ww_semiLeptonic_00001.slcio

echo ${inputlcio}

#inputlcio=${jobPath}/${dir}/RecoJettwo/Reco_${i}.slcio
#inputlcio=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/MCPnngg/twelve/${par}/${i}.slcio
#inputlcio=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/zPole/MCPSLCIO/${i}.slcio
outputroot=$RecoWorkDir/$OUTPUTDATA/${par}.e0.p0.${i}.root
#outputroot=$RecoWorkDir/$OUTPUTDATA/240GeV2ferminosdd.root
#outputroot=$RecoWorkDir/$OUTPUTDATA/ww_semiLeptonic_00001.root 

sed -i "s#LCIOINPUT#${inputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${par}_${i}.xml
sed -i "s#OUTROOT#${outputroot}#g" $RecoWorkDir/$OUTPUTDATA/${par}_${i}.xml


echo "$RecoWorkDir/$OUTPUTDATA/${par}_${i}.xml"
export ILC_HOME=/afs/ihep.ac.cn/soft/common/gcc/v01-17-05
# the following code used to write a sh file to be hep_sub
echo \
"#! /bin/bash
unset MARLIN_DLL
source /afs/ihep.ac.cn/users/z/zengh/zhuyf/analysis/ArborGeneralLICH/env.sh
Marlin  $RecoWorkDir/$OUTPUTDATA/${par}_${i}.xml
#> $RecoWorkDir/$OUTPUTDATA/${par}_${i}.log

" > $RecoWorkDir/$OUTPUTDATA/${par}_${i}_ana.sh
echo "$RecoWorkDir/$OUTPUTDATA/${par}_${i}_ana.sh"
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
chmod +x $RecoWorkDir/$OUTPUTDATA/${par}_${i}_ana.sh
hep_sub $RecoWorkDir/$OUTPUTDATA/${par}_${i}_ana.sh 
#sh $RecoWorkDir/$OUTPUTDATA/${par}_${i}.sh 
let "j+=1"
done

#let "idir+=1"
#done

#let "iGpara+=1"
#done

#let "iPpara+=1"
#done

#let "iRpara+=1"
#done

# let "idir+=1"
#done


#let "ipar+=1"
#done
