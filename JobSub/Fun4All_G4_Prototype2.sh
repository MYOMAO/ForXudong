#!/bin/csh
setenv HOME /star/u/$LOGNAME
#setenv HOME /sphenix/user/$LOGNAME

source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
 source $i
end

source $HOME/.login
#source /direct/star+u/zshi/.login

#source /cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.csh -n ana.141

#source /cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.csh -n
#source /cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_root6.csh

#source /opt/sphenix/core/bin/sphenix_setup.csh -n
#source /opt/sphenix/core/bin/setup_root6.csh

#source /opt/sphenix/core/bin/setup_root6_include_path.csh


echo "START PRINT ENV"

#printenv


echo "DONE PRINt ENV"

set NEvent=$argv[1]
set Name=$argv[2]



#source /opt/sphenix/core/bin/sphenix_setup.csh -n





#source Build.sh


echo "Now PWD"

pwd

ls

echo "DONE CHECK"

cd workdir

mkdir ${Name}

#source Reconnect.sh


cp -r ../macros/ ${Name}/ 
cp ../Reconnect.sh  ${Name}/


cd ${Name}

echo "NowList"



#source BuildJob.sh


#setenv ROOT_INCLUDE_PATH /sphenix/user/zshi/FastMLWork7/JobSub/workdir/${Name}/macros/common:$ROOT_INCLUDE_PATH





#cd HFMLTrigger_LANL 

source Reconnect.sh



#setenv ROOT_INCLUDE_PATH /sphenix/user/zshi/EvtGenTestJobSub/workdir/${Name}/macros/common:$ROOT_INCLUDE_PATH


echo "DONE BUILD"

cd macros/detectors/sPHENIX/

ls *


#rm MyQAFile.root



root -b -l -q Fun4All_G4_sPHENIX.C'('${NEvent}')'

ls *root


#mv MyQAFile.root  ../../../../../OutFiles/MyQAFile_${Name}.root
#mv Test2.json ../../../../../OutFiles/Background/D0Background_${Name}.json


#mv Test2.json ../../../../../OutFiles/Signal/D0Signal_${Name}.json

#mv Test2.json ../../../../../OutFiles/NewSignal/D0Signal_${Name}.json

#mv BRFile.root ../../../../../OutFiles/BRFile_${Name}.root

#mv G4EICDetector.root_g4femc_eval.root ../../../ERECO/${Material}/${Rad}/${Energy}/G4EICDetector.root_g4femc_eval_${Name}.root


#mv G4sPHENIX.root /sphenix/tg/tg01/hf/zshi/BhadronAna/DST/G4sPHENIX_${Name}.root 

mv G4sPHENIX_g4svtx_eval.root /sphenix/tg/tg01/hf/zshi/BhadronAna/Eval/G4sPHENIX_g4svtx_eval_${Name}.root 

cd KFParticle/

source Reconnect.sh


sh run_MDC2reco.sh KF.list

mv ../G4sPHENIX.root /sphenix/tg/tg01/hf/zshi/BhadronAna/DST/G4sPHENIX_${Name}.root 

mv MyBPTest/outputTruthDecayTest_MyBPTest_KF.root  /sphenix/tg/tg01/hf/zshi/BhadronAna/KFOutPut/BPKF_${Name}.root 

cd ../../../../../

rm -rf ${Name}



