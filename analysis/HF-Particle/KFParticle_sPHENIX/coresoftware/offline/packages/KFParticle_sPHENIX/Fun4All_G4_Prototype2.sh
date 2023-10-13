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

cp -r ../../macros/ ${Name} 
cp ../Reconnect.sh  ${Name}/


cd ${Name}

echo "NowList"

ls *








#cd HFMLTrigger_LANL 

source Reconnect.sh



echo "DONE BUILD"

cd macros/detectors/sPHENIX/

#rm MyQAFile.root



root -b -l -q Fun4All_G4_sPHENIX.C'('${NEvent}')'

ls *root


cd KFParticle 

source ReconnectLocal.sh


sh run_MDC2reco.sh KF.list

mv mypipiReco/outputTruthDecayTest_mypipiReco_KF.root  ../../../../../../OutFiles/Solution/ZZ/KSKF_${Name}.root 

echo "DONE ZZ Codes"

source ReconnectKF.sh

sh run_MDC2reco.sh KF.list

mv mypipiReco/outputTruthDecayTest_mypipiReco_KF.root  ../../../../../../OutFiles/Solution/Cameron/KSKF_${Name}.root 


cd ../../../../../

rm -rf ${Name}



