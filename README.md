# X(3872) Simulation with sPHENIX


## Running the Codes

Follow the steps below to run a few events interactively:

git clone https://github.com/MYOMAO/ForXudong.git

source Build.sh

cd macros/detectors/sPHENIX

Look for the code 

https://github.com/MYOMAO/ForXudong/blob/master/macros/detectors/sPHENIX/Fun4All_G4_sPHENIX.C#L150

So this set the 1 B+ in 1 event with kinematics defined as

https://github.com/MYOMAO/ForXudong/blob/master/macros/detectors/sPHENIX/Fun4All_G4_sPHENIX.C#L151-L167

So basically the B+ is generated with a Gaussian distribution centered at (0,0,0) with width (0.01,0.01,5) cm 

The B+ is then handled by EvtGen. The EvtGen forces B+ to decay into D0bar pi+

https://github.com/MYOMAO/ForXudong/blob/master/macros/detectors/sPHENIX/BPD0Pi.DEC

The default settings allow tracking systems: MVTX + INTT + TPC + TPOT and tracking algorithm to run

You can try to run 100 events:

root -b -l -q Fun4All_G4_sPHENIX.C'(100)'


The output of the tracking files is saved at:


root -l G4sPHENIX_g4svtx_eval.root

You will have the following branches:

KEY: TNtuple	ntp_info;1	event info
KEY: TNtuple	ntp_vertex;1	vertex => max truth
KEY: TNtuple	ntp_gpoint;1	g4point => best vertex
KEY: TNtuple	ntp_cluster;1	svtxcluster => max truth
KEY: TNtuple	ntp_g4cluster;1	g4cluster => max truth
KEY: TNtuple	ntp_gtrack;1	g4particle => best svtxtrack
KEY: TNtuple	ntp_track;1	svtxtrack => max truth


ntp_vertex: the reconstructed vertex with tracking (Primary vertex only)
ntp_gpoint: the truth (primary) vertex 
ntp_gtrack: the truth primary particles generated in the simulation 
ntp_track: all reconstructed tracks with truth particle matched. Track parameters are evaluated at the nearest reconstructed primary vertex.


Play with codes and look at the TBranch in the TTree and tell me what you learn. :)
