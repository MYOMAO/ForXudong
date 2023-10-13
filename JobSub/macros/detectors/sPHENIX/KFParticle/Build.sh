source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.354

#source /cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/bin/sphenix_setup.csh -n ana.354

#source /opt/sphenix/core/bin/sphenix_setup.csh -n new

rm -r install

mkdir install
setenv MYINSTALL $PWD/install/
setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
set path = ( $MYINSTALL/bin $path )


#cd coresoftware/offline/framework/fun4all/

cd coresoftware/offline/packages/KFParticle_sPHENIX

make clean
sh autogen.sh --prefix=$MYINSTALL
make -j8 install


#cd ../trackreco/

#make clean
#sh autogen.sh --prefix=$MYINSTALL
#make -j8 install

cd ../../../../

echo "DONE BRO"
