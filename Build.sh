#source /opt/sphenix/core/bin/sphenix_setup.csh -n new

source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.358

rm -r install

mkdir install
setenv MYINSTALL $PWD/install/
setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
set path = ( $MYINSTALL/bin $path )

cd coresoftware/simulation/g4simulation/


cd g4decayer/

make clean
sh autogen.sh --prefix=$MYINSTALL
make -j8 install

cd ../

cd g4main/

make clean
sh autogen.sh --prefix=$MYINSTALL
make -j8 install

cd ../

cd ../../../

setenv ROOT_INCLUDE_PATH /sphenix/user/zshi/LightFlavor/macros/common:$ROOT_INCLUDE_PATH
