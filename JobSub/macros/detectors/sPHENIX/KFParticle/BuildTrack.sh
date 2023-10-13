
#cd coresoftware/offline/framework/fun4all/ 

cd coresoftware/offline/packages/trackreco/

make clean
sh autogen.sh --prefix=$MYINSTALL
make -j8 install

cd ../../../../

echo "DONE BRO"
