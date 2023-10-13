
cd coresoftware/offline/framework/fun4all/ 

make clean
sh autogen.sh --prefix=$MYINSTALL
make -j8 install

cd ../../../../

echo "DONE BRO"
