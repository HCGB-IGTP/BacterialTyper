# configure, make, and install
configure --prefix=$PREFIX --with-zlib=$PREFIX
make -j${CPU_COUNT}
make install