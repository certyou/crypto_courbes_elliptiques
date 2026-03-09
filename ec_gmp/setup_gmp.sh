sudo apt install build-essential
sudo apt install gcc m4
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz 
unxz gmp-6.2.1.tar.xz 
tar xvf gmp-6.2.1.tar 
cd gmp-6.2.1 
./configure 
make 
sudo make install
