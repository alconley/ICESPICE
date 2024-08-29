#! /bin/bash

############# prerequisites ######################
sudo apt update
sudo apt upgrade -y

sudo apt install -y build-essential
sudo apt install -y git git-man
sudo apt install -y curl
sudo apt install -y cmake

sudo apt install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
sudo apt install -y libxerces-c-dev

sudo apt install -y gfortran
sudo apt install -y python3-dev python3-numpy-dev

sudo apt install -y libxpm-dev libxft-dev libxmu-dev
sudo apt install -y libssl-dev

sudo apt-get install -y libpcre3-dev libglu1-mesa-dev \
                libftgl-dev libfftw3-dev libcfitsio-dev \
                libgraphviz-dev libxml2-dev libgsl-dev

#######################################################

echo "pull source code from github"  
git clone https://github.com/Geant4/geant4.git
cd geant4/
git checkout v11.2.0
cd ..

echo "make the build directory"
mkdir geant4_build/
cd geant4_build/

echo "make the install directory"
sudo mkdir /usr/local/geant4/11.2.0

# # mac os

# export PKG_CONFIG_PATH=/opt/homebrew/Cellar/qt@5/5.15.2_1/lib/pkgconfig/
# export PATH=/opt/homebrew/Cellar/qt@5/5.15.2_1/bin:$PATH
# export PATH="/opt/homebrew/opt/qt@5/bin:$PATH"

# cmake -DCMAKE_INSTALL_PREFIX=/usr/local/geant4/11.2.0 \
#     -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#     -DGEANT4_USE_GDML=ON \
#     -DGEANT4_BUILD_MULTITHREADED=ON \
#     -DXERCESC_ROOT_DIR=/opt/homebrew/Cellar/xerces-c/3.2.3 \
#     -DGEANT4_USE_QT=ON \
#     -DGEANT4_INSTALL_EXAMPLES=ON \
#     -DGEANT4_INSTALL_DATA=ON \
#     -DGEANT4_USE_SYSTEM_EXPAT=OFF \
#     -DGEANT4_BUILD_TLS_MODEL=auto \
#     ../geant4

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/geant4/11.2.0 \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_BUILD_MULTITHREADED=ON \
    -DGEANT4_USE_QT=ON \
    ../geant4

echo "installing"

sudo mkdir /usr/local/geant4/11.2.0/share/Geant4-11.2.0/data
sudo cp -r data/ /usr/local/geant4/11.2.0/share/Geant4-11.2.0/data/

sudo cmake --build . --target install -- -j10
echo "build done"