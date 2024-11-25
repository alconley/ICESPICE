#! /bin/bash

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

export PKG_CONFIG_PATH=/opt/homebrew/Cellar/qt@5/5.15.2_1/lib/pkgconfig/
export PATH=/opt/homebrew/Cellar/qt@5/5.15.2_1/bin:$PATH
export PATH="/opt/homebrew/opt/qt@5/bin:$PATH"

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/geant4/11.2.0 \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_BUILD_MULTITHREADED=ON \
    -DXERCESC_ROOT_DIR=/opt/homebrew/Cellar/xerces-c/3.2.3 \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_INSTALL_EXAMPLES=ON \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_USE_SYSTEM_EXPAT=OFF \
    -DGEANT4_BUILD_TLS_MODEL=auto \
    ../geant4

echo "installing"
sudo cmake --build . --target install -- -j10
echo "build done"