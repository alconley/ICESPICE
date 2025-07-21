#!/bin/bash

set -e

echo "📥 Cloning Geant4 source code..."
git clone https://github.com/Geant4/geant4.git
cd geant4/
git checkout v11.3.2
cd ..

echo "📁 Creating build directory..."
mkdir -p geant4_build/
cd geant4_build/

echo "📁 Creating install directory..."
INSTALL_DIR="/usr/local/geant4/11.3.2"
if [ ! -d "$INSTALL_DIR" ]; then
    sudo mkdir -p "$INSTALL_DIR"
fi

echo "🔍 Finding latest Homebrew Qt5 and Xerces-C..."
QT_PREFIX=$(brew --prefix qt@5)
XERCESC_PREFIX=$(brew --prefix xerces-c)
SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
CLANG=$(xcrun --find clang)
CLANGXX=$(xcrun --find clang++)

echo "📦 Qt prefix:        $QT_PREFIX"
echo "📦 Xerces-C prefix:  $XERCESC_PREFIX"
echo "🧭 macOS SDK root:   $SDKROOT"
echo "🛠  Using clang:     $CLANG"
echo "🛠  Using clang++:   $CLANGXX"

# Set compilers and paths
export CC="$CLANG"
export CXX="$CLANGXX"
export SDKROOT="$SDKROOT"
export CMAKE_FRAMEWORK_PATH="/System/Library/Frameworks"
export PKG_CONFIG_PATH="$QT_PREFIX/lib/pkgconfig:$XERCESC_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export PATH="$QT_PREFIX/bin:$PATH"

echo "🔧 Configuring Geant4 with CMake..."
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_OSX_SYSROOT="$SDKROOT" \
      -DCMAKE_FRAMEWORK_PATH="$CMAKE_FRAMEWORK_PATH" \
      -DCMAKE_C_COMPILER="$CC" \
      -DCMAKE_CXX_COMPILER="$CXX" \
      -DGEANT4_USE_GDML=ON \
      -DGEANT4_BUILD_MULTITHREADED=ON \
      -DXERCESC_ROOT_DIR="$XERCESC_PREFIX" \
      -DGEANT4_USE_QT=ON \
      -DGEANT4_INSTALL_EXAMPLES=ON \
      -DGEANT4_INSTALL_DATA=ON \
      -DGEANT4_USE_SYSTEM_EXPAT=OFF \
      -DGEANT4_BUILD_TLS_MODEL=auto \
      ../geant4

echo "⚙️  Building and installing Geant4..."
sudo cmake --build . --target install -- -j$(sysctl -n hw.logicalcpu)

echo "✅ Geant4 v11.3.2 build and install complete!"