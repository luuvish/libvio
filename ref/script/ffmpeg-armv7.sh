#!/bin/sh

IOSSDK_ROOT='/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer'
IOSSDK_SDK="$IOSSDK_ROOT/SDKs/iPhoneOS6.1.sdk"
IOSSDK_CC="$IOSSDK_ROOT/usr/bin/arm-apple-darwin10-llvm-gcc-4.2"
IOSSDK_AS="./gas-preprocessor/gas-preprocessor.pl $IOSSDK_CC"

SCRIPT_DIR=$( (cd -P $(dirname $0) && pwd) )
FFMPEG_SRC="$SCRIPT_DIR/ffmpeg-1.2"
FFMPEG_BIN="$SCRIPT_DIR/ffmpeg-1.2-armv7"

# FFmpeg source download from Git Repo.

#if [ ! -d $FFMPEG_SRC ]; then
#    git clone git://source.ffmpeg.org/ffmpeg.git $FFMPEG_SRC
#fi
#[ ! -d $FFMPEG_SRC ] && exit 1
#
#(cd $FFMPEG_SRC; git pull; cd $SCRIPT_DIR)

# Create Build Directory
mkdir -p $FFMPEG_BIN/dist

cd $FFMPEG_BIN

# GAS preprocessor download
git clone https://github.com/yuvi/gas-preprocessor.git gas-preprocessor

# FFmpeg configure
$FFMPEG_SRC/configure \
    --prefix="$FFMPEG_BIN/dist" \
    --enable-cross-compile --target-os=darwin --arch=arm --cpu=cortex-a8 \
    --cc="$IOSSDK_CC" --as="$IOSSDK_AS" --sysroot="$IOSSDK_SDK" \
    --extra-cflags="-mfpu=neon -pipe -Os -gdwarf-2 -miphoneos-version-min=5.0" \
    --extra-ldflags="-arch armv7 -isysroot $IOSSDK_SDK -miphoneos-version-min=5.0" \
    --enable-optimizations --enable-pic --disable-yasm --disable-asm \
    --disable-doc --disable-logging --disable-debug --disable-stripping \
    --disable-ffmpeg --disable-ffplay --disable-ffprobe --disable-ffserver \
    --disable-postproc

# Compile and Install
make && make install

cd $SCRIPT_DIR
