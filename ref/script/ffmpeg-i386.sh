#!/bin/sh

SCRIPT_DIR=$( (cd -P $(dirname $0) && pwd) )
FFMPEG_SRC="$SCRIPT_DIR/ffmpeg-1.2"
FFMPEG_BIN="$SCRIPT_DIR/ffmpeg-1.2-i386"

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

# FFmpeg configure
$FFMPEG_SRC/configure \
    --prefix="$FFMPEG_BIN/dist" \
    --enable-optimizations --enable-pic \
    --disable-doc --disable-logging --disable-debug --disable-stripping \
    --disable-postproc

# Compile and Install
make && make install

cd $SCRIPT_DIR
