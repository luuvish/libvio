#!/bin/sh
# ==============================================================================
#
#   This confidential and proprietary software may be used only
#  as authorized by a licensing agreement from Thumb o'Cat Inc.
#  In the event of publication, the following notice is applicable:
#
#       Copyright (C) 2013 - 2014 Thumb o'Cat
#                     All right reserved.
#
#   The entire notice above must be reproduced on all authorized copies.
#
# ==============================================================================
#
#  File      : bootstrap.sh
#  Author(s) : Luuvish
#  Version   : 1.0
#  Revision  :
#      1.0 May 18, 2014    3rd party build setup
#
# ==============================================================================

SCRIPT_DIR=$( (cd -P $(dirname $0) && pwd) )
SOURCE_DIR="$SCRIPT_DIR/../tool/tarball"
TARGET_DIR="$SCRIPT_DIR/../tool/3rd-party"

[ ! -d $TARGET_DIR ] && mkdir -p $TARGET_DIR


# ================================
# FFmpeg source download and build

function bootstrap_ffmpeg {
  FFMPEG_URL="git://source.ffmpeg.org/ffmpeg.git"
  FFMPEG_VER="2.2.1"
  FFMPEG_SRC="$TARGET_DIR/ffmpeg-$FFMPEG_VER.git"
  FFMPEG_BIN="$TARGET_DIR/ffmpeg-$FFMPEG_VER.bin"

  if [ ! -d $FFMPEG_SRC ]; then
    git clone $FFMPEG_URL $FFMPEG_SRC
    [ ! -d $FFMPEG_SRC ] && exit 1
    (cd $FFMPEG_SRC; git checkout n$FFMPEG_VER; cd $SCRIPT_DIR)
  fi

  (
    mkdir -p $FFMPEG_BIN/dist;
    cd $FFMPEG_BIN;
    $FFMPEG_SRC/configure --prefix="$FFMPEG_BIN/dist";
    make && make install;
    cd $SCRIPT_DIR;
  )
}


# ================================
# libvpx source download and build

function bootstrap_libvpx {
  LIBVPX_URL="http://git.chromium.org/webm/libvpx.git"
  LIBVPX_VER="1.3.0"
  LIBVPX_SRC="$TARGET_DIR/libvpx-$LIBVPX_VER.git"
  LIBVPX_BIN="$TARGET_DIR/libvpx-$LIBVPX_VER.bin"

  if [ ! -d $LIBVPX_SRC ]; then
    git clone $LIBVPX_URL $LIBVPX_SRC
    [ ! -d $LIBVPX_SRC ] && exit 1
    (cd $LIBVPX_SRC; git checkout v$LIBVPX_VER; cd $SCRIPT_DIR)
  fi

  (
    mkdir -p $LIBVPX_BIN/dist;
    cd $LIBVPX_BIN;
    $LIBVPX_SRC/configure --prefix="$LIBVPX_BIN/dist";
    make && make install;
    cd $SCRIPT_DIR;
  )
}


# ==================================
# openhevc source download and build

function bootstrap_openhevc {
  OPENHEVC_URL="https://github.com/OpenHEVC/openHEVC"
  OPENHEVC_SRC="$TARGET_DIR/openhevc.git"
  OPENHEVC_BIN="$TARGET_DIR/openhevc.bin"

  if [ ! -d $OPENHEVC_SRC ]; then
    git clone $OPENHEVC_URL $OPENHEVC_SRC
    [ ! -d $OPENHEVC_SRC ] && exit 1
    (cd $OPENHEVC_SRC; git checkout hevc; cd $SCRIPT_DIR)
  fi

  (
    mkdir -p $OPENHEVC_BIN; cd $OPENHEVC_BIN;
    cmake $OPENHEVC_SRC; make; cd $SCRIPT_DIR;
  )
}


# ==================================
# openh264 source download and build

function bootstrap_openh264 {
  OPENH264_URL="https://github.com/cisco/openh264"
  OPENH264_SRC="$TARGET_DIR/openh264.git"
  OPENH264_BIN="$TARGET_DIR/openh264.git"

  if [ ! -d $OPENH264_SRC ]; then
    git clone $OPENH264_URL $OPENH264_SRC
    [ ! -d $OPENH264_SRC ] && exit 1
    (cd $OPENH264_SRC; git checkout master; cd $SCRIPT_DIR)
  fi

  (cd $OPENH264_BIN; make; cd $SCRIPT_DIR)
}


# ==================================
# h.265 HM source download and build

function bootstrap_h265_hm {
  H265HM_URL="https://hevc.hhi.fraunhofer.de/svn/svn_HEVCSoftware/tags/HM-14.0+RExt-7.2/"
  H265HM_VER="14.0"
  H265HM_SRC="$TARGET_DIR/hm-$H265HM_VER"
  H265HM_BIN="$TARGET_DIR/hm-$H265HM_VER"

  if [ ! -d $H265HM_SRC ]; then
    svn checkout $H265HM_URL $H265HM_SRC
    [ ! -d $H265HM_SRC ] && exit 1
    (svn update r3977)
  fi

  (cd $H265HM_BIN/build/linux; make; cd $SCRIPT_DIR)
}


# ==================================
# h.264 JM source download and build

function bootstrap_h264_jm {
  H264JM_URL="http://iphome.hhi.de/suehring/tml/download"
  H264JM_VER="18.6"
  H264JM_ZIP="$TARGET_DIR/jm-$H264JM_VER.zip"
  H264JM_SRC="$TARGET_DIR/jm-$H264JM_VER"
  H264JM_BIN="$TARGET_DIR/jm-$H264JM_VER"

  if [ ! -f $H264JM_ZIP ]; then
    curl -L -o $H264JM_ZIP $H264JM_URL/jm$H264JM_VER.zip
    [ ! -f $H264JM_ZIP ] && exit 1
  fi
  if [ ! -d $H264JM_SRC ]; then
    (unzip $H264JM_ZIP; mv JM $H264JM_SRC)
    [ ! -d $H264JM_SRC ] && exit 1
  fi

  (cd $H264JM_BIN; make; cd $SCRIPT_DIR)
}


# =======================
# SMPTE VC-1 source build

function bootstrap_smpte_vc1 {
  SMPTE_VC1_VER="2008"
  SMPTE_VC1_ZIP="$SOURCE_DIR/VC-1-SourceCode-$SMPTE_VC1_VER.zip"
  SMPTE_VC1_SRC="$TARGET_DIR/smpte-vc1-$SMPTE_VC1_VER"
  SMPTE_VC1_BIN="$TARGET_DIR/smpte-vc1-$SMPTE_VC1_VER"

  if [ ! -d $SMPTE_VC1_SRC ]; then
    unzip $SMPTE_VC1_ZIP -d $SMPTE_VC1_SRC
    [ ! -d $SMPTE_VC1_SRC ] && exit 1
  fi

  (
    cd $SMPTE_VC1_BIN;
    patch -p0 < $SCRIPT_DIR/3rd-party/smpte-vc1.diff;
    cd decoder; make; cd ..;
    cd encoder; make; cd ..;
    cd $SCRIPT_DIR;
  )
}


# ================================
# chips&media coda960 source build

function bootstrap_coda960 {
  CODA960_VER="1.0.0"
  CODA960_ZIP="$SOURCE_DIR/coda960-v$CODA960_VER.tar.gz"
  CODA960_SRC="$TARGET_DIR/coda960-v$CODA960_VER"
  CODA960_BIN="$TARGET_DIR/coda960-v$CODA960_VER"

  if [ ! -d $CODA960_SRC ]; then
    (tar zxf $CODA960_ZIP; mv CODA960-v$CODA960_VER $CODA960_SRC)
    [ ! -d $CODA960_SRC ] && exit 1
  fi

  (
    cd $CODA960_BIN/design/ref_c;
    patch -p0 < $SCRIPT_DIR/3rd-party/coda960.diff;
    make -f $SCRIPT_DIR/3rd-party/coda960.mk;
    cd $SCRIPT_DIR;
  )
}


# ================================
# chips&media sparrow source build

function bootstrap_sparrow {
  SPARROW_VER="2896"
  SPARROW_ZIP="$SOURCE_DIR/sparrow-r$SPARROW_VER.tar.gz"
  SPARROW_SRC="$TARGET_DIR/sparrow-r$SPARROW_VER"
  SPARROW_BIN="$TARGET_DIR/sparrow-r$SPARROW_VER"

  if [ ! -d $SPARROW_SRC ]; then
    (tar zxf $SPARROW_ZIP; mv sparrow $SPARROW_SRC)
    [ ! -d $SPARROW_SRC ] && exit 1
  fi

  (
    cd $SPARROW_BIN/branches; mkdir -p $SPARROW_BIN/dist;
    cd sparrow_avs; make; mv build/sparrow_avs $SPARROW_BIN/dist/sparrow_avs; cd ..;
    cd sparrow_mp2; make; mv build/sparrow_mp2 $SPARROW_BIN/dist/sparrow_mp2; cd ..;
    cd sparrow_rvx; make; mv build/sparrow_rv  $SPARROW_BIN/dist/sparrow_rvx; cd ..;
    cd sparrow_vc1; make; mv sparrow           $SPARROW_BIN/dist/sparrow_vc1; cd ..;
    cd $SCRIPT_DIR;
  )
}


if [ "$1" == "3rd-party" ]; then
  bootstrap_ffmpeg
  bootstrap_libvpx
  bootstrap_openhevc
  bootstrap_openh264
  bootstrap_h265_hm
  bootstrap_h264_jm
  bootstrap_smpte_vc1
  bootstrap_coda960
  bootstrap_sparrow
fi
