# libvio test

## bootstrap 3rd-party reference models

```bash
> ./script/bootstrap 3rd-party
```


diff -urwN test/tools/VC-1-SourceCode-2008 test/tools/VC-1-SourceCode-2008-fix > smpte-vc1.diff
patch -p0 < smpte-vc1.diff



### Tables for codecs, models and test vectors

  Codec | Model                                       | Streams
--------|---------------------------------------------|-------------
  h265  | ffmpeg-2.2.1, HM-14.0+RExt-7.2              | bbc, Ngcodec
  h264  | ffmpeg-2.2.1, JM-18.6, Allegro, CODA960     | JVT, Allegro
  vp9   | ffmpeg-2.2.1, libvpx-1.3.0                  | libvpx
  vp8   | ffmpeg-2.2.1, libvpx-1.3.0, CODA960         | libvpx
  vc1   | ffmpeg-2.2.1, SMPTE-VC1, WMFdecode, CODA960 | SMPTE, WMV9
  mpeg4 | ffmpeg-2.2.1, CODA960                       |
  mpeg2 | ffmpeg-2.2.1, CODA960                       |
  jpeg  | JPEG-V9, CODA960                            |


1. Allegro

  JVT, Allegro streams을 디코딩하고 yuv 저장.  
  각각의 yuv에 대한 checksum 저장.  
  각각의 yuv를 frame 단위로 나누고 checksum 저장.

2. FFMPEG-1.2, JM-18.4

  Allegro Model로 checksum을 덤프하는 방식으로 checksum을 저장.  
  위의 방식을 자동화하는 환경 구축.

3. CODA960

  JVT, Allegro streams을 디코딩하고 frame 마다 checksum 저장.

4. Setup Auto Compare Checksum

  dump checksum of frames of all streams  
  run models and compare checksums  
  dump logs of compare results


### HEVC

1. test vectors

  download streams at `ftp://ftp.kw.bbc.co.uk/hevc/hm-11.0-anchors/bitstreams/`


### H264

1. test vectors

   JVT  
   Allegro

2. make digests of bitstreams from CODA960

  stream source: ./test/stream/h264/{jvt,allegro}/\*  
  digest target: ./test/digest/h264/{jvt,allegro}/\*.yuv.md5

  ```bash
  > ./script/testsuite.py digest-h264-coda960
  ```

3. test decoding using CODA960

  stream source: ./test/stream/h264/{jvt,allegro}/\*  
  digest golden: ./test/digest/h264/{jvt,allegro}/\*.yuv.md5

  ```bash
  > ./script/testsuite.py compare-h264-coda960
  ```

jvt/bp/CVSE2_Sony_B.yuv
jvt/bp/MR2_TANDBERG_B.yuv
jvt/bp/NRF_MW_B.yuv
jvt/bp/SVA_Base_A.yuv
jvt/bp/SVA_Base_B.yuv
jvt/bp/SVA_FM1_A.yuv
jvt/bp/sp1_bt_a.yuv
jvt/mp/CVSE2_Sony_B.yuv
jvt/mp/vc714.yuv
allegro/mp/Allegro_Inter_CABAC_10_L42_HD2K@60Hz_10.0.yuv
allegro/mp/Allegro_Inter_CAVLC_10_L42_HD2K@60Hz_10.0.yuv


### VC-1

1. test vectors

   SMPTE-VC1  
   Microsoft WMV9, Artificials

2. make digests of bitstreams from CODA960

  stream source: ./test/stream/vc1/\*\*/\*  
  digest target: ./test/digest/vc1/\*\*/\*.yuv.md5

  ```bash
  > ./script/testsuite.py digest-vc1-coda960
  ```

3. test decoding using CODA960

  stream source: ./test/stream/vc1/\*\*/\*  
  digest golden: ./test/digest/vc1/\*\*/\*.yuv.md5

  ```bash
  > ./script/testsuite.py compare-vc1-coda960
  ```


### VP8

1. test vectors

  download streams at `http://downloads.webmproject.org/test_data/libvpx/`  
  or `cd ./tool/3rd-party/libvpx-1.3.0.bin; make testdata;`

2. make digests of bitstreams from libvpx

  stream source: ./test/stream/vp8/\*.ivf  
  digest target: ./test/digest/vp8/\*.yuv.md5

  ```bash
  > ./script/testsuite.py digest-vp8-libvpx
  ```

3. test decoding using libvpx

  stream source: ./test/stream/vp8/\*.ivf  
  digest golden: ./test/digest/vp8/\*.yuv.md5

  ```bash
  > ./script/testsuite.py compare-vp8-libvpx
  ```


### VP9

1. test vectors

  download streams at `http://downloads.webmproject.org/test_data/libvpx/`  
  or `cd ./tool/3rd-party/libvpx-1.3.0.bin; make testdata;`

2. make digests of bitstreams from libvpx

  stream source: ./test/stream/vp9/\*.{webm,ivf}  
  digest target: ./test/digest/vp9/\*.yuv.md5

  ```bash
  > ./script/testsuite.py digest-vp9-libvpx
  ```

3. test decoding using libvpx

  stream source: ./test/stream/vp9/\*.{webm,ivf}  
  digest golden: ./test/digest/vp9/\*.yuv.md5

  ```bash
  > ./script/testsuite.py compare-vp9-libvpx
  ```
