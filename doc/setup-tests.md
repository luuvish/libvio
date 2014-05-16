# Swallow Test Scenario

## 

### Test Tables

  Codec |                   Models              |  Streams
--------|---------------------------------------|-------------
  H265  | HM-10.1                               | Ngcodec
  H264  | FFMPEG-1.2, JM-18.4, Allegro, CODA960 | JVT, Allegro
  VP9   |                                       |
  VP8   | FFMPEG-1.2, CODA960                   | Google
  VC1   | FFMPEG-1.2, SMPTE-VC1, WMV9, CODA960  | SMPTE, WMV9
  MPEG4 | FFMPEG-1.2, CODA960                   |
  MPEG2 | FFMPEG-1.2, CODA960                   |
  JPEG  | JPEG-V9, CODA960                      |

### H264

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

4. Compare tables of all Models


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


diff -urwN test/tools/VC-1-SourceCode-2008 test/tools/VC-1-SourceCode-2008-fix > smpte-vc1.diff
patch -p0 < smpte-vc1.diff