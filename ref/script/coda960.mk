# ===========================================================================
#
#   This confidential and proprietary software may be used only
#  as authorized by a licensing agreement from Thumb o'Cat Inc.
#  In the event of publication, the following notice is applicable:
# 
#       Copyright (C) 2013 - 2013 Thumb o'Cat
#                     All right reserved.
# 
#   The entire notice above must be reproduced on all authorized copies.
#
# ===========================================================================
#
#  File      : coda960.mk
#  Author(s) : Luuvish
#  Version   : 1.0
#  Revision  :
#      1.0 May 12, 2013    first release
#
# ===========================================================================

BIN_DIR = bin/darwin

ifeq "$(wildcard $(BIN_DIR))" ""
-include $(shell mkdir -p $(BIN_DIR))
endif


TARGET := $(BIN_DIR)/avcdecoder $(BIN_DIR)/avsdecoder $(BIN_DIR)/vpxdecoder
TARGET += $(BIN_DIR)/vc1decoder $(BIN_DIR)/rv10decoder
TARGET += $(BIN_DIR)/divx3decoder $(BIN_DIR)/mpeg4decoder $(BIN_DIR)/mpeg2decoder
TARGET += $(BIN_DIR)/jpegdecoder
TARGET += $(BIN_DIR)/avcencoder $(BIN_DIR)/mpeg4encoder $(BIN_DIR)/jpegencoder
TARGET += $(BIN_DIR)/bin2txt $(BIN_DIR)/noise $(BIN_DIR)/raw2mb $(BIN_DIR)/txt2bin

CC := gcc

CCFLAGS := -Isrc/dec_h264 -Isrc/dec_avs -Isrc/dec_vpx
CCFLAGS += -Isrc/dec_vc1 -Isrc/dec_rv10
CCFLAGS += -Isrc/dec_divx3 -Isrc/dec_mpeg4 -Isrc/dec_mpeg2
CCFLAGS += -Isrc/dec_jpeg -Isrc/com_jpg
CCFLAGS += -Isrc/enc_h264 -Isrc/enc_mpeg4 -Isrc/enc_jpeg
CCFLAGS += -Isrc/tools
#CCFLAGS += -std=c99 -mllvm

LDFLAGS := -lm -O2


all: $(TARGET)

clean:
	$(MAKE) -C src/dec_vpx clean
	-rm -rf $(TARGET)

$(BIN_DIR)/avcdecoder: src/dec_h264/*.c src/dec_h264/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_h264/*.c

$(BIN_DIR)/avsdecoder: src/dec_avs/*.c src/dec_avs/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_avs/*.c

$(BIN_DIR)/vpxdecoder: src/dec_vpx/blackbird
	-cp src/dec_vpx/blackbird $(BIN_DIR)/vpxdecoder

src/dec_vpx/blackbird: src/dec_vpx/src/codecs/host/*.c
	$(MAKE) -C src/dec_vpx

$(BIN_DIR)/vc1decoder: src/dec_vc1/*.c src/dec_vc1/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_vc1/*.c

$(BIN_DIR)/rv10decoder: src/dec_rv10/*.c src/dec_rv10/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_rv10/*.c

$(BIN_DIR)/divx3decoder: src/dec_divx3/*.c src/dec_divx3/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_divx3/*.c

$(BIN_DIR)/mpeg4decoder: src/dec_mpeg4/*.c src/dec_mpeg4/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_mpeg4/*.c

$(BIN_DIR)/mpeg2decoder: src/dec_mpeg2/*.c src/dec_mpeg2/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_mpeg2/*.c

$(BIN_DIR)/jpegdecoder: src/dec_jpeg/*.c src/dec_jpeg/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/dec_jpeg/*.c src/com_jpg/*.c

$(BIN_DIR)/avcencoder: src/enc_h264/*.c src/enc_h264/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/enc_h264/*.c

$(BIN_DIR)/mpeg4encoder: src/enc_mpeg4/*.c src/enc_mpeg4/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/enc_mpeg4/*.c

$(BIN_DIR)/jpegencoder: src/enc_jpeg/*.c src/enc_jpeg/*.h
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ src/enc_jpeg/*.c src/com_jpg/*.c

$(BIN_DIR)/bin2txt: src/tools/bin2txt.c
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ $<

$(BIN_DIR)/noise: src/tools/noise.c
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ $<

$(BIN_DIR)/raw2mb: src/tools/raw2mb.c
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ $<

$(BIN_DIR)/txt2bin: src/tools/txt2bin.c
	$(CC) $(LDFLAGS) $(CCFLAGS) -o $@ $<
