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
#  File      : Makefile
#  Author(s) : Luuvish
#  Version   : 1.0
#  Revision  :
#      1.0 June 15, 2013    first release
#
# ===========================================================================


NAME=   h264dec

### include debug information: 1=yes, 0=no
DBG?= 0
### Generate 32 bit executable : 1=yes, 0=no
M32?= 0
### include O level optimization : 0-3
OPT?= 3
### Static Compilation
STC?= 0
### OPENMP support : 1=yes, 0=no
OPENMP?= 0


DEPEND= dependencies

BINDIR= ../../../bin
INCDIR= lcommon/inc ldecod/inc core parser decoder
SRCDIR= ldecod/src
OBJDIR= objects

ADDSRCDIR= lcommon/src

ifeq ($(STC),1)
ifeq ($(DBG),1)  ### Do not use static compilation for Debug mode
STC=0
STATIC=
else
STATIC= -static
endif
else
STATIC= 
endif

#check for LLVM and silence warnings accordingly
LLVM = $(shell $(CC) --version | grep LLVM)
ifneq ($(LLVM),)
	CFLAGS+=-Qunused-arguments
else
	CFLAGS+=-Wno-unused-but-set-variable
endif

LIBS=   -lm $(STATIC)
#CFLAGS+=  -std=gnu99 -pedantic -ffloat-store -fno-strict-aliasing -fsigned-char $(STATIC)
CFLAGS+=  -pedantic -ffloat-store -fno-strict-aliasing -fsigned-char $(STATIC)
FLAGS=  $(CFLAGS) -Wall $(INCDIR:%=-I%)

ifeq ($(M32),1)
FLAGS+=-m32
endif

ifeq ($(OPENMP),1)
  FLAGS+=-fopenmp
endif

OPT_FLAG = -O$(OPT)
ifeq ($(DBG),1)
SUFFIX= .dbg
FLAGS+= -g
else
SUFFIX=
FLAGS+= $(OPT_FLAG)
endif

OBJSUF= .o$(SUFFIX)

SRC=    $(wildcard $(SRCDIR)/*.cc)
ADDSRC= $(wildcard $(ADDSRCDIR)/*.cc)
PARSER= $(wildcard parser/*.cc)
DECODER= $(wildcard decoder/*.cc)
CORE= $(wildcard core/*.cc)
OBJ=    $(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o$(SUFFIX))
OBJ+=   $(ADDSRC:$(ADDSRCDIR)/%.cc=$(OBJDIR)/%.o$(SUFFIX))
OBJ+=   $(PARSER:parser/%.cc=$(OBJDIR)/%.o$(SUFFIX))
OBJ+=   $(DECODER:decoder/%.cc=$(OBJDIR)/%.o$(SUFFIX))
OBJ+=   $(CORE:core/%.cc=$(OBJDIR)/%.o$(SUFFIX))
BIN=    $(BINDIR)/$(NAME)$(SUFFIX)

.PHONY: default distclean clean tags depend

default: objdir_mk depend bin 

clean:
	@echo remove all objects
	@rm -rf $(OBJDIR)

distclean: clean
	@rm -f $(DEPEND) tags
	@rm -f $(BIN)

tags:
	@echo update tag table
	@ctags inc/*.h src/*.cc

bin:    $(OBJ)
	@echo
	@echo 'creating binary "$(BIN)"'
	@$(CC) $(FLAGS) -o $(BIN) $(OBJ) $(LIBS)
	@echo '... done'
	@echo

depend:
	@echo
	@echo 'checking dependencies'
	@$(SHELL) -ec '$(CC) $(FLAGS) -MM $(CFLAGS) $(INCDIR:%=-I%) $(SRC) $(ADDSRC) $(PARSER) $(DECODER) $(CORE) \
         | sed '\''s@\(.*\)\.o[ :]@$(OBJDIR)/\1.o$(SUFFIX):@g'\'' \
         >$(DEPEND)'
	@echo

$(OBJDIR)/%.o$(SUFFIX): $(SRCDIR)/%.cc
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): $(ADDSRCDIR)/%.cc
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): parser/%.cc
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): decoder/%.cc
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): core/%.cc
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

objdir_mk:
	@echo 'Creating $(OBJDIR) ...'
	@mkdir -p $(OBJDIR)

-include $(DEPEND)
