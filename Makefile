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
INCDIR= ldecod/inc
SRCDIR= ldecod/src
OBJDIR= ldecod/obj

ADDSRCDIR= lcommon/src
ADDINCDIR= lcommon/inc

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
CFLAGS+=  -std=gnu99 -pedantic -ffloat-store -fno-strict-aliasing -fsigned-char $(STATIC)
FLAGS=  $(CFLAGS) -Wall -I$(INCDIR) -I$(ADDINCDIR) -D __USE_LARGEFILE64 -D _FILE_OFFSET_BITS=64

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

SRC=    $(wildcard $(SRCDIR)/*.c) 
ADDSRC= $(wildcard $(ADDSRCDIR)/*.c)
OBJ=    $(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o$(SUFFIX)) $(ADDSRC:$(ADDSRCDIR)/%.c=$(OBJDIR)/%.o$(SUFFIX)) 
BIN=    $(BINDIR)/$(NAME)$(SUFFIX)

.PHONY: default distclean clean tags depend

default: messages objdir_mk depend bin 

messages:
ifeq ($(M32),1)
	@echo 'Compiling with M32 support...'
endif
ifeq ($(DBG),1)
	@echo 'Compiling with Debug support...'
	@echo 'Note static compilation not supported in this mode.'
endif
ifeq ($(STC),1)
	@echo 'Compiling with -static support...'
endif
ifeq ($(OPENMP),1)
	@echo 'Compiling with -fopenmp support...'
endif

clean:
	@echo remove all objects
	@rm -rf $(OBJDIR)

distclean: clean
	@rm -f $(DEPEND) tags
	@rm -f $(BIN)

tags:
	@echo update tag table
	@ctags inc/*.h src/*.c

bin:    $(OBJ)
	@echo
	@echo 'creating binary "$(BIN)"'
	@$(CC) $(FLAGS) -o $(BIN) $(OBJ) $(LIBS)
	@echo '... done'
	@echo

depend:
	@echo
	@echo 'checking dependencies'
	@$(SHELL) -ec '$(CC) $(FLAGS) -MM $(CFLAGS) -I$(INCDIR) -I$(ADDINCDIR) $(SRC) $(ADDSRC)                  \
         | sed '\''s@\(.*\)\.o[ :]@$(OBJDIR)/\1.o$(SUFFIX):@g'\''               \
         >$(DEPEND)'
	@echo

$(OBJDIR)/%.o$(SUFFIX): $(SRCDIR)/%.c
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): $(ADDSRCDIR)/%.c
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

objdir_mk:
	@echo 'Creating $(OBJDIR) ...'
	@mkdir -p $(OBJDIR)

-include $(DEPEND)

