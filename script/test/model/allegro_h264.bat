@echo off

@rem ===========================================================================
@rem
@rem   This confidential and proprietary software may be used only
@rem  as authorized by a licensing agreement from Thumb o'Cat Inc.
@rem  In the event of publication, the following notice is applicable:
@rem 
@rem       Copyright (C) 2013 - 2013 Thumb o'Cat
@rem                     All right reserved.
@rem 
@rem   The entire notice above must be reproduced on all authorized copies.
@rem
@rem ===========================================================================
@rem
@rem  File      : allegro_h264.bat
@rem  Author(s) : Luuvish
@rem  Version   : 1.0
@rem  Revision  :
@rem      1.0 May 19, 2013    first release
@rem
@rem ===========================================================================


set streams_dir=..\..\streams\h264
set images_dir=..\..\images\h264
set log_file=allegro.log

echo. > %log_file%

goto :decode


:generate_option
	echo %2                       ........H.264/AVC coded bitstream > %1
	echo %3                       ........Output file, YUV/RGB >> %1
	echo %~n3.rec                 ........Ref sequence (for SNR) >> %1
	echo 1                        ........Write 4:2:0 chroma components for monochrome streams >> %1
	echo 0                        ........NAL mode (0=Annex B, 1: RTP packets) >> %1
	echo 0                        ........SNR computation offset >> %1
	echo 2                        ........Poc Scale (1 or 2) >> %1
	echo 500000                   ........Rate_Decoder >> %1
	echo 104000                   ........B_decoder >> %1
	echo 73000                    ........F_decoder >> %1
	echo leakybucketparam.cfg     ........LeakyBucket Params >> %1
	echo 0                        ........Err Concealment(0:Off,1:Frame Copy,2:Motion Copy) >> %1
	echo 2                        ........Reference POC gap (2: IPP (Default), 4: IbP / IpP) >> %1
	echo 2                        ........POC gap (2: IPP /IbP/IpP (Default), 4: IPP with frame skip = 1 etc.) >> %1
	echo 0                        ........Silent decode >> %1
	echo. >> %1
	echo This is a file containing input parameters to the JVT H.264/AVC decoder. >> %1
	echo The text line following each parameter is discarded by the decoder. >> %1
	echo. >> %1
	goto :eof

:decode_file
	call :generate_option %~n1.opt %1 %~dpn2.yuv
	call ldecod_Allegro.exe %~n1.opt
	del %~n1.opt %~n2.rec dataDec.txt log.dec
	goto :eof


:decode_echo
	echo ldecod -i %~nx1 -o %~n2.yuv
	goto :eof

:compare_file
	echo compare %1 %2
	call cmp %1 %2

	if not errorlevel 1 (
	 	echo. >> %log_file%
		echo Two file ^(%1, %2^) is same. >> %log_file%
		echo %date% %time% >> %log_file%
		echo. >> %log_file%
		del %1 %2
	) else (
	 	echo. >> %log_file%
	 	echo Two file ^(%1, %2^) is different. >> %log_file%
		echo %date% %time% >> %log_file%
		echo. >> %log_file%
	)
	set /a errorlevel = 0
	goto :eof

:decode_folder
	if not exist %~dpn2 (mkdir %~dpn2)
	for /f %%p in ('dir /b /o:n %1') do (
		call :decode_echo %~dp1%%p %~dpn2%%p
		call :decode_file %~dp1%%p %~dpn2%%p > %~n1%%p.log
		type %~n1%%p.log >> %log_file%
		del %~n1%%p.log
	)
	goto :eof

:decode
	call :decode_folder %streams_dir%\jvt\bp\*     %images_dir%\jvt\bp\
	call :decode_folder %streams_dir%\jvt\mp\*     %images_dir%\jvt\mp\
	call :decode_folder %streams_dir%\jvt\hp\*     %images_dir%\jvt\hp\
	call :decode_folder %streams_dir%\allegro\bp\* %images_dir%\allegro\bp\
	call :decode_folder %streams_dir%\allegro\mp\* %images_dir%\allegro\mp\
	call :decode_folder %streams_dir%\allegro\hp\* %images_dir%\allegro\hp\
	goto :eof
