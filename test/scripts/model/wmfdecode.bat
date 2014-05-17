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
@rem  File      : wmfdecode.bat
@rem  Author(s) : Luuvish
@rem  Version   : 1.0
@rem  Revision  :
@rem      1.0 May 20, 2013    first release
@rem
@rem ===========================================================================


set streams_dir=..\..\streams\vc1
set images_dir=..\..\images\vc1
set log_file=microsoft.log

echo. > %log_file%

goto :decode


:decode_file
	call wmfdecode_vc.exe -v IYUV_WMV -cut1stframe 1 -ignoreaudio 1 -i %1 -o %~dpn2.yuv
	goto :eof


:decode_echo
	echo wmfdecode_vc -i %~nx1 -o %~n2.yuv
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
	call :decode_folder %streams_dir%\artificial6\* %images_dir%\artificial6\
	call :decode_folder %streams_dir%\artificial7\* %images_dir%\artificial7\
	call :decode_folder %streams_dir%\smpte-vc1\*   %images_dir%\smpte-vc1\
	call :decode_folder %streams_dir%\conf-vc1\*    %images_dir%\conf-vc1\
	call :decode_folder %streams_dir%\conf-wmv9\*   %images_dir%\conf-wmv9\
	call :decode_folder %streams_dir%\stress\*      %images_dir%\stress\
	goto :eof
