
rc alview.rc
cl -Ox /MD -c -I . StdAfx.cpp
copy alvmisc.c alvmisc.cpp
cl -Ox /EHsc /MD -c -DWIN32=1 -UUNIX -I sam01832/ -I zlib32/ -I . alvmisc.cpp 
cl -Ox /MD -c -DWIN32=1 -UUNIX -I sam01832/ -I zlib32/ -I . alviewcore.cpp 
cl -Ox /MD -c -I . alvwin32.cpp 
:: cl /Wall /MD -c -I . alvwin32.cpp 

:: Windows 10 now needs shell32.lib ?
:: next line works with vs2012 32 bit
::link /MACHINE:x86 alview.res StdAfx.obj alvmisc.obj alviewcore.obj alvwin32.obj user32.lib gdi32.lib comdlg32.lib sam01832/my.lib zlib32/zlib.lib

:: next line works with Windows 10 
link /MACHINE:x86 alview.res StdAfx.obj alvmisc.obj alviewcore.obj alvwin32.obj user32.lib gdi32.lib comdlg32.lib shell32.lib sam01832/my.lib zlib32/zlib.lib


