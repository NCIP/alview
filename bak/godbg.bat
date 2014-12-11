

:: cl /EHsc -c -I c:/rich/H/ app_four.c
::odbc32.lib odbccp32.lib 

rc alview.rc
cl /Zi /MD -c -I . StdAfx.cpp
cp alvmisc.c alvmisc.cpp
cl /Zi /EHsc /MD -c -DWIN32=1 -UUNIX -I sam01832/ -I zlib32/ -I . alvmisc.cpp 
cl /Zi /MD -c -DWIN32=1 -UUNIX -I sam01832/ -I zlib32/ -I . alviewcore.cpp 
cl /Zi /MD -c -I . alvwin32.cpp 
:: cl /Wall /MD -c -I . alvwin32.cpp 
link /DEBUG /pdb:alview.pdb alview.res StdAfx.obj alvmisc.obj alviewcore.obj alvwin32.obj user32.lib gdi32.lib comdlg32.lib sam01832/my.lib zlib32/zlib.lib


:: -L/usr/lib -L/home/rfinney/samtools-0.1.18/ -lbam -lm -lz 

