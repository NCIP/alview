cl -Ox /MD -c -DWIN32=1 -DCMD_LINE=1 -I sam01832/ -I zlib32/ -I . alviewcore.cpp 
cl -Ox /MD -c -DWIN32=1 -DCMD_LINE=1 -I sam01832/ -I zlib32/ -I . alvmisc.cpp
link /MACHINE:x86  /OUT:alview_cmd_win.exe Alvmisc.obj alviewcore.obj sam01832/my.lib zlib32/zlib.lib

     
:: example:  alview_cmd_win.exe c:\rich\BAMS\RW1.BWA-aln.RG.MKDUP.bam temptest.png chr17:7512444-7519536 hg19
