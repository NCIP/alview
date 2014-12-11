
#"c:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\cl.exe" -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  adler32.c adler32.obj

"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  zutil.c zutil.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  adler32.c adler32.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  compress.c compress.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  crc32.c crc32.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  deflate.c deflate.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  gzclose.c gzclose.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  gzlib.c gzlib.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  gzread.c gzread.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  gzwrite.c gzwrite.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  infback.c infback.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  inffast.c inffast.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  inflate.c inflate.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  inftrees.c inftrees.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  trees.c trees.obj
"cl.exe" -MD -c -I. -D_USE_KNETFILE -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32  uncompr.c uncompr.obj

#"cl.exe" /out:zlib.lib adler32.obj compress.obj crc32.obj deflate.obj gzclose.obj gzlib.obj gzread.obj gzwrite.obj infback.obj inffast.obj inflate.obj inftrees.obj trees.obj uncompr.obj zutil.obj

lib adler32.obj compress.obj crc32.obj deflate.obj gzclose.obj gzlib.obj gzread.obj gzwrite.obj infback.obj inffast.obj inflate.obj inftrees.obj trees.obj uncompr.obj zutil.obj -OUT:zlib.lib
