:: MD Applications compiled with this option are statically linked to MSVCRT.lib _MT and causes the compiler to place the library name LIBCMT.lib 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o bcf.obj  bcftools/bcf.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o bcf2qcall.obj  bcftools/bcf2qcall.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o bcfutils.obj  bcftools/bcfutils.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o call1.obj  bcftools/call1.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o em.obj  bcftools/em.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o fet.obj  bcftools/fet.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o index.obj  bcftools/index.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o kfunc.obj  bcftools/kfunc.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o kmin.obj  bcftools/kmin.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o main.obj  bcftools/main.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o mut.obj  bcftools/mut.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o prob1.obj  bcftools/prob1.c 
cl -c -I. -D_CURSES_LIB=2 -D_MSC_VER -I. -Iwin32 /MD -o vcf.obj  bcftools/vcf.c 
exit /B
