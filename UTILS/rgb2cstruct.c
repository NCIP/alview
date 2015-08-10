
#if 0
vi rgb2cstruct.c ; gcc -Wall   -o rgb2cstruct rgb2cstruct.c 
 utility to convert fro rgb to cstruct to put image in static memory (i.e. a "C" array )
 usage : cat tmp.rgb | ./rgb2cstruct > outfile
#endif

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>


int main(int argc, char *argv[])
{
    int i = 0;
    int c;

    printf("unsigned char namehere[] = {\n"); 
    while ((c = getchar())  !=EOF)
    {
        printf("0x%02x,",c&0xFF); 
        i++;
        if ((i%15) == 0) printf("\n"); 
// unsigned char ncilogo[] = {
// 0xa9,0x01,0x01,0xa9,0x01,0x01,0xa9,0x01,0x01,0xa9,0x01,0x01,0xa9,0x01,0x01,0xa9,0x01,0x01,
    }
    printf("};\n"); 
fprintf(stderr,"%d\n",i); 
    return 0;
}

