

/*
vi +44 flatcheck.c ; gcc -Wall -o flatcheck flatcheck.c
cl -Wall -o flatcheck flatcheck.c 

./flatcheck hg18.refflat.dat hg18.refflat.dao
./flatcheck hg19.refflat.dat hg19.refflat.dao

./flatcheck.exe  ../release/GENOMEDATA/hg18/hg18.refflat.dat ../release/GENOMEDATA/hg18/hg18.refflat.da0

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>



#define MAX_ACC 48
#define MAX_CHRNAME 48

         /* some more info on refflat data structure is at http://genome.ucsc.edu/goldenPath/gbdDescriptions.html */
struct flattype           /* "Alignment" type */
{
    char geneName[MAX_ACC];
    char name[MAX_ACC];
    char chrom[MAX_CHRNAME];       /* "Target sequence name" */                             /* 14 */
    unsigned int txStart;
    unsigned int txEnd;
    unsigned int cdsStart;
    unsigned int cdsEnd;
    unsigned int exonCount;
    unsigned int *exonStarts; // tricky , really a ptr to exonStart and and exonEnds in another file
    unsigned int *exonEnds;
    unsigned int strand; // '+' = positive, '-' = negative   as integer !!!
          // just want a nice alignment
    unsigned int dork; // to give it 8 byte alignment (just add a dead field of 4 bytes)
};


int main(int argc , char *argv[])
{
    unsigned int start;
    unsigned int end;
    size_t read_status;
    int jj;
    int cnt;
    struct flattype r;
    FILE *fp;
    FILE *fp2;

    cnt = 0;
    fp = fopen(argv[1],"rb"); 
    fp2 = fopen(argv[2],"rb"); 
fprintf(stderr,"fopened %s %s \n", argv[1],argv[2]); fflush(stderr);
    while (1)
    {
        memset(&r,0,sizeof(struct flattype));
        if (fread (&r.geneName,MAX_ACC,1,fp) != 1) break;
        if (fread (&r.name,MAX_ACC,1,fp) != 1 ) break;
        if (fread (&r.chrom,MAX_ACC,1,fp) != 1 ) break;
        if (fread (&r.txStart,4,1,fp) != 1 ) break;
        if (fread (&r.txEnd,4,1,fp) != 1 ) break;
        if (fread (&r.cdsStart,4,1,fp) != 1 ) break;
        if (fread (&r.cdsEnd,4,1,fp) != 1 ) break;
        if (fread (&r.exonCount,4,1,fp) != 1 ) break;
        if (fread (&r.exonStarts,4,1,fp) != 1 ) break;
        if (fread (&r.exonEnds,4,1,fp) != 1 ) break;
        if (fread (&r.strand,4,1,fp) != 1 ) break;
// xxx
printf("szof=%ld |%s| %s %u %u xons=%d\n",sizeof(struct flattype),r.geneName,r.chrom,r.txStart,r.txEnd,r.exonCount); 
        fseek(fp2, (long int)(r.exonStarts), SEEK_SET); // really an offset into an "overflow" file
        for (jj=0 ; jj<r.exonCount ; jj++) 
        {
            end = start = (long int) 0;
            read_status = fread(&start,sizeof(unsigned int),1,fp2); if (read_status != 1) break;
            read_status = fread(&end,sizeof(unsigned int),1,fp2); if (read_status != 1) break;
printf("%u %u\n",start,end);
        }
printf("\n"); 
        cnt++;
    }
    fclose(fp); 
    return 0;
}

