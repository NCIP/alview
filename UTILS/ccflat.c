

/*
vi ccflat.c ; gcc -Wall   -o ccflat ccflat.c 
cl -Wall   -o ccflat ccflat.c 

windows ..
./ccflat refFlat.april.2014.hg19.txt hg19.refflat
./ccflat refFlat.april.2014.hg18.txt hg18.refflat

linux
... fill this in

*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

FILE *fp_refflat_d; // fixed size data
FILE *fp_refflat_o; // fixed size data

#define MAX_ACC 48

// 48 + 48 + 48 ( 8 * 4)  = 176 bytes
         /* some more info on refflat data structure is at http://genome.ucsc.edu/goldenPath/gbdDescriptions.html */
struct flattype           /* "Alignment" type */
{
    char geneName[MAX_ACC];
    char name[MAX_ACC];
    char chrom[MAX_ACC];       /* "Target sequence name" */                             /* 14 */
    unsigned int txStart;
    unsigned int txEnd;
    unsigned int cdsStart;
    unsigned int cdsEnd;
    unsigned int exonCount;
    unsigned int *exonStarts; // tricky , really a ptr to exonStart and exonEnds in another file
    unsigned int *exonEnds;
    unsigned int strand; // '+' = positive, '-' = negative   as integer !!!
};

static int numflats = 0;
static int maxflats = 0;
static struct flattype *flats; // used to be hardwired to array of MAX_FLATS, but now must be alloced.

/*
Field   Type    Null    Key     Default Extra
chrom   varchar(31)     NO      MUL
chromStart      int(10) unsigned        NO
chromEnd        int(10) unsigned        NO
name    varchar(15)     NO      MUL
score   smallint(5) unsigned    NO
strand  enum('+','-')   NO
refNCBI blob    NO
refUCSC blob    NO
observed        varchar(255)    NO
molType enum('unknown','genomic','cDNA')        NO
class   enum('single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion')     NO
valid   set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes')        NO
avHet   float   NO
avHetSE float   NO
func    set('unknown','coding-synon','intron','near-gene-3','near-gene-5','ncRNA','nonsense','missense','stop-loss','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5')NO
...
*/

struct snptype           /* "snp" type */
{
    char name[MAX_ACC];
    char chrom[MAX_ACC];       /* "Target sequence name" */                             /* 14 */
    int chromStart;
    int chromEnd;
    char refNCBI[MAX_ACC];
    char refUCSC[MAX_ACC];
    char observed[MAX_ACC];
    char func[MAX_ACC];
    double avHet;
    char *other[MAX_ACC];
};

static int numsnps = 0;
static int maxsnps = 0;
static struct snptype *snps; // used to be hardwired to array of MAX_FLATS, but now must be alloced.

int cmpfunc_snp(const void *arg1, const void *arg2)
{
    int i;
    struct snptype *z1;
    struct snptype *z2;

    z1 = (struct snptype *)arg1;
    z2 = (struct snptype *)arg2;

    i = strcmp(z1->chrom,z2->chrom);
    if (i != 0) return i;

    if (z1->chromStart < z2->chromStart) return -1;
    else if (z1->chromStart > z2->chromStart) return 1;

    if (z1->chromEnd < z2->chromEnd) return -1;
    else if (z1->chromEnd > z2->chromEnd) return 1;

    return 0;
}

#if 0
int cmpfunc_flats_by_id(const void *arg1, const void *arg2)
{
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;

    if (z1->egid < z2->egid) return -1;
    else if (z1->egid > z2->egid) return 1;
    return 0;
}


int cmpfunc_flats2(const void *arg1, const void *arg2)
{
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;

    if (z1->egid < z2->egid) return -1;
    else if (z1->egid > z2->egid) return 1;
    if (z1->uniq < z2->uniq) return -1;
    else if (z1->uniq > z2->uniq) return 1;
    return 0;
}


int cmpfunc_flats3(const void *arg1, const void *arg2)
{
    struct flattype *z1;
    struct flattype *z2;


    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;

    if (z1->khri < z2->khri) return -1;
    else if (z1->khri > z2->khri) return 1;

    if (z1->txStart < z2->txStart) return -1;
    else if (z1->txStart > z2->txStart) return 1;

    return 0;
}
#endif


int cmpfunc_refflat_by_chrstartend(const void *arg1, const void *arg2)
{
    int i;
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;

    i = strcmp(z1->chrom,z2->chrom);
    if (i != 0) return i;

    if (z1->txStart < z2->txStart) return -1;
    else if (z1->txStart > z2->txStart) return 1;

    if (z1->txEnd < z2->txEnd) return -1;
    else if (z1->txEnd > z2->txEnd) return 1;

    return 0;
}

int cmpfunc_flats5(const void *arg1, const void *arg2)
{
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;
    return strcmp(z1->chrom,z2->chrom);
}


int cmpfunc_flats_by_gene(const void *arg1, const void *arg2)
{
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;
    return strcmp(z1->geneName,z2->geneName);
}


int cmpfunc_flats_chr_start_end(const void *arg1, const void *arg2)
{
    int i;
    struct flattype *z1;
    struct flattype *z2;

    z1 = (struct flattype *)arg1;
    z2 = (struct flattype *)arg2;

    i = strcmp(z1->chrom,z2->chrom);
    if (i != 0) return i;

    if (z1->txStart < z2->txStart) return -1;
    else if (z1->txStart > z2->txStart) return 1;

    if (z1->txEnd < z2->txEnd) return -1;
    else if (z1->txEnd > z2->txEnd) return 1;

    return 0;
}

static size_t stolen_from_the_internet_strlcpy(dst, src, siz)
        char *dst;
        const char *src;
        size_t siz;
{
        register char *d = dst;
        register const char *s = src;
        register size_t n = siz;

        /* Copy as many bytes as will fit */
        if (n != 0 && --n != 0) {
                do {
                        if ((*d++ = *s++) == 0)
                                break;
                } while (--n != 0);
        }

        /* Not enough room in dst, add NUL and traverse rest of src */
        if (n == 0) {
                if (siz != 0)
                        *d = '\0';                /* NUL-terminate dst */
                while (*s++)
                        ;
        }
        return(s - src - 1);        /* count does not include NUL */
}

int chr2khri(char *chrom)
{
    int k;


    k = -1;

         if (strcmp(chrom,"chrX")==0) {  k = 22; }
    else if (strcmp(chrom,"chrY")==0) {  k = 23; }
    else if (strcmp(chrom,"chrM")==0) {  k = 24; }
    else
    {
        k = atoi(&chrom[3]);
        k = k - 1;
    }
    return k;
}


int load_ref_flats(char filename[])
{
    int numflds;
    int uniq;
    int exonCount = 0;
    register long int *liptr = (void *)0;
    register int i,j;
    register long int *z = (void *)0;
    int lineno = 0;
    int mallocsize = 0;
    FILE *fp;
    char file_name[2000];
    char s[16000];
    char t[16000];
    char tmp[16000];
static char data_input[22][5000];    /* record data from "alignment" text file gets read in (sscanf'ed) into here */


    if (flats == (void *)0) { fprintf(stderr,"ERROR: bad load flats(), no malloc.\n");  exit(0); }

    strcpy(file_name,filename);

    fp = fopen(file_name,"rb");
    if (fp == (FILE *)0)
    {
        fprintf(stderr,"Can not open file [%s] in load_ref_flats\n",file_name);
        exit(0);
    }

    numflats = uniq = 0;
    lineno = 0;
    memset(s,0,sizeof(s));
    while (fgets(s,14998,fp))
    {
        lineno++;
        for (i=0;s[i];i++) { if ((s[i] == '\n') || (s[i] == '\r') ) { s[i] = (char)0; break; } }
        memset(data_input,0,sizeof(data_input));
        numflds = sscanf(s,"%s %s %s %s %s %s %s %s %s %s %s %s",
               data_input[ 0],    /* geneName */
               data_input[ 1],    /* name */
               data_input[ 2],    /* chrom */
               data_input[ 3],    /* strand */
               data_input[ 4],    /* txStart*/
               data_input[ 5],    /* txEnd*/
               data_input[ 6],    /* cdsStart */
               data_input[ 7],    /* cdsEnd */
               data_input[ 8],    /* exonCount */
               data_input[ 9],    /* exonStarts */
               data_input[10],    /* exonEnds*/
               data_input[11]) ;  /* egid*/

        if (strstr(data_input[2],"hap") ) continue;
        if (strstr(data_input[2],"chrUn") ) continue;
        if (strstr(data_input[2],"random") ) continue;
        if (strstr(data_input[2],"chrom") ) continue;
        if (strcmp(data_input[0],"geneName") == 0) continue;

        strcpy(flats[numflats].geneName,data_input[0]);
        strcpy(flats[numflats].name,data_input[1]);
        strcpy(flats[numflats].chrom,data_input[2]);
        flats[numflats].strand = (int)data_input[3][0];
        flats[numflats].txStart  = atoi(data_input[4]);
        flats[numflats].txEnd    = atoi(data_input[5]);
        flats[numflats].cdsStart = atoi(data_input[6]);
        flats[numflats].cdsEnd   = atoi(data_input[7]);
        exonCount = flats[numflats].exonCount = atoi(data_input[8]);

/* --- */
        mallocsize = ( sizeof(long int) * exonCount );
        z = malloc( mallocsize );
        memset(z,0,mallocsize);
        if (z==(void *)0) { fprintf(stderr,"ERROR: no more space for exonCount %d %d\n",exonCount,mallocsize); exit(0); }
        stolen_from_the_internet_strlcpy(t,data_input[9],sizeof(t));  /* block sizes */
        for ( liptr=(long int *)z, i=j=0 , tmp[0] = (char)0 ; t[i] ; i++)
        {
            if (t[i] == ',' ) { *liptr++  = atol(tmp); tmp[0] = (char)0; j = 0; }
            else {tmp[j++] = t[i]; tmp[j] = (char)0; }
        }
        flats[numflats].exonStarts = (unsigned int *)z;

/* --- */
        mallocsize = ( sizeof(long int) * exonCount );
        z = malloc( mallocsize );
        memset(z,0,mallocsize);
        if (z==(void *)0) { fprintf(stderr,"ERROR: no more space for exonCount %d %d\n",exonCount,mallocsize); exit(0); }
        stolen_from_the_internet_strlcpy(t,data_input[10],sizeof(t));  /* block sizes */
        for ( liptr=(long int *)z, i=j=0 , tmp[0] = (char)0 ; t[i] ; i++)
        {
            if (t[i] == ',' ) { *liptr++  = atol(tmp); tmp[0] = (char)0; j = 0; }
            else {tmp[j++] = t[i]; tmp[j] = (char)0; }
        }
        flats[numflats].exonEnds = (unsigned int *)z;
/* --- */
#if 0
        flats[numflats].egid   = atoi(data_input[11]);
        flats[numflats].khri = chr2khri(flats[numflats].chrom);
        if (flats[numflats].khri == -1) 
        {
            fprintf(stderr,"invalid chrom = %s in reflflat input \n",data_input[2]); 
            exit(0);
        }
#endif

        if (++numflats == maxflats)
        {
            fprintf(stderr,"ERROR: very unexpected count of numflats %d == %d\n",numflats,maxflats);
            exit(0);
        }
    }
    fclose(fp);

    qsort((void *)&flats[0],numflats, sizeof(struct flattype), cmpfunc_refflat_by_chrstartend);

/* debug 
    char m[16000];
    for (i=0;(i<numflats)&&(i<1000);i++)
    {
sprintf(m,"chrom = %s:%d-%d, khri = %d ",flats[i].chrom,flats[i].txStart,flats[i].txEnd,flats[i].khri);  mydebug(m);
    }
*/

    return 0;
}


int load_snps(char filename[])
{
    register int i;
    register int numflds;
    int lineno = 0;
    FILE *fp;
    char file_name[2000];
    char s[16000];
static char data_input[22][5000];    /* record data from "alignment" text file gets read in (sscanf'ed) into here */
//     register long int *z = (void *)0;
//     int mallocsize = 0;
//    char t[16000];
//    char tmp[16000];

printf("here in load_snps(%s)\n",filename); 

    if (snps == (void *)0) { fprintf(stderr,"ERROR: bad load snps(), no malloc.\n");  exit(0); }

    strcpy(file_name,filename);

    fp = fopen(file_name,"rb");
    if (fp == (FILE *)0)
    {
        fprintf(stderr,"Can not open file [%s] in load_ref_snps\n",file_name);
        exit(0);
    }

    numsnps = 0;
    lineno = 0;
    memset(s,0,sizeof(s));
    while (fgets(s,14998,fp))
    {
        lineno++;
        for (i=0;s[i];i++) { if ((s[i] == '\n') || (s[i] == '\r') ) { s[i] = (char)0; break; } }
        memset(data_input,0,sizeof(data_input));

        numflds = sscanf(s,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ",
               data_input[ 0], // bin
               data_input[ 1], // chrom
               data_input[ 2], // chromstart
               data_input[ 3], // chromEned
               data_input[ 4], // name
               data_input[ 5], // score
               data_input[ 6], // strand
               data_input[ 7], // refncbi
               data_input[ 8], // refucsc
               data_input[ 9], // observed
               data_input[10], // moltype
               data_input[11], // class
               data_input[12], // valid
               data_input[13], // avHet
               data_input[14], // mavHetSe
               data_input[15], // func
               data_input[16]);

        if (strstr(data_input[1],"hap") ) continue;
        if (strstr(data_input[1],"chrUn") ) continue;
        if (strstr(data_input[1],"random") ) continue;
        if (strstr(data_input[1],"chrom") ) continue; // don't want header

/* yyy 
    char name[MAX_ACC];
    char chrom[MAX_CHRNAME];
    char chromStart;
    char chromEnd;
    char refNCBI[MAX_ACC];
    char refUCSC[MAX_ACC];
    char observed[MAX_ACC];
    char func[MAX_ACC];
    char other[MAX_ACC];
*/
        strcpy(snps[numsnps].chrom,data_input[1]);
        snps[numsnps].chromStart  = atoi(data_input[2]);
        snps[numsnps].chromEnd    = atoi(data_input[3]);
        strncpy(snps[numsnps].name,      data_input[4],MAX_ACC);
        strncpy(snps[numsnps].refNCBI, data_input[7],MAX_ACC);
        strncpy(snps[numsnps].refUCSC, data_input[8],MAX_ACC);
        strncpy(snps[numsnps].observed,data_input[9],MAX_ACC);
        strncpy(snps[numsnps].func,data_input[15],MAX_ACC);
        snps[numsnps].avHet = atof(data_input[13]);

// printf("%d %s\n",numsnps,snps[numsnps].name);

        if (numsnps++ >= maxsnps)
        {
            fprintf(stderr,"ERROR: very unexpected count of numsnps %d (%d)\n",numsnps,maxsnps);
            exit(0);
        }
    }
    fclose(fp);

    qsort((void *)&snps[0],numsnps, sizeof(struct snptype), cmpfunc_snp);

#if 0
    qsort((void *)&snps[0],numsnps, sizeof(struct snptype), cmpfunc_snps4);
    prev = -1;
    for (i=0;i<numsnps;i++)
    {
        if (prev != snps[i].egid)
           uniq = 1;
        snps[i].uniq = uniq++;
        prev = snps[i].egid;
    }

/* debug 
    char m[16000];
    for (i=0;(i<numflats)&&(i<1000);i++)
    {
sprintf(m,"chrom = %s:%d-%d, khri = %d ",flats[i].chrom,flats[i].txStart,flats[i].txEnd,flats[i].khri);  mydebug(m);
    }
*/
#endif

printf("end load_snps()\n"); 
    return 0;
}



int is_in_exon(FILE *fp, long int spot, int cnt, int loc)
{
    int i = 0;
    int stat = 0;
    long int start = 0;
    long int end = 0;

    fseek(fp, spot, SEEK_SET);
    for (i=0;i<cnt;i++)
    {
        stat = fread(&start,sizeof(long int),1,fp); if (stat != 1) break;
        stat = fread(&end,sizeof(long int),1,fp); if (stat != 1) break;
// printf("debug in is_in_exon, spot = %ld , cnt=%d loc=%d start=%ld end=%ld\n",spot,cnt,loc,start,end);
        if ((loc >= start) && (loc < end)) return 1;
    }
    return 0;
}


   

int count_lines_in_file(char *txt_file_name)
{
    int cnt = 0; // count
    char s[12000];
    FILE *fp;
    fp = fopen(txt_file_name,"r");
    if (!fp) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order to count the number of lines [ in count_lines_in_file() ] \n",txt_file_name); exit(0); }
    while (fgets(s,11999,fp)) cnt++;
    fclose(fp);
    return cnt;
}


               // used for annotation
void compile_refflat(char *input_txt_refflat_file_name, char *outname)
{
    size_t spot;
    char output_file_name_flat_d[512];
    char output_file_name_flat_o[512];
    int i,j;
    FILE *fp;
    int cnt = 0;           // count

    cnt = count_lines_in_file(input_txt_refflat_file_name);
    if (cnt < 10) { fprintf(stderr,"ERROR: lines of refflat file named \"%s\" is too low, probably it's a bad file.  %d ?\n",input_txt_refflat_file_name,cnt); exit(0); }
    fp = fopen(input_txt_refflat_file_name,"r");
    if (!fp) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order read refflat data [ in compile_refflat() ] \n",input_txt_refflat_file_name); exit(0); }

    maxflats = cnt;
    flats = (void *)malloc((cnt + 1) * sizeof(struct flattype));
    memset(flats,0,((cnt + 1) * sizeof(struct flattype)));
    load_ref_flats(input_txt_refflat_file_name);
// flats now loaded and sorted

    sprintf(output_file_name_flat_d,"%s.dat",outname);
    sprintf(output_file_name_flat_o,"%s.dao",outname);

    fp_refflat_d = fopen(output_file_name_flat_d,"wb");
    if (!fp_refflat_d) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order refflat data file[in compile_refflat() 2]\n",output_file_name_flat_d); exit(0); }

    fp_refflat_o = fopen(output_file_name_flat_o,"wb");
    if (!fp_refflat_o) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order open refflat exon overflow file [ in compile_refflat() 3]3 \n",output_file_name_flat_o); exit(0); }

    for (i=0 ; i<numflats ; i++)
    {
        spot = ftell(fp_refflat_o);
        for (j=0 ; j<flats[i].exonCount ; j++)
        {
            fwrite(&flats[i].exonStarts[j],sizeof(unsigned int),1,fp_refflat_o);
            fwrite(&flats[i].exonEnds[j],sizeof(unsigned int),1,fp_refflat_o);
        }
        flats[i].exonStarts = (unsigned int *)spot;                 // HACK OVER RIDE  - this is now pointer into file
// sprintf(m,"chrom = %s:%d-%d, khri = %d ",flats[i].chrom,flats[i].txStart,flats[i].txEnd,flats[i].khri);  mydebug(m);
        fwrite (&flats[i].geneName,MAX_ACC,1,fp_refflat_d);
        fwrite (&flats[i].name,MAX_ACC,1,fp_refflat_d);
        fwrite (&flats[i].chrom,MAX_ACC,1,fp_refflat_d);
        fwrite (&flats[i].txStart,4,1,fp_refflat_d);
        fwrite (&flats[i].txEnd,4,1,fp_refflat_d);
        fwrite (&flats[i].cdsStart,4,1,fp_refflat_d);
        fwrite (&flats[i].cdsEnd,4,1,fp_refflat_d);
        fwrite (&flats[i].exonCount,4,1,fp_refflat_d);
        fwrite (&flats[i].exonStarts,4,1,fp_refflat_d);
        fwrite (&flats[i].exonEnds,4,1,fp_refflat_d);
        fwrite (&flats[i].strand,4,1,fp_refflat_d);
// xxx
    }

    fclose(fp_refflat_d);
    fclose(fp_refflat_o);
    if (flats) free(flats);

    return;
}


#if 0
void compile_snp(char *input_txt_snp_file_name, char *outname)
{
    char output_file_name_snp_d[512];
    int i;
    FILE *fp;
    int cnt = 0;           // count

printf("here in compile_snp()\n"); 
    cnt = count_lines_in_file(input_txt_snp_file_name);
    if (cnt < 10) { fprintf(stderr,"ERROR: lines of snp file named \"%s\" is too low, probably it's a bad file.  %d ?\n",input_txt_snp_file_name,cnt); exit(0); }
    fp = fopen(input_txt_snp_file_name,"r");
    if (!fp) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order read snp data [ in compile_snp() ] \n",input_txt_snp_file_name); exit(0); }

    maxsnps = cnt;
printf("malloc %ld for snps ,cnt = %d\n", (cnt + 1) * sizeof(struct snptype),cnt);
    snps = (void *)malloc((cnt + 1) * sizeof(struct snptype));
printf("malloced  %ld for snps \n", (cnt + 1) * sizeof(struct snptype));
    memset(snps,0,((cnt + 1) * sizeof(struct snptype)));
printf("memset  %ld for snps \n", (cnt + 1) * sizeof(struct snptype));
    load_snps(input_txt_snp_file_name);

    sprintf(output_file_name_snp_d,"%s.dat",outname);
printf("snp loaded , file name is  %s (for output)\n", output_file_name_snp_d);

    fp_snp_d = fopen(output_file_name_snp_d,"w");
    if (!fp_snp_d) { fprintf(stderr,"ERROR: can't open file named \"%s\" \n",output_file_name_snp_d); exit(0); }

    for (i=0 ; i<numsnps ; i++)
    {
if ((i % 100000)==0) { printf("%d %s %s %s \n",i,snps[i].chrom,snps[i].name,snps[i].observed); fflush(stdout); }
        fwrite (&snps[i],sizeof(struct snptype),1,fp_snp_d);
    }
    fclose(fp_snp_d);

    return;
}
#endif


#define SIZE_REFFLAT_REC sizeof(struct flattype)
long int last_refflat_fseek_place;
struct flattype temp_refflat_space;

struct flattype *binary_search_refflat_file(struct flattype *match_me, long int lo, long int hi)
{
    long int seekto;
    long int mid;
    int k;

    if (fp_refflat_d == (void *)0) { fprintf(stderr,"ERROR: in binary_search_refflat_file() - null fp_refflat_d\n"); exit(0); } 
// printf("dbg in binary_search_refflat_file() chrom=%s loc=%d last_refflat_fseek_place=%ld lo=%ld hi=%ld\n",match_me->chrom,match_me->txStart,last_refflat_fseek_place,lo,hi); 

    if (lo>hi) 
    {
//        fprintf(stderr,"ERROR: can't find cuz lo > hi %ld and %ld in binary_search_refflat_file\n",lo,hi);
        return (void *)0;
    }
    if (lo==hi)
    {
        seekto = (SIZE_REFFLAT_REC * lo);
        if (fseek(fp_refflat_d,(size_t) seekto , SEEK_SET) < 0)
        {
            fprintf(stderr,"SEEK to %ld ERROR in binary_search_refflat_file lo=%ld\n",seekto,lo);
            exit(0);
        }
        last_refflat_fseek_place = seekto; 
        fread(&temp_refflat_space,1,SIZE_REFFLAT_REC,fp_refflat_d);
        if (cmpfunc_flats_chr_start_end(&temp_refflat_space,match_me) == 0) 
        {
// fprintf(stderr,"GOT binary_search_refflat_file lo=%ld hi=%ld \n",lo,hi); fflush(stderr); 
            return &temp_refflat_space;
        }
        else 
        {
            return (void *)0;
        }
    }
    mid = (lo + hi) / 2 ;

    seekto = (SIZE_REFFLAT_REC * mid);
//printf("in search seekto=%ld [ which is mid=%d]\n",seekto,mid); fflush(stdout);

    if (fseek(fp_refflat_d,(size_t) seekto , SEEK_SET) < 0)
    {
        fprintf(stderr,"SEEK to %ld ERROR in binary_search_refflat_file, mid=%ld\n",seekto,mid);
        exit(0);
    }

    last_refflat_fseek_place = seekto; 
    fread( &temp_refflat_space, 1 , SIZE_REFFLAT_REC , fp_refflat_d );

// printf("dbg read in %s:%d-%d at last_refflat_fseek_place =%ld \n",temp_refflat_space.chrom,temp_refflat_space.txStart,temp_refflat_space.txEnd,last_refflat_fseek_place);
    k = cmpfunc_flats_chr_start_end(&temp_refflat_space,match_me);
// fprintf(stderr,"after cmpfunc_flats_chr_start_end() k=%d \n",k);
    if (k == 0)  // unlikely !!!
    {
// printf("in search found , mid = %ld\n",mid); fflush(stdout);
        return &temp_refflat_space;
    }
    if (strcmp(temp_refflat_space.chrom,match_me->chrom) == 0) 
    {
        if ((match_me->txStart >= temp_refflat_space.txStart)  && (match_me->txEnd < temp_refflat_space.txEnd))
            return &temp_refflat_space;
    }
    if (k > 0) 
    {
        return binary_search_refflat_file(match_me, lo, mid-1);
    }
    else 
    {
        return binary_search_refflat_file(match_me, mid+1, hi);
    }
}


long int refflat_fixed_hi;
void setup_refflat(char *input_flat_file_name_d, char *input_flat_file_name_o)
{
    fp_refflat_d = fopen(input_flat_file_name_d,"r");
    if (fp_refflat_d == (void *)0) { fprintf(stderr,"ERROR: %s\n",input_flat_file_name_d); fflush(stderr); exit(0); }

printf("in setup_refflat FILE=%p \n",fp_refflat_d); 
    fseek(fp_refflat_d,(size_t) 0 , SEEK_END);
    refflat_fixed_hi = ftell(fp_refflat_d) / SIZE_REFFLAT_REC;
printf("in setup_refflat refflat_fixed_hi=%ld \n",refflat_fixed_hi);  fflush(stdout); 

    fp_refflat_o = fopen(input_flat_file_name_o,"r");
    if (!fp_refflat_o) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order setup in setup_refflat() ] \n",input_flat_file_name_o); exit(0); }
printf("end  setup_refflat refflat_fixed_hi=%ld \n",refflat_fixed_hi);  fflush(stdout); 
}

void close_refflat(void)
{
    fclose(fp_refflat_d);
    fclose(fp_refflat_o);
}


#define SIZE_SNPTYPE sizeof(struct snptype)
long int last_snp_fseek_place;
struct snptype temp_snp_space;


#if 0
struct snptype *binary_search_snp_file(struct snptype *match_me, long int lo, long int hi)
{
    long int seekto;
    long int mid;
    int k;


    if (fp_snp_d == (void *)0) { fprintf(stderr,"ERROR: in binary_search_snp_file() - null fp_snp_d\n"); exit(0); } 
// printf("dbg in binary_search_snp_file() chrom=%s loc=%d last_snp_fseek_place=%ld lo=%ld hi=%ld\n",match_me->chrom,match_me->chromStart,last_snp_fseek_place,lo,hi); 

    if (lo > hi) 
    {
//        fprintf(stderr,"ERROR: can't find [%s] cuz lo > hi %ld and %ld in binary_search_snp_file\n",match_me,lo,hi);
        return (void *)0;
    }
    if (lo==hi)
    {
        seekto = (SIZE_SNPTYPE * lo);
        if (fseek(fp_snp_d,(size_t) seekto , SEEK_SET) < 0)
        {
            fprintf(stderr,"SEEK to %ld ERROR in binary_search_snp_file lo=%ld\n",seekto,lo);
            exit(0);
        }

        last_snp_fseek_place = seekto; 
        fread(&temp_snp_space,1,SIZE_SNPTYPE,fp_snp_d);
        if (cmpfunc_snp(&temp_snp_space,match_me) == 0) 
        {
// fprintf(stderr,"GOT binary_search_snp_file lo=%ld hi=%ld %s\n",lo,hi,temp_snp_space); fflush(stderr); 
            return &temp_snp_space;
        }
        else 
        {
            return (void *)0;
        }
    }
    mid = (lo + hi) / 2 ;

    seekto = (SIZE_SNPTYPE * mid);

    if (fseek(fp_snp_d,(size_t) seekto , SEEK_SET) < 0)
    {
        fprintf(stderr,"SEEK to %ld ERROR in binary_search_snp_file, mid=%ld\n",seekto,mid);
        exit(0);
    }
    last_snp_fseek_place = seekto; 
    fread( &temp_snp_space, 1 , SIZE_SNPTYPE , fp_snp_d );

    k = cmpfunc_snp(&temp_snp_space,match_me);
    if (k == 0)  // unlikely !!!
    {
        return &temp_snp_space;
    }
    if (strcmp(temp_snp_space.chrom,match_me->chrom) == 0) 
    {
        if ((match_me->chromStart >= temp_snp_space.chromEnd) && (match_me->chromEnd < temp_snp_space.chromEnd))
            return &temp_snp_space;
    }
    if (k > 0) 
    {
        return binary_search_snp_file(match_me, lo, mid-1);
    }
    else 
    {
        return binary_search_snp_file(match_me, mid+1, hi);
    }
}
#endif


#if 0
static long int snp_fixed_hi;
static void setup_snp(char *input_flat_file_name_d)
{
printf("in setup_snp %s\n",input_flat_file_name_d);
    fp_snp_d = fopen(input_flat_file_name_d,"r");
    if (fp_snp_d == (void *)0) { fprintf(stderr,"ERROR: file name \"%s\" in setup_snp\n",input_flat_file_name_d); fflush(stderr); exit(0); }
    fseek(fp_snp_d,(size_t) 0 , SEEK_END);
    snp_fixed_hi = ftell(fp_snp_d)/SIZE_SNPTYPE;

/*
    fp_snp_o = fopen(input_flat_file_name_o,"r");
    if (!fp_snp_o) { fprintf(stderr,"ERROR: can't open file named \"%s\" in order setup in setup_snp() ] \n",input_flat_file_name_o); exit(0); }
*/
}

static void close_snp(void)
{
    fclose(fp_snp_d);
}

void snp_annot( char khr[],int loc, char stuff[])
{
    long int spot;
    struct snptype f;
    struct snptype *z;


    strcpy(stuff,"unknown_snp");

    memset(&f,0,sizeof(struct snptype));
    strcpy(f.chrom,khr);
    f.chromStart = loc -1;
    f.chromEnd = loc;

// printf("in snp_annot %s:%d\n",khr,loc); 
//
    z = binary_search_snp_file(&f,0,snp_fixed_hi);
    if (z)
    {
// loc = 212723198 , chromStart = 212723197 
// fprintf(stderr,"loc = %d , chromStart = %d \n", loc,z->chromStart); 
        if ((strcmp(z->chrom,khr) == 0) && ((loc-1) == z->chromStart))
        {
sprintf(stuff,"%s %s %s %s %s %10.5f",z->refNCBI,z->refUCSC,z->observed,z->name,z->func,z->avHet);
            return;
        }
    }

    spot = last_snp_fseek_place;

#if 0
    int stat;
if (!z) 
{
    fseek(fp_snp_d,spot,SEEK_SET);
    stat = fread(&temp_snp_space,1,SIZE_SNPTYPE ,fp_snp_d);
fprintf(stderr,"MISSED %s %s %s %s %s %10.5f , spot = %ld\n",temp_snp_space.refNCBI,temp_snp_space.refUCSC,temp_snp_space.observed,temp_snp_space.name,temp_snp_space.func,temp_snp_space.avHet,spot);
}
// else fprintf(stderr,"MISSED loc=%d %s %s %s %s %s %10.5f",loc,z->refNCBI,z->refUCSC,z->observed,z->name,z->func,z->avHet);
    return;
#endif

#if 0
    if (spot < 0) spot = 0;
// printf("in snp_annot() after binary_search_snp_file %s %d spot=%ld\n",khr,loc,spot); 
    j = 0;
    while ((j<100)) 
    {
        stat = fread(&temp_snp_space,1,SIZE_SNP_REC,fp_snp_d);
        if (stat < 1) break;
// printf("wantchrloc=%s:%d  gotchr=%s gene?=%s gotloc=%d\n",khr,loc,temp_snp_space.chrom,temp_snp_space.geneName,temp_snp_space.txStart);
        if (strcmp(temp_snp_space.chrom,khr) == 0) 
        {
            if ((loc >= temp_snp_space.txStart) && (loc < temp_snp_space.txEnd)) 
            {
                strcpy(stuff,temp_snp_space.geneName);
                *ingene = 1;
// exonStarts is REALLY the seek spot for exonStart/exonEnds in overflow file
                *inexon = is_in_exon(fp_snp_o,(long int)(temp_snp_space.exonStarts),temp_snp_space.exonCount, temp_snp_space.txStart);
                if (*inexon) return; // kick out if got exon
            }
            if ((temp_snp_space.txStart > loc) && (temp_snp_space.txEnd > loc)) 
                break;
        }
        j++;
    }
    fseek(fp_snp_d,0,SEEK_SET);
#endif
    return;
}
#endif


long int refflat_fixed_hi;


#if 0
void gene_annot( char khr[],int loc, int *ingene, int *inexon, char *name) 
{
    int stat;
    long int spot;
    int j = 0;
    struct flattype f;
    struct flattype *z;

    *ingene = *inexon = 0;
    strcpy(name,"UNKNOWN");

    memset(&f,0,sizeof(struct flattype));
    strcpy(f.chrom,khr); 
    f.txStart = f.txEnd = loc+1;

// printf("in gene_annot %s:%d %d %ld\n",khr,loc,0,refflat_fixed_hi); 
    z = binary_search_refflat_file(&f,0,refflat_fixed_hi);
// not really expecting a match, but need last looked at location in the binary refflat file
// just badk up 100 and look for 200

    spot = last_refflat_fseek_place - (SIZE_REFFLAT_REC*50);
    if (spot < 0) spot = 0;
// printf("in gene_annot() after binary_search_refflat_file %s %d spot=%ld\n",khr,loc,spot); 
    fseek(fp_refflat_d,last_refflat_fseek_place,SEEK_SET);
    j = 0;
    while ((j<100)) 
    {
        stat = fread(&temp_refflat_space,1,SIZE_REFFLAT_REC,fp_refflat_d);
        if (stat < 1) break;
// printf("wantchrloc=%s:%d  gotchr=%s gene?=%s gotloc=%d\n",khr,loc,temp_refflat_space.chrom,temp_refflat_space.geneName,temp_refflat_space.txStart);
        if (strcmp(temp_refflat_space.chrom,khr) == 0) 
        {
            if ((loc >= temp_refflat_space.txStart) && (loc < temp_refflat_space.txEnd)) 
            {
                strcpy(name,temp_refflat_space.geneName);
                *ingene = 1;
// exonStarts is REALLY the seek spot for exonStart/exonEnds in overflow file
                *inexon = is_in_exon(fp_refflat_o,(long int)(temp_refflat_space.exonStarts),temp_refflat_space.exonCount, temp_refflat_space.txStart);
                if (*inexon) return; // kick out if got exon
            }
            if ((temp_refflat_space.txStart > loc) && (temp_refflat_space.txEnd > loc)) 
                break;
        }
        j++;
    }
    fseek(fp_refflat_d,0,SEEK_SET);
    return;
}
#endif


void usage(int argc)
{
    fprintf(stderr,"./somad COMPILE REFFLAT ../refFlat hg19.refflat\n"); 
}


int main(int argc, char *argv[])
{
    fprintf(stderr,"example : ./ccflat ../refFlat hg19.refflat\n"); 
    fprintf(stderr,"recsize = %d \n",sizeof(struct flattype)); 
    compile_refflat(argv[1], argv[2]); 
    return 0;
}


