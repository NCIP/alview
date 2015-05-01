

/* 
Written by Richard Finney, 2015, National Cancer Institute, National Institutes of Health 
This file "alviewcore.cpp" contains the core alview routines.  
This is *** MOSTLY ***  public domain ...  EXCEPTION: note GOTOH code by Peter Clote.  
Clote's GOTOH (optimized Smith-Watterman) code is only for non-commercial users.  
To deal with this : just keep keep #define GOTOH=0 to disable it OR re-engineer it if you are commercial outfit.
Otherwise, do whatever you want!  NIH assumes no liabilities or responsibilities.  Use at your own risk.

#Compile for Webserver on Linux ...
g++ -DUNIX -DWEB_SERVER=1 -UQT_GUI -UCOMMAND_LINE -DNATIVE -Wall -o alview alviewcore.cpp alvmisc.cpp -lz -Wall -I/h1/finneyr/samtools-0.1.18/ /h1/finneyr/samtools-0.1.18/libbam.a -lm -lz

# notes ...
#Target specific options:
#  -m64                      Generate 64bit x86-64 code
#  -minline-all-stringops    Inline all known string operations
#  -mno-align-stringops      Do not align destination of the string operations
#  -momit-leaf-frame-pointer Omit the frame pointer in leaf functions


Compile for command line Using NEW samtools 1.0 ( after august 2014)  on linux 
gcc -Wall -DSAMTOOLS1=1 -DUNIX=1 -DCMD_LINE=1 -o alview alviewcore.cpp  -I. \
-I/h1/finneyr/samtools-1.0/ \
-I/h1/finneyr/samtools-1.0/htslib-1.0/ \
-I/h1/finneyr/samtools-1.0/htslib-1.0/htslib/ \
-L/h1/finneyr/samtools-1.0/ \
/h1/finneyr/samtools-1.0/libbam.a \
/h1/finneyr/samtools-1.0//htslib-1.0/libhts.a \
-lbam -lm -lz -lgd -lpthread -lstdc++

Compile for command line Using old samtools  - compile for command line on linux ...
****** note, can't mix c and cpp effortlessley
cp alvmisc.c alvmisc.cpp
gcc -Wall -DUNIX=1 -DCMD_LINE=1 -o alview_cmdline_linux alvmisc.cpp alviewcore.cpp -I. -I/data/nextgen/finneyr/samtools-0.1.18/ -L/data/nextgen/finneyr/samtools-0.1.18/  -lbam -lm -lz  -lstdc++

Compile for command line cygwin (USE SAMTOOLS 0.18) 
cp alvmisc.c alvmisc.cpp
gcc -Wall -DUNIX=1 -DCMD_LINE=1  -o alview_cmdline_cygwin alvmisc.cpp alviewcore.cpp -I.  -Ic:/rich/samtools-0.1.18 -Lc:/rich/samtools-0.1.18 -lbam -lm -lz  -lstdc++ 

Example command line run (cygwin) :
G/alview_cmdline_cygwin.exe  c:/rich/BAMS/RW1.BWA-aln.RG.MKDUP.bam x.png chr17:7512444-7519536 hg19

Compile command line for Windows cl.exe compliter (visual C), note usage of hacked samtools library  
cl -Ox /MD -c -DWIN32=1 -DCMD_LINE=1 -I sam01832/ -I zlib32/ -I . alviewcore.cpp 
cl -Ox /MD -c -DWIN32=1 -DCMD_LINE=1 -I sam01832/ -I zlib32/ -I . alvmisc.cpp
link /MACHINE:x86  /OUT:alview_cmd_win.exe Alvmisc.obj alviewcore.obj sam01832/my.lib zlib32/zlib.lib

     




todo: (soon !) 
Prev Next
guard against too big mallocs
multi image operations
divorce from hg18, hg19

fix these messages 
[bam_parse_region] fail to determine the sequence name.

example web usage in Excel ...
=HYPERLINK("https://cgwb-test.nci.nih.gov:8443/cgi-bin/alview?position="&D1&":"&E1-100&"-"&E1+100&"&iw=1000&ih=400&file="&A1,"n")
*/

/*
***************** ---------- define these on compile line , example gcc -DWEB_SERVER=1
#define UNIX 0
#define WEB_SERVER 0
#define USE_MYSQL 0
#define CMD_LINE 0
#define QT_GUI 0
*/


#define GOTOH 0
// gotoh turns on Peter Clotes GOTOH code   It is limited to noncommercial use.  You can re-engineer it for commercial use if you really need it.
// You may TURN ON GOTOH if you are non-commerical 

// WxWidgets version , possible someday but not implemented for now ...
#define WX_GUI 0

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#if 0 // def WIN32 // i think some version of msvc didn't get this right
#define unlink _unlink
#endif


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "alview.h" 

#ifdef WIN32
// windows wacky-ness 
#include "WINWACKY/mman.h"
#include <time.h>
#include <WINWACKY/dirent.h>
#include <io.h> 
#include <sys/stat.h>
#include <sys/types.h>
#endif

#ifdef UNIX
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#endif

#if MAC
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <dirent.h>
#endif

#ifdef QT_GUI
#ifdef WIN32

#if 0
*** boutell gd graphics library had been removed ***
// needed for gd.h stuff ...
#define NONDLL 1
#define BGDWIN32 1
#include "../libgd-2.1.0_64/rpf/gd.h"
#include "../libgd-2.1.0_64/rpf/gdfontl.h"
#endif

#else
#undef WIN32 
#include <gd.h>
#include <gdfontl.h>
#endif
#include "interface.h" // routines to make work on a GUI
#else

#define USE_GD 0 // got rid of gd graphics library - using custom implementation now
#if USE_GD
#include <gd.h> // gd library
gdImagePtr im;
#else

#endif

#define MAXBUFF 20000
#define MAXSMALLBUFF 2000

#ifdef CMD_LINE
int snp_call_flag = 0;
int snp_call_spot = 0; // genomic location 
int snp_call_dnaat = 0;
char snp_call_referece = ' ';
char snp_call_chr[MAXBUFF];
int snp_call_A_cnt = 0;
int snp_call_C_cnt = 0;
int snp_call_G_cnt = 0;
int snp_call_T_cnt = 0;
int snp_call_Ins_cnt = 0;
int snp_call_Del_cnt = 0;
#else
int snp_call_flag = 0;
#endif

int width,height;
struct image_type image;
struct image_type *im; //  = &image;
#endif

#define NEEDSQL 0
// needsql is for my tcga and target stuff.  private access only.  sadly.  they really should make this data public 

#define FAVEBACKGROUNDCOLOR "F8F2E4"

int ih = 350; // Image Height, 350 picked as default
int iw = 1000; // Image Width, 1000 picked as default
char blds[512];   // "genomic build" string : example "hg19"
char fn_bam[FILENAME_MAX];    // bam file name
int did_content_type = 0; // for web services , "started output" , this means printed "Content-type: text/html\n\n" to output for html rendering
int inited_colors_flag = 0;  // used to init colors for graphics

char CONFIGFILE[] = "alview.conf";
char URL[512];
char HTDOCSDIR[512];
char CGIDIRIFNECESSARY[512];
// old char FULL_PATH_TO_TTF[512]; //   "/usr/X11R6/lib/X11/fonts/TTF/luximb.ttf" ?
char GENOMEDATADIR[512];  // eample /home/rfinney/MD/ALV/release/GENOMEDATA  will append blds (build) string 

void freedom_for_memory(void);
void free_filez(void);

#ifdef MAC
//strlcpy already defined
#else
size_t strlcpy(char *dst, char *src, size_t siz)
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
			*d = '\0';		/* NUL-terminate dst */
		while (*s++)
			;
	}

	return(s - src - 1);	/* count does not include NUL */
}
#endif

#define NCI_RED "A90101"
   // nice red color associated with NCI brand

#ifdef UNIX
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h>
#endif

char cwd[FILENAME_MAX]; 
char ced[FILENAME_MAX]; 

#ifdef WIN32
char hacked_ced[FILENAME_MAX];  //  windows - for current executing directory
char hacked_ced_for_sprintf[FILENAME_MAX];       // windows - for URL 
#endif

void full_path_to_current_and_executing_directory(char cwd[], char ced[])
{
#ifdef WIN32
    // Windows:
    int i,j;
    int bytes = 0;
    TCHAR cExePath[FILENAME_MAX];
    char t[FILENAME_MAX];
    TCHAR cCurrentPath[FILENAME_MAX];

    strcpy(cwd,".");
    strcpy(ced,".");
    cExePath[0] = (DWORD)0;
    cCurrentPath[0] = (DWORD)0;

#if 0
// worked on qt
    bytes = GetModuleFileName(NULL, cExePath, FILENAME_MAX);
    if(bytes == 0) { } // return -1;
    else 
    {
                // dos   size_t wcstombs( char *mbstr, const wchar_t *wcstr, size_t count); 
                // unix  size_t wcstombs(char *S, const wchar_t *PWC, size_t N); #include <stdlib.h>
        bytes = wcstombs(ced,cExePath,bytes);
    }
    bytes = GetCurrentDirectory(FILENAME_MAX,cCurrentPath);
    bytes = wcstombs(cwd,cCurrentPath,bytes);
#endif

    bytes = GetModuleFileName(NULL, t, FILENAME_MAX);
    if(bytes == 0) { } // return -1;
    strcpy(ced,t);
    bytes = GetCurrentDirectory(FILENAME_MAX,t);
    strcpy(cwd,t);


    char drive[FILENAME_MAX];
    char dir[FILENAME_MAX];
    char fname[FILENAME_MAX];
    char ext[FILENAME_MAX];

    _splitpath(ced,drive,dir,fname,ext);
#if QT_GUI
    sprintf(ced,"%s\%s",drive,dir);
#else
    sprintf(ced,"%s\\%s",drive,dir);
#endif

    strcpy(hacked_ced,ced);
    for (i=0;ced[i];i++)
    {
        if (ced[i] == '\\') hacked_ced[i] = '/';
        else hacked_ced[i] = ced[i];
    }
    hacked_ced[i] = (char)0;

    for (j=i=0;ced[i];i++)
    {
        if (ced[i] == '\\') { hacked_ced_for_sprintf[j++] = '\\'; hacked_ced_for_sprintf[j++] = '\\'; }
        else hacked_ced_for_sprintf[j++] = ced[i];
    }
    hacked_ced_for_sprintf[j] = (char)0;


    return;
#endif

#ifdef UNIX
// Linux:
    char szTmp[32];
    char tmps[FILENAME_MAX];
    int bytes;
    int len = FILENAME_MAX;
    sprintf(szTmp, "/proc/%d/exe", getpid());
    bytes = readlink(szTmp, ced, len);
    if (bytes >= 0)
    	ced[bytes] = '\0';
    strcpy(tmps,ced);
    char *z = dirname(tmps);
    if (z) strcpy(ced,z);
    getcwd(cwd,FILENAME_MAX);
    return;
#endif
}
#if 0 // needed for really old compilers
long ftell(FILE *stream);
off_t ftello(FILE *stream);
#endif

#ifdef WIN32
off_t ftello(FILE *stream);
#endif

 /// --- samtools includes ...
#ifdef WIN32
#include "winwacky.h"
#endif


#if NATIVE
// #include "minibam.h" 
#include "alview.h"
#else
#endif
#include "bam.h"
#include "sam.h"
#include "kstring.h"


#if USE_MYSQL
// dont worry about this unless you want to hook up your own database.  
// Regardless you probably don't want to deal with UCSC cruft and their pathological fear of the stanadard library and endless piles of abstration . 

#include "/usr/local/mysql-standard-5.0.27-linux-x86_64-icc-glibc23/include/mysql.h"

// #define NEEDCONFIG 1
#include "cbiitcommon.c"
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "cheapcgi.h"
#include "htmshell.h"
#include "obscure.h"
#include "web.h"
#include "cart.h"
#include "hdb.h"
#include "dbDb.h"
#include "hgFind.h"
#include "hCommon.h"
#include "hui.h"
#include "customTrack.h"
#include "../lib/cbiit.h"

#endif


// #define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)

struct disease_type
{ // tcga stuff 
    char *ln;
    char *sn;
} diseases[] = 
{
    {(char *)"Acute Myeloid Leukemia", (char *)"LAML" } , 
    {(char *)"Bladder Urothelial Carcinoma", (char *)"BLCA" } , 
    {(char *)"Brain Lower Grade Glioma", (char *)"LGG" } , 
    {(char *)"Breast invasive carcinoma", (char *)"BRCA" } , 
    {(char *)"Cervical Squamous Cell Carcinoma", (char *)"CESC" } , 
    {(char *)"Colon adenocarcinoma", (char *)"COAD" } , 
    {(char *)"Glioblastoma multiforme", (char *)"GBM" } , 
    {(char *)"Head and Neck squamous cell carcinoma", (char *)"HNSC" } , 
    {(char *)"Kidney renal clear cell carcinoma", (char *)"KIRC" } , 
    {(char *)"Kidney renal papillary cell carcinoma", (char *)"KIRP" } , 
    {(char *)"Liver hepatocellular carcinoma", (char *)"LIHC" } , 
    {(char *)"Lung adenocarcinoma", (char *)"LUAD" } , 
    {(char *)"Lung squamous cell carcinoma", (char *)"LUSC" } , 
    {(char *)"Ovarian serous cystadenocarcinoma", (char *)"OV" } , 
    {(char *)"Pancreatic adenocarcinoma", (char *)"PAAD" } , 
    {(char *)"Prostate adenocarcinoma", (char *)"PRAD" } , 
    {(char *)"Rectum adenocarcinoma", (char *)"READ" } , 
    {(char *)"Stomach adenocarcinoma", (char *)"STAD" } , 
    {(char *)"Thyroid carcinoma", (char *)"THCA" } , 
    {(char *)0, (char *)0 } 
};


#if WIN32
// jdebug is in alvwin32.cpp file
void jdebug(const char *s)         // for use in debuging 
{
 // V--- if this set to zero then off, 
#if 0
 // ^--- this thing here  - if it is set to 1, then it writes to a file called "x"

    static FILE *debugfp = (FILE *)0;
    static int cnt = 0;
//     int i = 0;
    char fn[2048];
//     char m[2048];

    strcpy(fn,"x"); 
    if (!debugfp) debugfp = fopen(fn,"a");
    if (debugfp == (FILE *)0) { return; }
    fprintf(debugfp,"%d:%s\n",cnt++,s);
    fflush(debugfp); 
//     fclose(debugfp);
#endif
    return;
}
#else

void jdebug(const char *s) // for use in debuging on internet
{
#if 1
// ^^^^ TO TURN OFF IF SET THIS TO ONE 
// fprintf(stderr,"%s\n",s); fflush(stderr);
#else
    register int i;
    time_t timer;
    char fn[256];
    char m[512];
    FILE *fp;

/* printf("DEBUG: in jdebug %s <br>\n",s); printf("%s\n",s); */

#ifdef WEB_SERVER
    strcpy(fn,"/tmp/x");
#else
    strcpy(fn,"x");
#endif

    fp = fopen(fn,"a");
    if (fp == (FILE *)0) { return; }

    timer=time(NULL);
    sprintf(m,"%s",asctime(localtime(&timer)));
    for (i=0;m[i];i++) if (m[i] == '\n') {m[i]=(char)0; break; }
    fprintf(fp,"[%s]%s\n",m,s);
    fclose(fp);

    return;
#endif
}
#endif


int samoption = 0; 
int global_code = 0;
int global_fastafastq_option = 0; // 0=fasta,1=fastq
char global_khr[126];


// cat /h1/finneyr/WIGZ/bamf | awk '{FS="/"}{print $NF" "$0}'  |sort | awk '{print $2}' | awk '{print "    \""$0"\" , "}'

struct filez_type
{
    char *fullpath;
    char *info;
// **** IF fullpath is null, it's just an information line  
    int file_id;
};

struct filez_type *filez;
int num_filez;
int filez_id; // actual id, not the index anymore , can't call it "fileno" which is in stdio.h,  cgi argument "filenum" now refers to this, is >= 1
int filez_index; // index into filez
int num_ids;

int get_filez_index_for_id(int file_id_arg)
{
    int i;
    for (i=0;i<num_filez;i++)
       if (filez[i].file_id == file_id_arg) return i;
    return -1;
}


void load_filez() // only called by webserver
{
    size_t len = 0;
    int i = 0;
    int j = 0;
    int bamcnt = 0;
    FILE *fp = (FILE *)0;
    char s[2000];
    char t[2000];
    char junque[2000];

    s[0] = t[0] = junque[0] = (char)0;

    fp = fopen("list.public","r"); // this file is a list of full path names to bam files
    if (!fp) 
    {
        jdebug("ERROR: can't open list.public in cgi-bin directory to read bam files to present in the front menu."); 
        return;
    }
    num_filez = 0;
    while ( fgets(s,1024,fp) )
       num_filez++;
    rewind(fp);
    filez = (struct filez_type *)malloc(sizeof(struct filez_type) * (num_filez+1) );
    memset(filez,0 ,sizeof(struct filez_type) * (num_filez+1)); // filez 
    i = 0;
    bamcnt = 0;
    while ( fgets(s,1024,fp) )
    {
// sprintf(m,"in load_filez s=%s ",s);jdebug(m);  
        filez[i].fullpath = filez[i].info = (char *)0;
        for (j=0;s[j];j++) { if ((s[j] == '\n') || (s[j] == '\r')) {s[j] = (char)0; break; }}
        if (s[0] == '#') 
        {
// sprintf(m,"in load_filez --- comment ");  jdebug(m); 
            filez[i].fullpath = (char *)0;
            filez[i].info = strdup(s);
            filez[i].file_id = -1;
        }
        else
        {
            sscanf(s,"%s %s",t,junque);    // put path in "t"
// sprintf(m,"in load_filez --- t is %s ",t);  jdebug(m); 
            len = strlen(t);
            if (len) filez[i].fullpath = strdup(t);
            else     filez[i].fullpath = (char *)0;
            if ((s[len+1])&&(s[len+2]))   // does have whitespace and (maybe) something beyond a "white" char
                filez[i].info = strdup(&s[len+1]);
            else 
                filez[i].info = (char *)0;
            filez[i].file_id = bamcnt+1;
            bamcnt++;
        }
        memset(s,0,sizeof(s));
        i++;
    }
    num_ids = bamcnt;
}


void free_filez(void)
{
    int i;

    if (filez == (struct filez_type *)0) return;
    for (i=0 ; i<num_filez ; i++)
    {
        if (filez[i].fullpath) free(filez[i].fullpath); 
        if (filez[i].info) free(filez[i].info); 
    }
    free(filez);
    filez = (struct filez_type *)0;
    return;
}



int chrsizes_hg18[] = 
{
/* chr1 */        247249719 , 
/* chr2 */        242951149 , 
/* chr3 */        199501827 , 
/* chr4 */        191273063 , 
/* chr5 */        180857866 , 
/* chr6 */        170899992 , 
/* chr7 */        158821424 , 
/* chr8 */        146274826 , 
/* chr9 */        140273252 , 
/* chr10 */        135374737 , 
/* chr11 */        134452384 , 
/* chr12 */        132349534 , 
/* chr13 */        114142980 , 
/* chr14 */        106368585 , 
/* chr15 */        100338915 , 
/* chr16 */        88827254 , 
/* chr17 */        78774742 , 
/* chr18 */        76117153 , 
/* chr19 */        63811651 , 
/* chr20 */        62435964 , 
/* chr21 */        46944323 , 
/* chr22 */        49691432 , 
/* chrX */        154913754 , 
/* chrY */        57772954 , 
/* chrM */        16571  
/* error */        -1 , 
};

int chrsizes_hg19[] =  
{
/* chr1 */        249250621 , 
/* chr2 */        243199373 , 
/* chr3 */        198022430 , 
/* chr4 */        191154276 , 
/* chr5 */        180915260 , 
/* chr6 */        171115067 , 
/* chr7 */        159138663 , 
/* chr8 */        146364022 , 
/* chr9 */        141213431 , 
/* chr10 */        135534747 , 
/* chr11 */        135006516 , 
/* chr12 */        133851895 , 
/* chr13 */        115169878 , 
/* chr14 */        107349540 , 
/* chr15 */        102531392 , 
/* chr16 */        90354753 , 
/* chr17 */        81195210 , 
/* chr18 */        78077248 , 
/* chr20 */        63025520 , 
/* chr19 */        59128983 , 
/* chr21 */        48129895 , 
/* chr22 */        51304566 , 
/* chrX */        155270560 , 
/* chrY */        59373566 , 
/* chrM */        16571  
/* error */        -1 , 
};


FILE *fp_refflat_d; // refflat fixed size data
FILE *fp_refflat_o; // data pointed to in fp_refflat_d   _o is for "overflow" 


#define MAX_CHROMSIZE 249250621
 // ^^^^ size of hg19 chr1 

#define MAX_GENENAME 48

         /* some more info on refflat data structure is at http://genome.ucsc.edu/goldenPath/gbdDescriptions.html */
struct flattype           /* "Alignment" type */
{
    char geneName[MAX_GENENAME];
    char name[MAX_GENENAME];
    char chrom[MAX_GENENAME];       /* "Target sequence name" */                             /* 14 */
    unsigned int txStart;
    unsigned int txEnd;
    unsigned int cdsStart;
    unsigned int cdsEnd;
    unsigned int exonCount;
    unsigned int *exonStarts; // tricky , really a ptr to exonStart and and exonEnds in another file
    unsigned int *exonEnds;
    unsigned int strand; // '+' = positive, '-' = negative   as integer !!!
};

char refflat_d_fn[FILENAME_MAX];
char refflat_o_fn[FILENAME_MAX];
void fix_up_support_file_paths(void)
{
    char m[2048];
    char binary_refflat_prefix[512];
    static int fixed_up_support_flag = 0;


sprintf(m,"xxx in fix_up_support_file_paths(), blds = %s",blds); jdebug(m);
    if (fixed_up_support_flag  == 1) return;
    if (blds[0] == 0)
    {
sprintf(m,"ERROR: xxx in fix_up_support_file_paths(), blds is null "); jdebug(m);
          return;
    }

    strcpy(refflat_d_fn,""); 
    strcpy(refflat_o_fn,""); 

    if (GENOMEDATADIR[0] ) // let user specify the directory 
    {
        strcpy(binary_refflat_prefix,GENOMEDATADIR);
#ifdef WIN32
        strcat(binary_refflat_prefix,"\\"); // does windows work this way ?
        strcat(binary_refflat_prefix,blds);
        strcat(binary_refflat_prefix,"\\"); 
#else
// linux and osx
        strcat(binary_refflat_prefix,"/");
        strcat(binary_refflat_prefix,blds);
        strcat(binary_refflat_prefix,"/"); 
#endif
        strcat(binary_refflat_prefix,blds);
        strcat(binary_refflat_prefix,".refflat");
        strcpy(refflat_d_fn,binary_refflat_prefix); 
        strcat(refflat_d_fn,".dat"); 
        strcpy(refflat_o_fn,binary_refflat_prefix); 
        strcat(refflat_o_fn,".dao");   // o not t
sprintf(m,"in fix_up_support_file_paths fns = %s %s",refflat_d_fn,refflat_o_fn);  jdebug(m); 
    }
    else
    {
#ifdef WIN32
    sprintf(binary_refflat_prefix,"%s\\%s.refflat",blds,blds);
    sprintf(refflat_d_fn,"%s\\GENOMEDATA\\%s.dat",hacked_ced_for_sprintf,binary_refflat_prefix);
    sprintf(refflat_o_fn,"%s\\GENOMEDATA\\%s.dao",hacked_ced_for_sprintf,binary_refflat_prefix);

#else
    sprintf(binary_refflat_prefix,"%s/%s.refflat",blds,blds);
    sprintf(refflat_d_fn,"%s/GENOMEDATA/%s.dat",ced,binary_refflat_prefix);
    sprintf(refflat_o_fn,"%s/GENOMEDATA/%s.dao",ced,binary_refflat_prefix);
#endif
    }

#if WEB_SERVER
    sprintf(binary_refflat_prefix,"GENOMEDATA/%s/%s.refflat",blds,blds);
    sprintf(refflat_d_fn,"%s.dat",binary_refflat_prefix);
    sprintf(refflat_o_fn,"%s.dao",binary_refflat_prefix);
#endif

#if CMD_LINE
    sprintf(binary_refflat_prefix,"%s/%s.refflat",blds,blds);
    sprintf(refflat_d_fn,"%s.dat",binary_refflat_prefix);
    sprintf(refflat_o_fn,"%s.dao",binary_refflat_prefix);
#endif

    fixed_up_support_flag = 1;
sprintf(m,"END fix_up_support_file_paths fns = %s %s",refflat_d_fn,refflat_o_fn);  jdebug(m); 
    return;
}

// (48*3)+(4*8)
#define SIZE_REFFLAT_REC 176

long int refflat_fixed_hi;
int numflats = 0;
int maxflats = 0;
struct flattype *flats; // used to be hardwired to array of MAX_FLATS, but now must be alloced.
long int last_refflat_fseek_place;


int setup_refflat(char *input_flat_file_name_d, char *input_flat_file_name_o)
{
    int error;
    int status;
    char m[1024];

sprintf(m,"in setup_refflat() d=%s o=%s",input_flat_file_name_d, input_flat_file_name_o);  jdebug(m); 
    if (fp_refflat_d == (FILE *)0)
    {
        fp_refflat_d = fopen(input_flat_file_name_d,"rb");
        if (fp_refflat_d == (FILE *)0)
        {
             sprintf(m,"ERROR: cant open reflatbinaryname=[%s] in setup_refflat() ",input_flat_file_name_d); jdebug(m);
             return -1;
        }
    }

sprintf(m,"in setup_refflat FILE=%p ",fp_refflat_d);  jdebug(m);
    status = fseek(fp_refflat_d,(size_t) 0 , SEEK_END);
    error = errno;
    if (status != 0) 
    {
       sprintf(m,"ERROR: can't fseek to end in file \"%s\" -- in setup in setup_refflat() ] - errno=%d",input_flat_file_name_d,errno); jdebug(m);
       fclose(fp_refflat_d);
       fp_refflat_d = (FILE *)0;
       return -2;
    }
/// --- SETUP the starting max for a binary search 
    refflat_fixed_hi = ftell(fp_refflat_d) / SIZE_REFFLAT_REC;
sprintf(m,"in setup_refflat refflat_fixed_hi=%ld ",refflat_fixed_hi); jdebug(m);

    if (fp_refflat_o == (FILE *)0) 
    {
        fp_refflat_o = fopen(input_flat_file_name_o,"rb");
        if (!fp_refflat_o) 
        { 
            sprintf(m,"ERROR: can't open file named \"%s\" in order setup in setup_refflat() ] \n",input_flat_file_name_o); jdebug(m);
            fclose(fp_refflat_d);
            fp_refflat_d = (FILE *)0;
            return -3;
        }
    }
sprintf(m,"end  setup_refflat refflat_fixed_hi=%ld , returning 0(==good)",refflat_fixed_hi);  jdebug(m);
    return 0;
}


int do_by_gene_name_from_refflat(char gene[],char chr[],int *start,int *end)
{
    char m[5120];
    struct flattype f;
    int rec_cnt = 0;


// -- just blitz through flat file, not using index (cuz not sorted by name) to find match on "gene" parameter
sprintf(m,"in do_by_gene_name_from_refflat - START , gene=[%s] blds=[%s]",gene,blds);  jdebug(m); 

    fix_up_support_file_paths();
    if (setup_refflat(refflat_d_fn,refflat_o_fn) < 0)  // 0 == good
    {
        strcpy(chr,"chr1");
        *start = 1 ; 
        *end = 4042;
sprintf(m,"ERROR: after setup_refflat() - in do_by_gene_name_from_refflat");  jdebug(m);
        return -1;
    }
sprintf(m,"in do_by_gene_name_from_refflat refflat_d_fn=%s, fp_refflat_d=%p ",refflat_d_fn,fp_refflat_d);  jdebug(m); 
    if (fp_refflat_d == (FILE *)0)
    {
// try and open it 
        fp_refflat_d = fopen(refflat_d_fn,"rb");
        if (!fp_refflat_d) 
        { 
           strcpy(chr,"chr1");
           *start = 1 ; 
           *end = 4042;
           sprintf(m,"ERROR: can not open refflat_d_fn = %s ",refflat_d_fn); jdebug(m); 
           return -1;
        }
    }
 
    fseek(fp_refflat_d, 0, SEEK_SET); // rewind 
// xxx rpf 
sprintf(m,"do_by_gene_name_from_refflat,  opened [%s], fseeked to start ",refflat_d_fn);  jdebug(m); 

    rec_cnt = 0;
    while (1)
    {
// sprintf(m,"do_by_gene_name_from_refflat here 2");  jdebug(m); 
        if (fread (&f.geneName,MAX_GENENAME,1,fp_refflat_d) != 1) break;
        if (fread (&f.name,MAX_GENENAME,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.chrom,MAX_GENENAME,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.txStart,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.txEnd,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.cdsStart,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.cdsEnd,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.exonCount,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.exonStarts,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.exonEnds,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&f.strand,4,1,fp_refflat_d) != 1 ) break;

// sprintf(m,"do_by_gene_name_from_refflat CHECK :  %s %s",f.geneName,gene);  jdebug(m); 
        if (strcmp(f.geneName,gene) == 0) 
        {
// sprintf(m,"matched in do_by_gene_name_from_refflat() %s %s",f.geneName,gene);  jdebug(m); 
            *start = (int)f.txStart;
            *end = (int)f.txEnd;
            strcpy(chr,f.chrom); 
/*
---- NO !!!  keep it open
            fclose(fp_refflat_d);
            fp_refflat_d = (FILE *)0;
*/
            return 0;
        }
        rec_cnt++;
    }
sprintf(m,"do_by_gene_name_from_refflat here 999 , rec_cnt = %d",rec_cnt);  jdebug(m); 

/* 
---- NO !!!  keep it open
    fclose(fp_refflat_d);
    fp_refflat_d = (FILE *)0;
*/
    strcpy(chr,"chr1");
    *start = 1; 
    *end = 4043;

    return -1;
}


static int cmpfunc_refflat_by_chrstartend(const void *arg1, const void *arg2)
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


static struct flattype temp_refflat_space;

static int flat_search_count = 0;

static struct flattype *binary_search_refflat_file(struct flattype *match_me, long int lo, long int hi)
{
    int okay = 1;
    long int seekto;
    long int mid;
    int k;
char m[5120];

    if (fp_refflat_d == (void *)0) 
    { 
        sprintf(m,"ERROR: in binary_search_refflat_file() - null fp_refflat_d"); jdebug(m); 
        return (struct flattype *)0; 
    }

#if 0
sprintf(m,"dbg in binary_search_refflat_file() [%s:%d-%d] last_refflat_fseek_place=%ld lo=%ld hi=%ld, cnt=%d",
         match_me->chrom, match_me->txStart, match_me->txEnd, last_refflat_fseek_place,lo,hi,flat_search_count);  jdebug(m);
#endif

    flat_search_count++;

    if (lo>hi) 
    {
//        fprintf(stderr,"ERROR: can't find cuz lo > hi %ld and %ld in binary_search_refflat_file",lo,hi);
        return (struct flattype *)0;
    }
    if (lo==hi)
    {
        seekto = (SIZE_REFFLAT_REC * lo);
        if (fseek(fp_refflat_d,(size_t) seekto , SEEK_SET) < 0)
        {
            sprintf(m,"SEEK to %ld ERROR in binary_search_refflat_file lo=%ld",seekto,lo); jdebug(m);
            exit(0);
        }
        last_refflat_fseek_place = seekto; 

           // no - struct size in gcc != MSVC, use gcc's idea of it : fread(&temp_refflat_space,1,SIZE_REFFLAT_REC,fp_refflat_d);
        okay = 1;
        if (okay) { if (fread (&temp_refflat_space.geneName,MAX_GENENAME,1,fp_refflat_d) != 1) {okay = 0 ; } }  
        if (okay) { if (fread (&temp_refflat_space.name,MAX_GENENAME,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.chrom,MAX_GENENAME,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.txStart,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.txEnd,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.cdsStart,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.cdsEnd,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.exonCount,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.exonStarts,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.exonEnds,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (okay) { if (fread (&temp_refflat_space.strand,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
        if (cmpfunc_refflat_by_chrstartend(&temp_refflat_space,match_me) == 0) 
        {
// fprintf(stderr,"GOT binary_search_refflat_file lo=%ld hi=%ld \n",lo,hi); fflush(stderr); 
            return &temp_refflat_space;
        }
        else 
        {
            return (struct flattype *)0;
        }
    }
    mid = (lo + hi) / 2 ;

    seekto = (SIZE_REFFLAT_REC * mid);
//printf("in search seekto=%ld [ which is mid=%d]\n",seekto,mid); fflush(stdout);

    if (fseek(fp_refflat_d,(size_t) seekto , SEEK_SET) < 0)
    {
        fprintf(stderr,"SEEK to %ld - ERROR in binary_search_refflat_file, mid=%ld\n",seekto,mid);
        exit(0);
    }

    last_refflat_fseek_place = seekto; 
#if 0 
NO! CANT GURANANTEE FIELD ALIGNMENT BETWEEN VARIOUS C COMPILERS !!!  be explicit
    fread( &temp_refflat_space, 1 , SIZE_REFFLAT_REC , fp_refflat_d );
#endif
    okay = 1;
    if (okay) { if (fread (&temp_refflat_space.geneName,MAX_GENENAME,1,fp_refflat_d) != 1) {okay = 0 ; } }  
    if (okay) { if (fread (&temp_refflat_space.name,MAX_GENENAME,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.chrom,MAX_GENENAME,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.txStart,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.txEnd,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.cdsStart,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.cdsEnd,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.exonCount,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.exonStarts,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.exonEnds,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }
    if (okay) { if (fread (&temp_refflat_space.strand,4,1,fp_refflat_d) != 1 )  {okay = 0 ; } }

// sprintf(m,"dbg read in [%s] %s:%d-%d ", temp_refflat_space.geneName, temp_refflat_space.chrom, temp_refflat_space.txStart, temp_refflat_space.txEnd);

// printf("dbg read in %s:%d-%d at last_refflat_fseek_place =%ld \n",temp_refflat_space.chrom,temp_refflat_space.txStart,temp_refflat_space.txEnd,last_refflat_fseek_place);
    k = cmpfunc_refflat_by_chrstartend(&temp_refflat_space,match_me);
// fprintf(stderr,"after cmpfunc_refflat_by_chrstartend() k=%d \n",k);

    if (k == 0)  // less likely !!!
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



void close_refflat(void)
{
    fclose(fp_refflat_d);
    fp_refflat_d = (FILE *)0;
    fclose(fp_refflat_o);
    fp_refflat_o = (FILE *)0;
}



int is_in_exon(FILE *fp, long int spot, int cnt, int loc)
{
    int i = 0;
    size_t read_status = 0;
    long int start = 0;
    long int end = 0;

    fseek(fp, spot, SEEK_SET);
    for (i=0;i<cnt;i++)
    {
        read_status = fread(&start,sizeof(long int),1,fp); if (read_status != 1) break;
        read_status = fread(&end,sizeof(long int),1,fp); if (read_status != 1) break;
        if ((loc >= start) && (loc < end)) return 1;
    }
    return 0;
}



int chr2khri_alternate(char *chrom)
{
    int k;

    k = -1;
         if (strcmp(chrom,"X")==0) {  k = 22; }
    else if (strcmp(chrom,"Y")==0) {  k = 23; }
    else if (strcmp(chrom,"MT")==0) {  k = 24; }
    else if (strcmp(chrom,"M")==0) {  k = 24; }
    else
    {
        k = atoi(&chrom[0]);
        k = k - 1;
    }
    return k;
}



struct chridxtype
{
    int chridx;
    int chrlast;
    int size;
    char chrom[40];
};

struct chridxtype chridexes[25];

void dbg_chridx(void)
{
   int i;
   for (i=0;i<25;i++)
   {
char m[512];
       sprintf(m,"chridexes i=%d idx=%d last=%d sz=%d %s",i,chridexes[i].chridx,chridexes[i].chrlast,chridexes[i].size,chridexes[i].chrom); 
       jdebug(m); 
   }
}



long int totalcnt = 0;
int dcnt = 0;
int icnt = 0;
int nocigarcount = 0;
int unmapcnt = 0;
int opticaldupecount = 0;
int badqualcount = 0;
char global_lpgcookie[2048];

/* --- globals ... */
#define PREFIX ""
#define OPERATION_COMPARE 0
#define OPERATION_GENES 1



char REQUEST_METHOD[512];
char QUERY_STRING[512];
char REMOTE_ADDR[512];
char HTTP_USER_AGENT[512];

#if USE_MYSQL
#else
static int MOZILLA = 0;
static int MSIE = 1;
static int WINDOWS = 1;
#endif

#if 0
static char FECONFFILE[256];     /* = "/usr/local/apache/fe/fe.conf"; */
static char FEUSERDATADIR[256];  /* = "/usr/local/apache/fe/userdata/"; */
static char FESYSDATADIR[256];   /* = "/usr/local/apache/fe/sysdata/"; */
static char FEHTMLDIR[256];      /*  "/usr/local/apache/htdocs/fe/data/"; */
static char FEURL[256];          /* http://127.0.0.1/ */
static char FEBBURL[256];        /* http://127.0.0.1/ */
static char FEHTMLRELDIR[256];   /* fe/ */
static char FECGI[256];          /* cgi-bin/ */
static char FEBASEDIR[256];
#endif



int gs; // global start genomic cooordinate 
int ge; // global end genomic cooordinate 
int gdiff;   // start - end  = gs - ge
int menu; // main menu
int pickmenu;
char pickmenustring[MAXBUFF];
char menustring[MAXBUFF];
char forcefilecgiarg[MAXBUFF];
char pickq[MAXBUFF];
char bamname[MAXBUFF];
char filecgiarg[MAXBUFF];
char jutlfn[MAXBUFF]; // "just use this long file name"
char position[MAXBUFF];
char op[MAXBUFF];
char chr[MAXBUFF];




    /* bowtie stuff */
static int bi1 , bi2;

#define  QUALCUTOFF 6
int spliceonly_flag =0;  // "splice on"
int qual_flag = 0; // quality (not sequence) 
int mmspot_flag;  // flag for whether to draw the red up/down bar at mm spot
int extrapolate_flag;  // extrapolate flag
int nso_flag;  // "non specifics on"
// int altonly_flag =0;  // "non ref splice only on"
int geneannot_flag =1;
int basecolors_flag =0;
int uniq_flag =0;
int tds_val = 3;  // data set id 

#define JPG_IMAGE 0
#define PNG_IMAGE 1
#define GIF_IMAGE 2


void alview_load_config(void)
{
char m[1024];
    int i,j,k,mode;
    FILE *fp;
    char s[1026];
    char t1[500];
    char t2[500];

// jdebug("in alview_load_config"); 

    t1[0]= t2[0] = (char)0;
    fp = fopen(CONFIGFILE,"r");
    if (fp == (FILE *)0) { return; }
    while ( fgets(s,1024,fp) )
    {
        mode = i = j = k = 0;
        while (s[i])
        {
            if (s[i] < ' ') break;
            if (s[i] == '=') { s[i] = ' '; mode = 1;}
            else
            {
                if (mode == 0)
                {
                    t1[j++] = s[i];
                    t1[j] = (char)0;
                }
                else if (mode == 1)
                {
                    t2[k++] = s[i];
                    t2[k] = (char)0;
                }
            }
            i++;
        }
/* add global configs here ... */
        if (strcmp(t1,"URL") == 0)                    strcpy(URL,t2);        // URL=lpgws501.nci.nih.gov
        else if (strcmp(t1,"HTDOCSDIR") == 0)         strcpy(HTDOCSDIR,t2);
        else if (strcmp(t1,"CGIDIRIFNECESSARY") == 0) strcpy(CGIDIRIFNECESSARY,t2); 
// old        else if (strcmp(t1,"FULL_PATH_TO_TTF") == 0) strcpy(FULL_PATH_TO_TTF,t2); 
        else if (strcmp(t1,"GENOMEDATADIR") == 0) 
        {
             strcpy(GENOMEDATADIR,t2); 
sprintf(m,"in alview_load_config  GENOMEDATADIR %s ",GENOMEDATADIR);  jdebug(m);
        }
    }
    fclose(fp);
}


void alview_init(void)
{
    alview_load_config();
    full_path_to_current_and_executing_directory(cwd, ced);
}


//---- 2bit format from UCSC 
struct n_type // 'N' as in not in (ACGT, the 4 nucleotides)
{
    unsigned int nStarts;
    unsigned int nSizes;
};


static char *do_2bit_dna(FILE *fp, unsigned int fragStart, unsigned int fragEnd, struct n_type *Ns, unsigned int NBlockCount)
{
    char *dna_original = (char *)0;
    char *dna = (char *)0;
    unsigned char uc;
    unsigned char uc2;
    unsigned int val2bit;
    unsigned int spot;
    int modder = 0;
    unsigned int size = 0;
    unsigned int i,j;
    unsigned int k;
    unsigned int fsm1 = fragStart - 1; // Frag Start Minus 1


    i = j = 0;
    k = 0;
    if (fp == (FILE *)0) return (char *)0;

/* Skip to bits we need and read them in. */
    spot = (fsm1) / 4;
    modder = (fsm1) % 4;
    size = fragEnd - fragStart + 1;
    dna_original = (char *)malloc(size+8); // plus slack for now
    if (!dna_original) return((char *)0); 
    dna = dna_original;

    fseek(fp,spot,SEEK_CUR);

#define BIT2VAL(VAL2BIT) \
        if (VAL2BIT == 0) uc2 = 'T'; \
        else if (VAL2BIT == 1) uc2 = 'C'; \
        else if (VAL2BIT == 2) uc2 = 'A'; \
        else uc2 = 'G'; \
        *dna++ = uc2;

    i = 0;
// handler first 4 2bits  (byte) 
    if (modder != 0) 
    {
        if (fread(&uc,1,1,fp) != 1) { free(dna_original); return (char *)0; }
// may not print out all 4 twobits
        if ((modder <= 0) && (i < size))
        { val2bit = (uc >> 6) & 3; BIT2VAL(val2bit); i++; }
        if ((modder <= 1) && (i < size))
        { val2bit = (uc >> 4) & 3; BIT2VAL(val2bit); i++; }
        if ((modder <= 2) && (i < size))
        { val2bit = (uc >> 2) & 3; BIT2VAL(val2bit); i++; }
        if ((modder <= 3) && (i < size))
        { val2bit = (uc >> 0) & 3; BIT2VAL(val2bit); i++; }
    }

// handle remainder of bytes 
    for ( ; i<size ; )
    {
        if (fread(&uc,1,1,fp) != 1) { free(dna_original); return (char *)0; }
        val2bit = (uc >> 6) & 3; BIT2VAL(val2bit) if (++i >= size) break;
        val2bit = (uc >> 4) & 3; BIT2VAL(val2bit) if (++i >= size) break;
        val2bit = (uc >> 2) & 3; BIT2VAL(val2bit) if (++i >= size) break;
        val2bit =        uc & 3; BIT2VAL(val2bit) if (++i >= size) break;
    }

    for (i=0 ; i<NBlockCount ; i++) 
    {
/*
    |    |      <-- nstart,nend
... fragstart,fragends cases ...
| |                      no    - pre
|    |                   yes - overlap before 
      ||                 yes  - inclusive
       |    |            yes  - overlap after
  |             |        yes  - before and after
          |    |         no    -after 
*/
        if (((Ns+i)->nStarts+1) > fragEnd)
            break;
        if (
            ((fragStart >= ((Ns+i)->nStarts+1)) && (fragStart < ((Ns+i)->nStarts+1 + (Ns+i)->nSizes))) || 
            ((fragEnd   >= ((Ns+i)->nStarts+1)) && (fragEnd   < ((Ns+i)->nStarts+1 + (Ns+i)->nSizes)))
           )
        {
// mask off N's
            for (j=(Ns+i)->nStarts+1; j < ( (Ns+i)->nStarts + (Ns+i)->nSizes + 1); j++)
            {
                k = j - fragStart;
                if (k < 0) continue;
                if (k >= (fragEnd - fragStart+1)) 
                {
                    break;
                }
                if (k >= 0) 
                {
                    *(dna_original+k) = 'N';
                }
            }
        }
        else continue;
    }
    return dna_original;
}


char *rd2bit(char *fn2bit,char *chr,unsigned int fragStart, unsigned int fragEnd, int *ret)
{
    int error;
    char *dna;
    char name[256];
    FILE *fp = (FILE *)0;
    struct n_type *Ns;
    unsigned int NBlockCount;
    unsigned int size;
    unsigned int ui;
    unsigned int signature;
    unsigned int version;
    unsigned int sequenceCount;
    unsigned int maskBlockCount;
    unsigned char namesize;
    unsigned int reserved;
    unsigned int offset;
    unsigned int i;
char m[1024];


    *ret = 0;
    fp = fopen(fn2bit,"rb"); 
    error = errno;
    if (fp == (FILE *)0) 
    { 
        sprintf(m,"ERROR: in rd2bit, fopen failed, errno = %d",error); jdebug(m);  
        *ret= -1; 
        return (char *)0; 
    }  

    if ((fread(&signature,4,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -2; return((char *)0); }
    if ((fread(&version,4,1,fp)) != 1)       { fclose(fp); fp = (FILE *)0; *ret= -3; return((char *)0); }
    if ((fread(&sequenceCount,4,1,fp)) != 1) { fclose(fp); fp = (FILE *)0; *ret= -4; return((char *)0); } // number of sequences in the file
    if ((fread(&reserved,4,1,fp)) != 1)      { fclose(fp); fp = (FILE *)0; *ret= -5; return((char *)0); }
    for (i=0 ; i<sequenceCount ; i++)
    {
        if ((fread(&namesize,1,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -6; return((char *)0); }
        if ((fread(name,(size_t)namesize,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -7; return((char *)0); }
        name[namesize] = (char)0;
        if ((fread(&offset,4,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -8; return((char *)0); } // SEEK HERE if we match the name 
        if (strcmp(chr,name) == 0)  // found it
            break;
    }
    if (i == sequenceCount) { fclose(fp); fp = (FILE *)0; *ret= -9; return((char *)0); } // cant find name

    fseek(fp,offset,SEEK_SET);

    if (fread(&size,4,1,fp) != 1) { fclose(fp); fp = (FILE *)0; *ret= -10; return((char *)0); }
    if (fread(&NBlockCount,4,1,fp) != 1) { fclose(fp); fp = (FILE *)0; *ret= -11; return((char *)0); }

    Ns = (struct n_type *)malloc(NBlockCount * sizeof (struct n_type));
    if (!Ns) { fclose(fp); fp = (FILE *)0; *ret= -12; return((char *)0);}

    for (i=0;i<NBlockCount;i++)
    {
        if (fread(&ui,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -13; return((char *)0); }
        (Ns+i)->nStarts = ui;
    }
    for (i=0;i<NBlockCount;i++)
    {
        if (fread(&ui,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -14; return((char *)0); }
        (Ns+i)->nSizes = ui;
    }
    if (fread(&maskBlockCount,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -15; return((char *)0); }
    for (i=0;i<maskBlockCount;i++)
    {              //mask starts
        if (fread(&ui,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -16; return((char *)0); }
    }

    for (i=0;i<maskBlockCount;i++)
    {              //mask sizes
        if (fread(&ui,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -17; return((char *)0); }
    }
    if (fread(&ui,4,1,fp) != 1) { free(Ns); fclose(fp); fp = (FILE *)0; *ret= -18; return((char *)0); } // reserved

// okay -so freegan A, we are at the DNA ...
    dna = do_2bit_dna(fp,fragStart,fragEnd,Ns,NBlockCount);

    if (Ns) free(Ns); 
    if (fp) fclose(fp); 
    fp = (FILE *)0;
    if (!dna) *ret= -19;

    return dna;
}

void setup_2bit_file(char fn[], char blds_arg[])
{
   int win=0;

#ifdef WIN32
   win=1;
#endif

    if (GENOMEDATADIR[0]) // let user specify the directory 
    {
        strcpy(fn,GENOMEDATADIR); // user configureable 
        if (win) strcat(fn,"\\"); else strcat(fn,"/");
        strcat(fn,blds_arg);
        if (win) strcat(fn,"\\"); else strcat(fn,"/");
        strcat(fn,blds_arg);
        strcat(fn,".2bit");
    }
    else 
    {
        strcpy(fn,"GENOMEDATA"); // look in default place
        if (win) strcat(fn,"\\"); else strcat(fn,"/");
        strcat(fn,blds_arg);
        if (win) strcat(fn,"\\"); else strcat(fn,"/");
        strcat(fn,blds_arg);
        strcat(fn,".2bit");
    }
}


unsigned int get_chr_size_from_2bit(char *chr, char blds[],int *ret)
{
    char name[256]; // example value : "chr1"
    char fn2bit[FILENAME_MAX];
    FILE *fp = (FILE *)0;
//     char *dna;
//     struct n_type *Ns;
//     unsigned int ui;
//     unsigned int maskBlockCount;
//     unsigned int NBlockCount;
    unsigned int size;
    unsigned int signature;
    unsigned int version;
    unsigned int sequenceCount;
    unsigned char namesize;
    unsigned int reserved;
    unsigned int offset;
    unsigned int i;


    setup_2bit_file(fn2bit,blds); 

    *ret = 0;
    fp = fopen(fn2bit,"rb"); 
    if (fp == (FILE *)0) { *ret= -1; return 0; }  

    if ((fread(&signature,4,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -2; return 0; }
    if ((fread(&version,4,1,fp)) != 1)       { fclose(fp); fp = (FILE *)0; *ret= -3; return 0; }
    if ((fread(&sequenceCount,4,1,fp)) != 1) { fclose(fp); fp = (FILE *)0; *ret= -4; return 0; } // number of sequences in the file
    if ((fread(&reserved,4,1,fp)) != 1)      { fclose(fp); fp = (FILE *)0; *ret= -5; return 0; }
    for (i=0 ; i<sequenceCount ; i++)
    {
        if ((fread(&namesize,1,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -6; return 0; }
        if ((fread(name,(size_t)namesize,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -7; return 0; }
        name[namesize] = (char)0;
        if ((fread(&offset,4,1,fp)) != 1)     { fclose(fp); fp = (FILE *)0; *ret= -8; return 0; } // SEEK HERE if we match the name 
        if (strcmp(chr,name) == 0)  // found it
            break;
    }
    if (i == sequenceCount) { fclose(fp); fp = (FILE *)0; *ret= -9; return 0; } // cant find name

    fseek(fp,offset,SEEK_SET);

    if (fread(&size,4,1,fp) != 1) { fclose(fp); fp = (FILE *)0; *ret= -10; return 0; }
    *ret = 0;
    return size;
}



#if 0 
static int getnib(char *nibfn, int start, int size, char puthere[])
{
    long pos;
    unsigned char ch;
    int i;
    int bVal;
    int byteSize;
    FILE *fp;
    char m[512];
// DNA *valToNtTbl = ((options &  NIB_MASK_MIXED) ? valToNtMasked : valToNt);


    fp = fopen(nibfn,"rb"); 
    if (!fp) 
    {
sprintf(m,"in getnib(), failed to open file %s ",nibfn); jdebug(m); 
        return -1;
    }
    pos = (start>>1) + 8;
// sprintf(m,"in getnib(), fseek to %ld  (from start=%d) ",pos,start); jdebug(m); 
    if (fseek(fp, pos, SEEK_SET) == -1)
    {
        sprintf(m,"ERROR: in getnib(), fseek failed "); jdebug(m); 
        fclose(fp); 
        return -2;
    }
    i = 0;
    if (start & 1) 
    {
       bVal = fgetc(fp);
       puthere[i++] = (bVal&0xf);
    }
    byteSize = (size>>1);
    while (--byteSize >= 0)
    {
        bVal = fgetc(fp);
        if (bVal < 0) { fclose(fp); return -3; }
        puthere[i++] = (bVal>>4);
        puthere[i++] = (bVal&0xf);
    }
    if (size&1) 
    {
        bVal = fgetc(fp);
        puthere[i++] = (bVal>>4);
        if (bVal < 0) { fclose(fp); return -4; }
    }
    i--;
    while (i>=0)
    {
         ch = puthere[i] & 0x7; // silly mask
         if (ch ==  0) ch = 'T';                      //  00
         else if (ch ==  1) ch = 'C';                 //  01
         else if (ch ==  2) ch = 'A';                 //  10
         else if (ch ==  3) ch = 'G';                 //  11
         else if (ch ==  4) ch = 'N';                 // 100
         puthere[i] = ch;
         i--;
    }
    fclose(fp);
sprintf(m,"end getnib()"); jdebug(m); 
    return 0;
}
#endif


char *dnaspace; // reference genome data 
unsigned short int *dnacnts;
unsigned short int *dnamms;
int dnasize;

void setup_dnacnts_and_dnamms(int size)
{
    char m[1024];
    unsigned short int *z;

jdebug("in setup_dnacnts_and_dnamms() "); 
    if (dnacnts) free(dnacnts);
    dnacnts = (unsigned short int *)0;
    dnasize = size;
    z = (unsigned short int *)malloc((size*sizeof(unsigned short int))+32); // dnacnts - coverage count
    if (!z) 
    {
// NB: This will actually happen on win32
sprintf(m,"ERROR: can not malloc space in setup_dnacnts_and_dnamms(%d) 1",size); jdebug(m);
         return;
    }
    memset(z,0,((size*sizeof(unsigned short int)))+32);
    dnacnts = z;
    
    if (dnamms) free(dnamms);
    dnamms = (unsigned short int *)0;
    z = (unsigned short int *)malloc((size*sizeof(unsigned short int))+32); // dnamms - mis matches
    if (!z) 
    {
// NB: This will actually happen on win32
sprintf(m,"ERROR: can not malloc space in setup_dnacnts_and_dnamms(%d) 2",size); jdebug(m);
         free(dnacnts);
         dnacnts = (unsigned short int *)0;
         return;
    }
    memset(z,0,((size*sizeof(unsigned short int)))+32);
    dnamms = z;

jdebug("end setup_dnacnts_and_dnamms() "); 
    return;
}


#if 0 // old???? maybe delete 
static char *nib(char khr[], int s, int e,char blds_arg[] )
{
    int status;
    char *z;
    int size;
    char fn[512];
    char m[512];

/* examples 
/app1/gbdb/hg18/nib/chr8.nib
/app1/gbdb/hg18/nib/chr14.nib
*/
    size = e - s;
sprintf(m,"start in nib khr=%s s=%d e=%d size=%d",khr,s,e,size); jdebug(m);
#ifdef QT_GUI
#ifdef WIN32
    if (GENOMEDATADIR[0] ) // let user specify the directory 
    {
        strcpy(fn,GENOMEDATADIR);
        strcat(fn,"\\"); // does windows work this way ?
        strcat(fn,blds);
        strcat(fn,"\\");  // windows uses backslash \ not forward /
        strcat(fn,khr);
        strcat(fn,".nib");
    }
    else sprintf(fn,"%s\\GENOMEDATA\\%s\\%s.nib",ced,blds_arg,khr);
#else
    if (GENOMEDATADIR[0] ) 
    {
        strcpy(fn,GENOMEDATADIR);
        strcat(fn,"/");
        strcat(fn,blds);
        strcat(fn,"/");
        strcat(fn,khr);
        strcat(fn,".nib");
    }
    else sprintf(fn,"%s/GENOMEDATA/%s/%s.nib",ced,blds_arg,khr);
#endif
#else
    sprintf(fn,"/app1/gbdb/%s/nib/%s.nib",blds_arg,khr);
#endif

sprintf(m,"in nib, fn = %s",fn); jdebug(m);
    z = (char *)malloc(size+10); // nibgenome with slack  BE SURE TO FREE THIS, silly
    if (!z) 
    {
sprintf(m,"ERROR: can notmalloc space in nib(%d)",size); jdebug(m);
         return (char *)0;
    }

    memset(z,0,size+10);
    status = getnib(fn, s-1,e-s-1,z);
sprintf(m,"status from getnib() is %d ",status); jdebug(m);
    return z;
}
#endif

int hex2int(int c)
{
    if (c >= 'a' && c <= 'f') return(c - 'a' + 10);
    if (c >= 'A' && c <= 'F') return(c - 'A' + 10);
    if (c >= '0' && c <= '9') return(c - '0');
    return(-1);
}

void un_escape(unsigned char s[])
{
     unsigned char *l = s;
     unsigned char *o = s;
     unsigned char c;
 
     while (*l) 
     {
         switch(*l) 
         {
             case '+': *o = ' ';
                 break;
             case '%': if (l[1] != '%') 
             {
                 if (isxdigit(l[1]) && isxdigit(l[2])) 
                 {
                     c = (hex2int(l[1]) << 4) + hex2int(l[2]);
                     if (c < ' ' || c == 127 ) c = '_';
                     *o = c;
                     l += 2;
                     break;
                 }
                     else l++;
             } 
             default: *o = *l;
         }
         o++;
         l++;
     }
     *o = '\0';
     return;
}


#if USE_MYSQL
#else
void serious_problem(char *msg)
{
    printf("Content-type: text/html\n\n");
    printf("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN>\n");
    printf("<html>\n");
    printf("ERROR: %s (serious problem) <br>\n",msg);
    printf("<br>\n");
    printf("</html>");
    return;
}


void bark()
{
    printf("Content-type: text/html\n\n");
    printf("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN>\n");
    printf("<html>\n");
    printf("quit trying to hack this site.<br>\n");
    printf("no blind strcpy's here<br>\n");
    printf("plese hack here instead : <a href=\"http://www.iwillnothackthissite.com\">http://www.iwillnothackthissite.com</a><br>\n");
    printf("<br>\n");
    printf("happy hacking!!!<br>\n");
    printf("</html>");
    return;
}


void bomb(char *s)
{
    printf("Content-type: text/html\n\n");
    printf("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN>\n");
    printf("<html>\n");
    printf("bomb<br>\n");
printf("%s",s);
    printf("</html>");
fflush(stdout);
    return;
}


void debug1(char *s)
{
    printf("Content-type: text/html\n\n");
    printf("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN>\n");
    printf("<html>\n");
    printf("<h1>NOTE: %s</H1><br>\n",s);
    printf("hack here instead : <a href=\"http://www.iamanevilhacker.com\">http://www.iamanevilhacker.com</a><br>\n");
    printf("<br>\n");
    printf("happy hacking!!!<br>\n");
    printf("</html>");
    return;
}


char evars[20][80]={"SERVER_SOFTWARE", "SERVER_NAME", "SERVER_PROTOCOL", 
                    "SERVER_PORT",
                    "GATEWAY_INTERFACE", "REQUEST_METHOD", 
                    "PATH_INFO", "PATH_TRANSLATED", "SCRIPT_NAME", 
                    "QUERY_STRING", 
                    "REMOTE_HOST", "REMOTE_ADDR", "REMOTE_USER", 
                    "REMOTE_IDENT",
                    "AUTH_TYPE", "CONTENT_TYPE", "CONTENT_LENGTH", 
                    "HTTP_ACCEPT", "HTTP_USER_AGENT", "HTTP_REFERER"};

#define MAX_CONTENTS 500
#define MAX_LEN 4000
struct varval_type
{
    char var[MAX_LEN];   
    char val[MAX_LEN];   
};

struct varval_type cgi_contents[MAX_CONTENTS];
int num_varvals = 0;

void parse_cgi_contents(char *s)
{
    int i,j;
    int mode;
    unsigned char ch;
    char var[MAX_LEN+1];
    char val[MAX_LEN+1];
    char m[MAXBUFF];

    memset(var,0,sizeof(var));
    memset(val,0,sizeof(val));
    memset(cgi_contents,0,sizeof(cgi_contents));


sprintf(m,"in parse_cgi_contents[%s]",s); jdebug(m);

    num_varvals = 0;
    mode = 0; /* on var */
    j = i = 0;
    while (s[i])
    {
        if (s[i] == '+') s[i] = ' ';
        if (s[i] == '=')
        {
            var[j] = (char)0;
            mode = 1;
            j = 0;
            strncpy(cgi_contents[num_varvals].var,var, MAX_LEN-1);
            memset(var,0,sizeof(var));
        }
        else if (s[i] == '&')
        {
            val[j] = (char)0;
            mode = 0;
            j = 0;
            strncpy(cgi_contents[num_varvals].val,val, MAX_LEN-1);
            memset(val,0,sizeof(val));
            num_varvals++;
        }
        else if (s[i] == '%')
        {
            ch = ((s[i+1]-'0') * 16 ) + (s[i+2] - '0') ;
sprintf(m,"ch = %c ",ch); jdebug(m);
            if (mode == 0)
                var[j] = ch;
            else
                val[j] = ch;
            j++;
            i += 2;
        }
        else
        {
            ch = s[i];
            if (mode == 0)
                var[j] = ch;
            else
                val[j] = ch;
            j++;
            if (mode == 0)
                var[j] = (char)0;
            else
                val[j] = (char)0;
        }
        i++;
    }

    if ( strlen(cgi_contents[num_varvals].var) )
    {
        strncpy(cgi_contents[num_varvals].val,val, MAX_LEN-1);
        memset(val,0,sizeof(val));
        num_varvals++;
    }

#if 0
    for (i=0 ; i<num_varvals ; i++)
        printf("prsd: VAR=[%s] VAL=[%s]<br>\n", cgi_contents[i].var, cgi_contents[i].val);
#endif

    return;
}


void print_hide_cgi_contents(FILE *fp)
{
    register int i;

    for (i=0 ; i<num_varvals ; i++)
    {
        fprintf(fp,"<!-- VAR=[%s] VAL=[%s] -->\n", cgi_contents[i].var, cgi_contents[i].val);
    }
    return;
}


int get_cgi_param(char *match, char *retval,int skipamount)
{
    register int i;
    int skipped = 0;


#if 0 // some debug ...
char m[5120];
for (i=0 ; i<num_varvals ; i++)
{
    sprintf(m,"VAR=[%s] VAL=[%s]", cgi_contents[i].var, cgi_contents[i].val);
    jdebug(m);
}
#endif


    for (i=0 ; i<num_varvals ; i++)
    {
        if (strcmp(cgi_contents[i].var,match) == 0)
        {
            skipped++;
            if (skipped > skipamount)
            {
                strcpy(retval,cgi_contents[i].val); 
//sprintf(m,"in get_cgi_param(), num_varvals = %d , match=[%s] GOT IT [%s]",num_varvals,match,retval); jdebug(m);
                return 0; 
            }
        }
    }

    return -1;
}


void log_error(char *msg)
{
    printf("LOG : %s\n",msg);
}



void debug_cgi_environment_variables_or_whatever_they_are()
{
    int j;
    char *ptr2current_env_variable;
    long mylen;
    char *endptr;
    long contentlength;

    for (j=0 ; j<20 ; j++)
    {
        ptr2current_env_variable  = getenv(evars[j]);
        if (ptr2current_env_variable == (char *)0)
        {
            printf("%d NOTE: NO [%s] <br>\n",j+1,evars[j]);
            fflush(stdout);
            continue;
        }

        mylen = strtol(ptr2current_env_variable, &endptr, 10);
        if (strcmp("CONTENT_LENGTH",evars[j]) == 0) /* same ? */      
        {
            contentlength = mylen;
//char m[500];
//sprintf(m,"contentlength = %ld ",contentlength ); jdebug(m);
        }
        printf("%d %s = %s\n<br>", j+1,evars[j], getenv(evars[j]));
    }
}




void get_some_environment_variables()
{
    register char *z;
    char m[256];

    /* let's get some interesting cgi variables */

    REQUEST_METHOD[0] = QUERY_STRING[0] = REMOTE_ADDR[0] = HTTP_USER_AGENT[0] = (char)0;


    z = getenv("REQUEST_METHOD");
    if (z) strncpy(REQUEST_METHOD,z,sizeof(REQUEST_METHOD)-2);

    z = getenv("REMOTE_ADDR");
    if (z) strncpy(REMOTE_ADDR,z,sizeof(REMOTE_ADDR)-2);

    z = getenv("HTTP_USER_AGENT");
    if (z) strncpy(HTTP_USER_AGENT,z,sizeof(HTTP_USER_AGENT)-2);

sprintf(m,"HTTP_USER_AGENT:%s",HTTP_USER_AGENT); 
jdebug(m);

    z = getenv("HTTP_USER_AGENT");
    if (z) strncpy(HTTP_USER_AGENT,z,sizeof(HTTP_USER_AGENT)-2);

// need safari/chrome/ stuff
    if (strstr(HTTP_USER_AGENT,"MSIE"))                     /* figure out user's browser */
    {
/* eg. HTTP_USER_AGENT=Mozilla/4.0 (compatible; MSIE 6.0; Windows 98) */
        MSIE = 1;
        MOZILLA = 0;
    } else if ( (strstr(z,"Netscape")) || (strstr(z,"Gecko"   )) )
    {
       MSIE = 0;
       MOZILLA = 1;
    }

    if (strstr(z,"Windows"))    /* figure out user's OS */
    {
// linux, mac, what-ev-er 
    }
    else
    {
        WINDOWS = 0;
    }
#if 0 /* example useragent strings ... */
HTTP_USER_AGENT:Mozilla/5.0 (Windows; U; Win98; en-US; rv:0.9.4) Gecko/20011128 Netscape6/6.2.1
HTTP_USER_AGENT:Mozilla/5.0 (Windows; U; Win98; en-US; rv:0.9.4) Gecko/20011128 Netscape6/6.2.1
HTTP_USER_AGENT:Mozilla/4.0 (compatible; MSIE 6.0; Windows 98)
HTTP_USER_AGENT:Mozilla/4.0 (compatible; MSIE 5.0; AOL 6.0; Windows 98; DigExt)
HTTP_USER_AGENT:Mozilla/4.0 (compatible; MSIE 5.0; Windows 98; DigExt)
HTTP_USER_AGENT:Mozilla/4.0 (compatible; MSIE 5.5; Windows 98)
HTTP_USER_AGENT:Mozilla/4.0 (compatible; MSIE 6.0; Windows 98)
HTTP_USER_AGENT:Mozilla/5.0 (Windows; U; Win98; en-US; m18) Gecko/20010131 Netscape6/6.01
#endif

}

char *global_QUERY_STRING_ptr = (char *)0;
char buff[MAXBUFF];


#if 0 
int go_front_end(void)
{
char m[MAXBUFF];
    int i = 0;
    long contentlength;
    int devastating_error = 0;
    char *endptr = (char *)0;
    char *ptr2contentlength = (char *)0;
    char *ptr2querystring = (char *)0;
    char dbmsg[MAXBUFF];


           /* make some assumtions until proven otherwise ... */
    WINDOWS = 1;   /* where do you want to go today ??? */
    MSIE = 1;   /* microsoft internet explorer */
    MOZILLA = 0;  /* netscape or mozilla : bow to the beast !!! */

    ptr2querystring = (char *)0;
    ptr2contentlength = getenv("CONTENT_LENGTH");
    if (ptr2contentlength == (char *)0)
    {
        ptr2querystring  = getenv("QUERY_STRING");
        global_QUERY_STRING_ptr = ptr2querystring;
        if (ptr2querystring == (char *)0)
        {
             jdebug("no CONTENT_LENGTH in go_front_end");

             printf("Content-type: text/html\n\n");
             did_content_type = 1;
             printf("<html><body><H1>Invalid Input: No \"Content Length\" and No Query String for program<br>\n");
             debug_cgi_environment_variables_or_whatever_they_are();
             fflush(stdout);
             return -1;
        }
    }

    if (ptr2contentlength)
    {
        contentlength = strtol(ptr2contentlength, &endptr, 10);
        if ( contentlength > MAXBUFF - 10)
        {
            bark();
        }
        fread(buff, contentlength, 1, stdin);
sprintf(m,"zzz before un_escape %s ",buff);  jdebug(m); 
        un_escape((unsigned char *)buff);
sprintf(m,"zzz after  un_escape %s ",buff);  jdebug(m); 

#if 1 /* debug stuff */
printf("Content-type: text/html\n\n");
did_content_type = 1;
printf("<html><body><br>\n");
printf("hello. this is debug. there is a problem<br>\n");
debug_cgi_environment_variables_or_whatever_they_are();
printf("<br>");
printf("buff len=%ld content = [%s]<br>\n",contentlength, buff);
printf("</body></html>\n");
return -1;
#endif
    }
else {
  // sprintf(dbmsg,"ptr2contentlength = NULL "); jdebug(dbmsg); 
     }

    if (ptr2querystring)
    {
        memset(buff,0,sizeof(buff));
        i = strlen(ptr2querystring);
        if (i > (int)(sizeof(buff) - 2) )
            i = (int)(sizeof(buff) - 2); 

        strncpy(buff,ptr2querystring,i);

    }
    else { sprintf(dbmsg,"ptr2querystring = NULL\n"); jdebug(dbmsg); }

    get_some_environment_variables(); /* will also set some global variables */

    printf("Content-type: text/html\n\n");
    did_content_type = 1;

    un_escape((unsigned char *)buff); 
    parse_cgi_contents(buff);

    devastating_error = 0;

    if (devastating_error)
    {
        serious_problem((char *)"ERROR XXX");
        debug_cgi_environment_variables_or_whatever_they_are();
        printf("buff = [%s]",buff);
        printf("<br>");
        return -1;
    }

    jdebug(dbmsg);

    return 0;
}
#endif




void lowerize( char *str )
{
   while( *str != 0 )
   {
      if ( (*str >= 'A') && (*str <='Z') )
      {
         *str += 'a' - 'A';
      }
      str++;
   }
}


#endif // USE_MYSQL


int cutlohiflag = 0;
double cutlo;
double cuthi;

char hei[MAXBUFF]; // height
char wid[MAXBUFF]; // width
char fl1[MAXBUFF]; // file 1
char fl2[MAXBUFF]; // file 2
char los[MAXBUFF];
char his[MAXBUFF];
// OLD char typ[MAXBUFF]; // type "jpg" or "png" or "gif"


static const char *progname = "alview";
static int aflag = 0;
static int lflag = 1;

/// this is cuz bamwhateverlib seems to depend on open output files which we do not care about
/// they are still there if we abend so zap dem old files
#if 0 // OLD - not used, used /dev/null for unix OR "nul" for windows instead 
void zap_old_okay2_files(void)
{
    DIR *d;
    struct dirent *e;
    time_t now;
    time_t timer;
    time_t ptimer;        // past
    char filename2zap[2048];
    char directory[2048];


    strcpy(directory, "/tmp/");

    timer = time(NULL);         // returns the time since 00:00:00 GMT, Jan. 1, 1970, measured in seconds.
//    ptimer = timer - (60 * 3 * 1);           // sec min hours -- cut off is 30 minutes
    ptimer = timer - (60 );           // sec min hours -- cut off is 30 minutes
    
    d = opendir(directory);
    if (! d) {
        perror(progname);
        return;
    }

#ifdef WIN32
#else 
    if (chdir(directory) == -1) {
        perror(progname);
        return;
    }
#endif

    now = time(NULL);
    
    while (1) 
    {
        e = readdir(d);
        if (! e) {
            break;
        }
        if (aflag || e->d_name[0] != '.') {
            if (lflag) 
            {
                struct stat s;
                if (stat(e->d_name, &s) == -1) {
                    fprintf(stderr, "%s: %s: ", progname, e->d_name);
                    perror(NULL);
                    continue;
                }
                           //123456789012
                if (strstr(e->d_name ,"okay2delete") != (char *)0)
                {
                    sprintf(filename2zap,"%s/%s",directory,e->d_name);
                    if (s.st_mtime < ptimer) 
                    {
                        printf("<!-- zapping %s -->\n",filename2zap); fflush(stdout);
#ifdef WIN32
                        _unlink(filename2zap);
#else
                        unlink(filename2zap);
#endif
                    }
                }
            }
        }
    }
    closedir(d);
}
#endif


static void zap_tmpalview_files(void)
{
    DIR *d;
    struct dirent *e;
    time_t now;
    time_t timer;
    time_t ptimer;        // past
    char filename2zap[2048];
    char directory[2048];


    strcpy(directory, HTDOCSDIR );

    timer = time(NULL);         // returns the time since 00:00:00 GMT, Jan. 1, 1970, measured in seconds.
//    ptimer = timer - (60 * 3 * 1);           // sec min hours -- cut off is 30 minutes
    ptimer = timer - (60 );           // sec min hours -- cut off is 30 minutes
    
    d = opendir(directory);
    if (! d) {
        perror(progname);
        return;
    }

#ifdef WIN32
#else
    if (chdir(directory) == -1) {
        perror(progname);
        return;
    }
#endif

    now = time(NULL);
    
    while (1) 
    {
        e = readdir(d);
        if (! e) {
            break;
        }
        if (aflag || e->d_name[0] != '.') {
            if (lflag) 
            {
                struct stat s;
                if (stat(e->d_name, &s) == -1) {
                    fprintf(stderr, "%s: %s: ", progname, e->d_name);
                    perror(NULL);
                    continue;
                }
                           //123456789012
                if (strstr(e->d_name ,"tmpalview") != (char *)0)
                {
                    sprintf(filename2zap,"%s/%s",directory,e->d_name);
                    if (s.st_mtime < ptimer) 
                    {
                        printf("<!-- zapping %s -->\n",filename2zap); fflush(stdout);
#ifdef WIN32
                        _unlink(filename2zap);
#else
                        unlink(filename2zap);
#endif
                    }
                }
            }
        }
    }
    closedir(d);
}




#ifdef WIN32
// OKAY, major windows nonsense here
//  fopen binary DOESN'T ALWAYS WORK ON SOME DRIVES !!!!!!!   --- NO.  it's a permissions thing.  icacls ?
#if 0 // dont need it after all
void write_windows_binary(char fn[], char *data, int datalen)
{
   HANDLE hFile = 0;
    DWORD dBytesWritten;

#if 0
    if ((hFile = CreateFile(LOG_FILE_PATH, GENERIC_WRITE, 0, NULL, 
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL)) == INVALID_HANDLE_VALUE
HANDLE WINAPI CreateFile(
  _In_      LPCTSTR lpFileName,
  _In_      DWORD dwDesiredAccess,
  _In_      DWORD dwShareMode,
  _In_opt_  LPSECURITY_ATTRIBUTES lpSecurityAttributes,
  _In_      DWORD dwCreationDisposition,
  _In_      DWORD dwFlagsAndAttributes,
  _In_opt_  HANDLE hTemplateFile
);
#endif


    size_t newsize = strlen(fn) + 1;
    wchar_t * wcstring = new wchar_t[newsize];
    size_t convertedChars = 0;
    mbstowcs_s(&convertedChars, wcstring, newsize, fn, _TRUNCATE);

    hFile = CreateFile(wcstring, GENERIC_WRITE|GENERIC_WRITE, FILE_SHARE_DELETE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) return;
    delete wcstring;
    if (!WriteFile(hFile, data, datalen, &dBytesWritten, NULL))
    { 
     // error            aff_error("WriteInLogFile");
    } 
    CloseHandle(hFile);
}
// 0
#endif
// windows
#endif


int savefile(char *outfn, struct image_type *im, int image_type)
{

// prototype next line ...
unsigned lodepng_encode24_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);

    lodepng_encode24_file(outfn, im->data, (unsigned) im->width, (unsigned) im->height);
    return 0;
}


#if 0
int spitfile(gdImagePtr im, int image_type)
{
  int i;
  int size;
  char *data;
char m[512];

// sprintf(m,"in spitfile im = %p  image_type = %d",im,image_type); jdebug(m); 
  size = 0;
  data = (char *)0;
  if (image_type == JPG_IMAGE)
      data = (char *) gdImageJpegPtr(im, &size, -1);
  else if (image_type == GIF_IMAGE)
      data = (char *) gdImageGifPtr(im, &size);
  else if (image_type == PNG_IMAGE)
  {
// sprintf(m,"in spitfile before gdImagePngPtr() "); jdebug(m); 
      data = (char *) gdImagePngPtr(im, &size);
// sprintf(m,"in spitfile after gdImagePngPtr() "); jdebug(m); 
  }
  else
  {
      printf("ERROR: invalid filetype %d <br>\n",image_type);
      return -1;
  }
// sprintf(m,"in spitfile after gdImagePngPtr() size = %d",size); jdebug(m); 
// sprintf(m,"in spitfile after gdImagePngPtr() data = %p",data); jdebug(m); 
sprintf(m,"in spitfile data = %p , size = %d ,before outputing ",data,size); jdebug(m);

  for (i=0 ; i<size ; i++)
  {
sprintf(m,"in spitfile data writing %d of %d ",i,size); jdebug(m);
      putchar(*(data+i));
  }
//  if (fwrite(data, 1, size, stdout) != size) 
  gdFree(data);  

sprintf(m,"end spitfile "); jdebug(m); 
  return 0;
}
#endif


int black, white, green, red, blue, orange, gray, yellow, purple,pink, lightgreen,lightgray, lightpink, lightblue, lightbrown,cyan;
int maroon;
int verylightgray;
int grays[10];
int grays2[10];
int schemered;
int schemecyan;
int schemeblue;
int schemegreen;
int schemeyellow;
int schemewhite;
int schemepurple;
int dna_a ;
int dna_c ;
int dna_g ;
int dna_t ;
char  dna_a_s[24];
char  dna_c_s[24];
char  dna_g_s[24];
char  dna_t_s[24];
char  dna_I_s[24]; // insert
char  dna_D_s[24]; // delete

void ImageColorAllocate(int *color,unsigned char R,unsigned char G,unsigned char B)
{
    unsigned char *z;
    z = (unsigned char *)color;
    *(z+0) = R;
    *(z+1) = G;
    *(z+2) = B;
    return;
}

void init_colors(void)
{
#if 0 // USE_GD 
    if (inited_colors_flag == 0) 
    {
        red   = gdImageColorAllocate(im, 255,    0,   0);
        green = gdImageColorAllocate(im,   0,  255,   0);
        blue =  gdImageColorAllocate(im,   0,    0, 255);
        cyan =  gdImageColorAllocate(im, 0x5, 0xeb,0xf9); // 05ebf9

        lightgreen = gdImageColorAllocate(im, 64, 125, 64);
        black = gdImageColorAllocate(im, 0, 0, 0);
        purple = gdImageColorAllocate(im, 150, 0, 255);
        orange = gdImageColorAllocate(im, 255, 100 , 0  );

// 9f005f ff5f00 bfff00 003fbf
strcpy(dna_a_s,"9f005f");  
strcpy(dna_c_s,"ff5f00");  
strcpy(dna_g_s,"bfff00");  
strcpy(dna_t_s,"003fbf"); 
strcpy(dna_I_s,"00ffff");  // cyan insert
strcpy(dna_D_s,"ffb000");  // yellow  delete

        dna_a = gdImageColorAllocate(im,0x9f,0x00,0x5f );
        dna_c = gdImageColorAllocate(im,0xff,0x5f,0x00 );
        dna_g = gdImageColorAllocate(im,0xbf,0xff,0x00 );
        dna_t = gdImageColorAllocate(im,0x00,0x3f,0xbf );

        pink = gdImageColorAllocate(im, 0xfd, 0x8b , 0xe9  ); // pink = gdImageColorAllocate(im, 200, 100 , 100  );
        lightblue = gdImageColorAllocate(im, 170, 170, 255);
        lightbrown = gdImageColorAllocate(im, 0xEE, 0xC5, 0x91); // A67D3D = 166 , 125, 61 //  EEC591
        white = gdImageColorAllocate(im, 255, 255, 255);
        gray = gdImageColorAllocate(im, 0xc0,0xc0,0xc0);
        lightgray = gdImageColorAllocate(im, 0xd0,0xd0,0xd0);
        lightpink = gdImageColorAllocate(im, 0xf1,0xcd,0xee); // f1cdee
        yellow = gdImageColorAllocate(im, 0xff, 0xb0 , 0x00  );//  fffb00
        maroon = gdImageColorAllocate(im, 0xd1, 0x34 , 0x00  );
        schemepurple =   gdImageColorAllocate(im, 0xea, 0xc7 , 0xeb  );
        schemewhite =   gdImageColorAllocate(im, 0xff, 0xff , 0xff  );
        schemered =    gdImageColorAllocate(im, 0xFF, 0x00 , 0x00  );
        schemecyan =   gdImageColorAllocate(im, 0x00, 0xff , 0xf0  );
        schemeblue =   gdImageColorAllocate(im, 0x00, 0xce , 0xf3  );
        schemegreen =  gdImageColorAllocate(im, 0x39, 0xe6 , 0x00  );
        schemeyellow = gdImageColorAllocate(im, 0xff, 0xb0 , 0x00  );//  fffb00

        verylightgray = gdImageColorAllocate(im, 0xf8,0xf8,0xf8);
        grays[0] = gdImageColorAllocate(im, 0xff,0xff,0xff);
        grays[1] = gdImageColorAllocate(im, 0xe0,0xe0,0xe0);
        grays[2] = gdImageColorAllocate(im, 0xc0,0xc0,0xc0);
        grays[3] = gdImageColorAllocate(im, 0xa0,0xa0,0xa0);
        grays[4] = gdImageColorAllocate(im, 0x80,0x80,0x80);

        grays[5] = gdImageColorAllocate(im, 0x60,0x60,0x60);
        grays[6] = gdImageColorAllocate(im, 0x40,0x40,0x40);
        grays[7] = gdImageColorAllocate(im, 0x20,0x20,0x20);
        grays[8] = gdImageColorAllocate(im, 0x10,0x10,0x10);
        grays[9] = gdImageColorAllocate(im, 0x00,0x00,0x00);

        grays2[0] = gdImageColorAllocate(im, 255,255,255);
        grays2[1] = gdImageColorAllocate(im, 245,245,245);
        grays2[2] = gdImageColorAllocate(im, 235,235,235);
        grays2[3] = gdImageColorAllocate(im, 225,225,225);
        grays2[4] = gdImageColorAllocate(im, 215,215,215);
        grays2[5] = gdImageColorAllocate(im, 205,205,205);
        grays2[6] = gdImageColorAllocate(im, 195,195,195);
        grays2[7] = gdImageColorAllocate(im, 185,185,185);
        grays2[8] = gdImageColorAllocate(im, 175,175,175);
        grays2[9] = gdImageColorAllocate(im, 165,165,165);

        inited_colors_flag = 1; 
    }
#endif

// our own custom implemention that does not use GD graphics library, but uses our custom implementation  ...
    if (inited_colors_flag == 0) 
    {
        ImageColorAllocate(&red, 255,    0,   0);
        ImageColorAllocate(&green,   0,  255,   0);
        ImageColorAllocate(&blue,   0,    0, 255);
        ImageColorAllocate(&cyan, 0x5, 0xeb,0xf9); // 05ebf9

        ImageColorAllocate( &lightgreen,64, 125, 64);
        ImageColorAllocate( &black,0, 0, 0);
        ImageColorAllocate( &purple,150, 0, 255);
        ImageColorAllocate( &orange, 255, 100 , 0  );

// 9f005f ff5f00 bfff00 003fbf
strcpy(dna_a_s,"9f005f");   // nucleotides acgt ...
strcpy(dna_c_s,"ff5f00");  
strcpy(dna_g_s,"bfff00");  
strcpy(dna_t_s,"003fbf"); 
strcpy(dna_I_s,"00ffff");  // cyan insert
strcpy(dna_D_s,"ffb000");  // yellow  delete

        ImageColorAllocate(&dna_a,0x9f,0x00,0x5f );
        ImageColorAllocate(&dna_c,0xff,0x5f,0x00 );
        ImageColorAllocate(&dna_g,0xbf,0xff,0x00 );
        ImageColorAllocate(&dna_t,0x00,0x3f,0xbf );

        ImageColorAllocate(&pink, 0xfd, 0x8b , 0xe9  ); // pink = ImageColorAllocate( 200, 100 , 100  );
        ImageColorAllocate(&lightblue, 170, 170, 255);
        ImageColorAllocate(&lightbrown, 0xEE, 0xC5, 0x91); // A67D3D = 166 , 125, 61 //  EEC591
        ImageColorAllocate(& white, 255, 255, 255);
        ImageColorAllocate(& gray, 0xc0,0xc0,0xc0);
        ImageColorAllocate(& lightgray, 0xd0,0xd0,0xd0);
        ImageColorAllocate(& lightpink, 0xf1,0xcd,0xee); // f1cdee
        ImageColorAllocate(& yellow, 0xff, 0xb0 , 0x00  );//  fffb00
        ImageColorAllocate(& maroon, 0xd1, 0x34 , 0x00  );
        ImageColorAllocate(& schemepurple, 0xea, 0xc7 , 0xeb  );
        ImageColorAllocate(&schemewhite , 0xff, 0xff , 0xff  );
        ImageColorAllocate(&schemered  ,0xFF, 0x00 , 0x00  );
        ImageColorAllocate(&schemecyan  ,0x00, 0xff , 0xf0  );
        ImageColorAllocate(&schemeblue  ,0x00, 0xce , 0xf3  );
        ImageColorAllocate(&schemegreen  ,0x39, 0xe6 , 0x00  );
        ImageColorAllocate(&schemeyellow  ,0xff, 0xb0 , 0x00  );//  fffb00

        ImageColorAllocate(&verylightgray,  0xf8,0xf8,0xf8);
        ImageColorAllocate(& grays[0] , 0xff,0xff,0xff);
        ImageColorAllocate(& grays[1] , 0xe0,0xe0,0xe0);
        ImageColorAllocate(& grays[2] , 0xc0,0xc0,0xc0);
        ImageColorAllocate(& grays[3] , 0xa0,0xa0,0xa0);
        ImageColorAllocate(& grays[4] , 0x80,0x80,0x80);
        ImageColorAllocate(& grays[5] , 0x60,0x60,0x60);
        ImageColorAllocate(& grays[6] , 0x40,0x40,0x40);
        ImageColorAllocate(& grays[7] , 0x20,0x20,0x20);
        ImageColorAllocate(& grays[8] , 0x10,0x10,0x10);
        ImageColorAllocate(& grays[9] , 0x00,0x00,0x00);

        ImageColorAllocate(&grays2[0],  255,255,255);
        ImageColorAllocate(&grays2[1],  245,245,245);
        ImageColorAllocate(&grays2[2],  235,235,235);
        ImageColorAllocate(&grays2[3],  225,225,225);
        ImageColorAllocate(&grays2[4],  215,215,215);
        ImageColorAllocate(&grays2[5],  205,205,205);
        ImageColorAllocate(&grays2[6],  195,195,195);
        ImageColorAllocate(&grays2[7],  185,185,185);
        ImageColorAllocate(&grays2[8],  175,175,175);
        ImageColorAllocate(&grays2[9],  165,165,165);

        inited_colors_flag = 1; 
    }
    return;
}



#if 0
void dumpgif(char fn[])
{
    int c;
    FILE *fp;

    fp = fopen(fn,"rb"); // debug routine
    while ((c = getc(fp) ) != -1 )
    {
        putchar(c);
    }
    fclose(fp);
    return;
}
#endif



#if 0
unsigned int chromoffsets[50];
int numchromoffsets;

int load_chroffsets(char *fn)
{
    int i;
    char m[512];
    int status;
    FILE *fp;


sprintf(m,"Loading chr offsets from %s",fn); jdebug(m); 
    fp = fopen(fn,"rb");
    if (fp == (FILE *)0)
    {
        fprintf(stderr,"Can not OPEN file %s \n",fn);
        return -1;
    }
    numchromoffsets = 0;
    while (1)
    {
        status = fread(chromoffsets+numchromoffsets,4,1,fp);
        if (status != 1)
            break;
        numchromoffsets++;
    }

    fclose(fp);
sprintf(m,"Loaded %d chrs from %s",numchromoffsets,fn); jdebug(m);

#if 1 // HACK
    for (i=0;i<numchromoffsets;i++)
        *(chromoffsets+i) = *(chromoffsets+i) + i;
#endif
    return 0;
}
#endif


#if 0
void ui2genomiccoords(unsigned int ui, char puthere[])
{
    int i;
    int lo = 0;
    for (i=0;i<(numchromoffsets-1);i++)
    {
        if ( (ui < *(chromoffsets+i))  && (ui >= lo) )
        {
             sprintf(puthere,"chr%d:%u",i+1,ui-lo);
             return;
        }
        lo = *(chromoffsets+i) ;
    }
    strcpy(puthere,"err:0");
    return;
}
#endif


int get_chr_spot(char chr[],int *spot_ptr,int *end_ptr,int *chrom_index_ptr)
{
    int k = -1;

    *chrom_index_ptr = k;
    if (strcmp("X",chr) == 0) k = 22;
    else if (strcmp("Y",chr) == 0) k = 23;
    else if (strcmp("M",chr) == 0) k = 24;
    else if (strcmp("chrX",chr) == 0) k = 22;
    else if (strcmp("chrY",chr) == 0) k = 23;
    else if (strcmp("chrM",chr) == 0) k = 24;
    else
    {
        if (strncmp(chr,"chr",3) == 0) k = atoi(&chr[3]); 
        else k = atoi(&chr[0]); 
        k = k - 1;
    }
    if (k>=0) 
    {
        *chrom_index_ptr = k;
        return 0;
    }
fprintf(stderr,"in get_chr_spot for [%s]chr is %d\n",chr,-1);  fflush(stderr);  
    return -1;
}


int get_chr_lo_hi(char chr[],int *start_ptr,int *end_ptr,int *chrom_index_ptr)
{
    int k = -1;

    *chrom_index_ptr = k;
    if (strcmp("X",chr) == 0) k = 22;
    else if (strcmp("Y",chr) == 0) k = 23;
    else if (strcmp("M",chr) == 0) k = 24;
    else if (strcmp("chrX",chr) == 0) k = 22;
    else if (strcmp("chrY",chr) == 0) k = 23;
    else if (strcmp("chrM",chr) == 0) k = 24;
    else
    {
        if (strncmp(chr,"chr",3) == 0) k = atoi(&chr[3]); 
        else k = atoi(&chr[0]); 
        k = k - 1;
    }
    if (k>=0) 
    {
        *chrom_index_ptr = k;
        return 0;
    }
  
    return -1;

// sprintf(m,"in get_chr_lo_hi() will search for khrom = %d (atoi()ed [and -1] from string %s)",k,chr); jdebug(m);

#if 0
char m[512];
    int i;
    int lo = 0;
    *chrom_index_ptr = k;
    for (i=0;i<(numchromoffsets);i++)
    {
        if (k == i)
        {
//             *start_ptr = *start_ptr + lo; // *(chromoffsets+i);
//             *end_ptr = *end_ptr + lo; // *(chromoffsets+i);
             return 0;
        }
//        lo = *(chromoffsets+i);
    }
sprintf(m,"ERROR: in get_chr_lo_hi() bad input "); jdebug(m);
    return -1;
#endif
}




#if 0
             // get and parse this -> "33:T>C,34:A>T"
int getbowtie_places(long int foff, int *i1, char *c1,char *c2, int *i2, char *c3, char *c4)
{
    int tokens = 0;
    int i;
    char s[1024];
     char m[1024];

// sprintf(m,"in getbowtie_places %ld %p",foff,fp_cigar); jdebug(m);
    if (fseek(fp_cigar,(size_t)foff,SEEK_SET) != 0)
    {
        sprintf(m,"ERROR in getbowtie_places(), fseek(%p) to %ld failed ",fp_cigar,foff); 
        jdebug(m); 
        return -1;
    }
    s[0] = (char)0;
    fgets(s,1000,fp_cigar);

    for (i=0 ; s[i] ; i++) 
    {
        if (s[i] == ':') s[i] = ' ';
        if (s[i] == '>') s[i] = ' ';
        if (s[i] == ',') s[i] = ' ';
    }
    tokens = sscanf(s,"%d %c %c %d %c %c",i1,c1,c2,i2,c3,c4);
    if (tokens == 6) return 2;
    return 1;
}
#endif


static int line[5000];
static int line2[5000];
static int mms[5000];
static int mmsLQ[5000];
static double pxwidth; // pixel width 

#define SLOTS 0
// slots is a new method - still working on this - coming soon
#if SLOTS
unsigned char *slots;   // slots is this y dimension is slots, x is pixels, set to one for is already covered.
int num_slots;
#else
unsigned char *ppp;
int ppp_firsttime = 1;
#endif

static int set_which_mm(int whichmm) 
{
    int k = orange; // error, should not get this
    char m[512];
// sprintf(m,"in whichmm , whichmm = %d ",whichmm);  jdebug(m);

    if (whichmm == 0)      k = dna_a;                    // a
    else if (whichmm == 1) k = dna_c;                    // c
    else if (whichmm == 2) k = dna_g;                    // g
    else if (whichmm == 3) k = dna_t;                    // t
    else
    {
sprintf(m,"ERROR: invalid param to set_which_mm = %d ",whichmm); jdebug(m);
    }
    return k;
}



#define MAX_CIG_OPS 50

struct cig_type
{
    char cmd;
    int len;
    int sofar;
} cig_ops[MAX_CIG_OPS];
int num_cig_ops = 0;


static void setup_cigar_stuff(char *s)
{
/*
A CIGAR string is comprised of a series of operation lengths plus the operations. The conventional CIGAR format allows
for three types of operations: M for match or mismatch, I for insertion and D for deletion. The extended CIGAR format
further allows four more operations, as is shown in the following table, to describe clipping, padding and splicing:
op Description
M Alignment match (can be a sequence match or mismatch)
I Insertion to the reference
D Deletion from the reference
N Skipped region from the reference
S Soft clip on the read (clipped sequence present in <seq>)
H Hard clip on the read (clipped sequence NOT present in <seq>)
P Padding (silent deletion from the padded reference sequence)
2.2.4. Format of optional fields
*/
    char nums[512];
    int i;
    int numi = 0;

    num_cig_ops = 0;
    i = 0;
    nums[0] = (char)0;
    while (1) 
    {
        if (s[i] == (char)0) break;
        if (num_cig_ops > MAX_CIG_OPS) break;
        if (isdigit(s[i])) { nums[numi++] = s[i]; nums[numi] = (char)0; }
        else 
        {
            cig_ops[num_cig_ops].cmd = s[i];
            cig_ops[num_cig_ops].len = atoi(nums);
            cig_ops[num_cig_ops].sofar = 0;
            num_cig_ops++;
            nums[0] = (char)0; numi = 0;
        }
        i++;
    }
#if 0
 char m[512];
 for (i=0;i<num_cig_ops;i++) { sprintf(m,"cig %d %c %d ",i,cig_ops[i].cmd, cig_ops[i].len); jdebug(m); }
#endif
    return;
}


static int extra_splice_space(void) 
{
    int i;
    int ret = 0;
    for (i=0;i<num_cig_ops;i++) 
    {
        if (cig_ops[i].cmd == 'N') ret += cig_ops[i].len; 
    }
    return ret;
}


long int time_out_every_once_a_while = 0;

static void aldetails(int diff, int offset, int len,  
        int fwdflag, int optical_dupe_flag ,char seq[], char *quality, int seqlen, char cigar[], int splice_source)
{
    char ch;
    int oset;
    int j,k,coff;
    int dnaat = 0;
    int kkolor;
    int gotn;
    double d;
    int x,xx;
    int jump = 1;
    int xcolor;
    int frcolor;
    int i;
    int x1,x2,y1,y2;
    int keepgoing = 1;

    if (!dnaspace) return;
    if (offset < 0) return;

    if (spliceonly_flag == 1)
    {
        for (gotn=i=0 ; cigar[i] ; i++)
        {
           if (cigar[i] == 'N')
           {
               gotn = 1;
               break;
           }
        }
       if (gotn == 0) return; 
    }


    dnaat = offset;

    setup_cigar_stuff(cigar);
    x1 = x2 = y1 = y2 = 0;
    frcolor = purple;


// fprintf(stderr,"in aldetails pxw=%f offset=%d",pxwidth,offset); fflush(stderr); 

#ifdef QT_GUI
 if ((++time_out_every_once_a_while % 50000) == 0)
 {
 // might need to cal this QCoreApplication::processEvents()QCoreApplication::processEvents()
 hack_timeout_to_process_events();
 time_out_every_once_a_while = 1; // reset count
 }
#endif

#if SLOTS
#else
// xxx
    if (ppp_firsttime == 1) 
    {
        if (ppp) {free(ppp); ppp = (unsigned char *)0; }
        ppp = (unsigned char *)malloc((iw*ih) * sizeof(unsigned char)); // ppp 
        if (ppp == (void *)0) 
        {
            fprintf(stderr,"ERROR- can't malloc space for image.  (bad news, we're out of memory) in aldetails()\n"); 
            fflush(stderr);
            exit(0); 
        }
        memset(ppp,0,(iw*ih) * sizeof(unsigned char));
        ppp_firsttime = 0;
    }

// FIND A PLACE TO PUT THE READ (alignment)
    if (cig_ops[0].cmd == 'S') // soft clip
    {
        offset  = offset - cig_ops[0].len;
    }
    jump = 1;
    while (jump < (ih-60) ) 
    {
        x1 = (int)(pxwidth * (double)offset); 
        x2 = x1 + (int)(pxwidth * (double)len); 
        x2 = (int)(x2 + extra_splice_space() * pxwidth);  // count the "N" (splice) in cigar
        if ( x2 > iw) x2=iw; 

        y1 = ih - 50 - jump;
        y2 = y1 - 3;
        if ((y2<=0)|| (y1<=0)) { y2 = 1; y1=4; break; }
        if (x1 < 0) x1 = 0;
        if (x2 < 0) x2 = 0;
        if (y1 < 0) y1 = 0;
        if (y2 < 0) y2 = 0;
        for (keepgoing = 0, i=x1 ; i<=(x2+1) ; i++) 
        {
             oset = (y2*iw)+i;
             if (oset < (iw*ih))
             {
                 if (*(ppp+oset) == 1) //  if (*(ppp[i][y2] == 1)
                 {   
                     jump += 5;
                     keepgoing = 1;
                     break;
                 }
             }
        }
        if (keepgoing == 0) break;
    }
  
    for (i=x1;i<=(x2+1);i++) 
    {
        oset = (y2*iw)+i;
        if (oset < (iw*ih))
        {
            *(ppp+oset) = (char)1;        // ppp[i][y2] = 1;
            // *(ppp+(y2*iw)+i) = (char)1;        // ppp[i][y2] = 1;
        }
    }
#endif

// sprintf(m,"in aldetailx x1= %d pxw=%f",x1,pxwidth); jdebug(m); 

    if (fwdflag == 0) xcolor = grays[2]; else  xcolor = grays[3];
    for (i=k=coff=j=0 ; j<num_cig_ops ; j++) 
    {
        if (i>= seqlen) break; // done
        if (cig_ops[j].cmd == 'N')  // skip 
        {
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            kkolor = yellow;
            if (splice_source == 1)  kkolor = black;  // refseq
            else if (splice_source == 2)  kkolor = red; // novel altsplice
            else if (splice_source == 3)  kkolor = green; // est
            ImageFilledRectangle(im,x,y2+1,xx,y1-2,kkolor); // make a thin line twixt alignments
            coff += cig_ops[j].len; 
            dnaat += cig_ops[j].len; 
            continue; 
        }
        else if (cig_ops[j].cmd == 'D')  // deletion 
        {
            dcnt++; // "d)elete count 
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            kkolor = yellow;
#if 0 // USE_GD
            gdImageFilledRectangle(im,x,y2-1,xx,y1+2,kkolor);
#endif
            ImageFilledRectangle(im,x,y2-1,xx,y1+2,kkolor);

            coff += cig_ops[j].len; 
            dnaat += cig_ops[j].len; 
            continue;
        }
        else if (cig_ops[j].cmd == 'I') // insert  
        {
            icnt++; // "i)nsert count 
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            xcolor = cyan; // mark insert here
#if 0 // USE_GD
            gdImageFilledRectangle(im,x,y2-2,xx,y1+2,xcolor);
#endif
            ImageFilledRectangle(im,x,y2-2,xx,y1+2,xcolor);

            coff -= cig_ops[j].len; 

            if ((dnaat>=0)&&(dnaat<dnasize))
            {
// potential integer overflow  - fix ?
                *(dnacnts+dnaat) += 1;
                *(dnamms+dnaat) += 1;
            }
            dnaat -= cig_ops[j].len; 
            continue;
        }
        else if (cig_ops[j].cmd == 'S') // soft clip
        {
            for (k=0 ; k<cig_ops[j].len ; k++)
            {
                xcolor = black;
                d = x1 + (pxwidth * (i+coff));
                x = (int)d;
                d = x1 + (pxwidth * (i+1+coff));
                x1 = (int)(pxwidth * (double)(offset+i+coff)); 
                x2 = (int)(x1 + pxwidth);
                x = x1;
                xx = x2;
                i++;
                xcolor = black;
#if 0 // USE_GD
                gdImageFilledRectangle(im,x,y2,xx,y1,xcolor);
#endif
                ImageFilledRectangle(im,x,y2,xx,y1,xcolor);
            }
            // dnaat += cig_ops[j].len; 
            continue;
        }
        else if (cig_ops[j].cmd == 'H') // ignore !
        {
// ignore
        }

            /* SO it's draw nucleotide "matches" = 'M' ------------------*/
        for (k=0 ; k<cig_ops[j].len ; k++)
        {
// fprintf(stderr,"in aldetails(), dnaat=%d\n", dnaat); fflush(stderr); 
            if (dnaat >= dnasize) break;
            ch = seq[i];
// fprintf(stderr,"in aldetails(), ch=%c\n", ch); fflush(stderr); 
            if (fwdflag == 0) xcolor = grays[2]; else  xcolor = grays[3];
            if ((dnaat >= 0) && (dnaat < dnasize))
            {
                if ((basecolors_flag == 1) || ( (toupper(ch)) != toupper(*(dnaspace+dnaat)) ))
                {
                    if      (ch == 'A') xcolor = dna_a; else if (ch == 'C') xcolor = dna_c;
                    else if (ch == 'G') xcolor = dna_g; else if (ch == 'T') xcolor = dna_t;
                }
                d = x1 + (pxwidth * (i+coff));
                x = (int)d;
                d = x1 + (pxwidth * (i+1+coff));
                x1 = (int)(pxwidth * (double)(offset+i+coff)); 
                x2 = (int)(x1 + pxwidth);
                x = x1;
                xx = x2;
                if (
                     (x>=0)  && (x  < iw) &&
                     (xx>=0) && (xx < iw) &&
                     (y2>=0) && (y2 < ih) &&
                     (y1>=0) && (y1 < ih) 
                   )
                {
#if 1
           if (qual_flag == 1)
           {                // override if qual flag on
// quality values co-incident with sequence 
               char qch;
                if (quality) qch = *(quality+i);
                else         qch = 0;
#if 0
/*
The character '!' represents the lowest quality while '~' is the highest. Here are the quality value characters in left-to-right increasing order of quality (ASCII):
  Vdoublequote
 !.#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz
 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
 0         1         2         3         4         5         6         7         8         9
*/
#endif
            int idx = qch;
// sprintf(m,"idx = %d from %c",idx,qch); jdebug(m);
// grays[0] = white
// grays[9] = black

            if (idx > 65)      xcolor = grays[1];
            else if (idx > 61) xcolor = grays[2];
            else if (idx > 57) xcolor = grays[3];
            else if (idx > 53) xcolor = grays[4];
            else if (idx > 49) xcolor = grays[5];
            else if (idx > 45) xcolor = grays[6];
            else if (idx > 41) xcolor = grays[7];
            else if (idx > 37) xcolor = grays[8];
            else               xcolor = grays[9];
            }
#endif
#if 0 // USE_GD
                    gdImageFilledRectangle(im,x,y2,xx,y1,xcolor);
#endif
                    ImageFilledRectangle(im,x,y2,xx,y1,xcolor);
// debug fprintf(stderr,"in aldetails pxw=%f offset=%d, before ImagFillRectangle %d %d %d %d",pxwidth,offset, x,y2,xx,y1); fflush(stderr);

                }
                if (dnacnts) *(dnacnts+dnaat) += 1;
                if ((dnaspace) && (dnamms))
                {
                    if ( (toupper(ch)) != toupper(*(dnaspace+dnaat)) ) { *(dnamms+dnaat) += 1; }
                }
            }
            i++;
            dnaat += 1;
        }
    }
// make a notch at start (fwd or rev) for read 
#if 0
    if (0) // pxwidth < 5
    {
        if (fwdflag == 0)
        {
            x1 = (int)(pxwidth * (double)offset); 
#if 0 // USE_GD
            gdImageFilledRectangle(im,x1  ,y2+1,x1+1,y1-1,yellow);
#endif
            ImageFilledRectangle(im,x1  ,y2+1,x1+1,y1-1,yellow);

        }
        else
#if 0 // USE_GD
            gdImageFilledRectangle(im,x2-2,y2+1,x2,  y1-1,yellow);
#endif
            ImageFilledRectangle(im,x2-2,y2+1,x2,  y1-1,yellow);

    }
#endif
    return;
}


void snp_call_aldetails(int diff, int offset, int len,  
        int fwdflag, int optical_dupe_flag ,char seq[], char *quality, int seqlen, char cigar[], int splice_source)
{
    char ch;
    int oset;
    int j,k,coff;
    int dnaat = 0;
    int kkolor;
    int gotn;
    double d;
    int x,xx;
    int jump = 1;
    int xcolor;
    int frcolor;
    int i;
    int x1,x2,y1,y2;
    int keepgoing = 1;


// debug fprintf(stderr,"in snp_call_aldetails(), snp_call_flag = %d %s offset=%d %p\n", snp_call_flag,chr,offset,dnaspace); 
fflush(stderr); 

    if (!dnaspace) return;
    if (offset < 0) return;

    if (spliceonly_flag == 1)
    {
        for (gotn=i=0 ; cigar[i] ; i++)
        {
           if (cigar[i] == 'N')
           {
               gotn = 1;
               break;
           }
        }
       if (gotn == 0) return; 
    }


    dnaat = offset;

    setup_cigar_stuff(cigar);
    x1 = x2 = y1 = y2 = 0;
    frcolor = purple;


// fprintf(stderr,"in snp_call_aldetails pxw=%f offset=%d",pxwidth,offset); fflush(stderr); 

#ifdef QT_GUI
 if ((++time_out_every_once_a_while % 50000) == 0)
 {
 // might need to cal this QCoreApplication::processEvents()QCoreApplication::processEvents()
 hack_timeout_to_process_events();
 time_out_every_once_a_while = 1; // reset count
 }
#endif

    if (ppp_firsttime == 1) 
    {
        if (ppp) {free(ppp); ppp = (unsigned char *)0; }
        ppp = (unsigned char *)malloc((iw*ih) * sizeof(unsigned char)); // ppp 
        if (ppp == (void *)0) 
        {
            fprintf(stderr,"ERROR- can't malloc space for image.  (bad news, we're out of memory) in snp_call_aldetails()\n"); 
            fflush(stderr);
            exit(0); 
        }
        memset(ppp,0,(iw*ih) * sizeof(unsigned char));
        ppp_firsttime = 0;
    }

// FIND A PLACE TO PUT THE READ(alignment)
    if (cig_ops[0].cmd == 'S') // soft clip
    {
        offset  = offset - cig_ops[0].len;
    }
    jump = 1;
    while (jump < (ih-60) ) 
    {
        x1 = (int)(pxwidth * (double)offset); 
        x2 = x1 + (int)(pxwidth * (double)len); 
        x2 = (int)(x2 + extra_splice_space() * pxwidth); 
        if ( x2 > iw) x2=iw; 

        y1 = ih - 50 - jump;
        y2 = y1 - 3;
        if ((y2<=0)|| (y1<=0)) { y2 = 1; y1=4; break; }
        if (x1 < 0) x1 = 0;
        if (x2 < 0) x2 = 0;
        if (y1 < 0) y1 = 0;
        if (y2 < 0) y2 = 0;
        for (keepgoing = 0, i=x1 ; i<=(x2+1) ; i++) 
        {
             oset = (y2*iw)+i;
             if (oset < (iw*ih))
             {
                 if (*(ppp+oset) == 1) //  if (*(ppp[i][y2] == 1)
                 {   
                     jump += 5;
                     keepgoing = 1;
                     break;
                 }
             }
        }
        if (keepgoing == 0) break;
    }
  
    for (i=x1;i<=(x2+1);i++) 
    {
        oset = (y2*iw)+i;
        if (oset < (iw*ih))
        {
            *(ppp+oset) = (char)1;        // ppp[i][y2] = 1;
            // *(ppp+(y2*iw)+i) = (char)1;        // ppp[i][y2] = 1;
        }
    }

    if (fwdflag == 0) xcolor = grays[2]; else  xcolor = grays[3];
    for (i=k=coff=j=0 ; j<num_cig_ops ; j++) 
    {
        if (i>= seqlen) break; // done
        if (cig_ops[j].cmd == 'N')  // skip 
        {
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            kkolor = yellow;
            if (splice_source == 1)  kkolor = black;  // refseq
            else if (splice_source == 2)  kkolor = red; // novel altsplice
            else if (splice_source == 3)  kkolor = green; // est
//             ImageFilledRectangle(im,x,y2+1,xx,y1-2,kkolor); // make a thin line twixt alignments
            coff += cig_ops[j].len; 
            dnaat += cig_ops[j].len; 
            continue; 
        }
        else if (cig_ops[j].cmd == 'D')  // deletion 
        {
            dcnt++; // "d)elete count 
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            kkolor = yellow;
//             ImageFilledRectangle(im,x,y2-1,xx,y1+2,kkolor);
            coff += cig_ops[j].len; 
            dnaat += cig_ops[j].len; 
#if CMD_LINE
            if (snp_call_flag == 1) if (dnaat == snp_call_dnaat) snp_call_Del_cnt++;
#endif
            continue;
        }
        else if (cig_ops[j].cmd == 'I') // insert  
        {
            icnt++; // "i)nsert count 
            d = x1 + (pxwidth * (i+coff));
            x = (int)d;
            d = x1 + (pxwidth * (i+1+coff));
            x1 = (int)(pxwidth * (double)(offset+i+coff)); 
            x2 = (int)(pxwidth * (double)(offset+i+coff+cig_ops[j].len)); 
            x = x1;
            xx = x2;
            xcolor = cyan; // mark insert here
//             ImageFilledRectangle(im,x,y2-2,xx,y1+2,xcolor);
            coff -= cig_ops[j].len; 

            if ((dnaat>=0)&&(dnaat<dnasize))
            {
// potential integer overflow  - fix ?
#if CMD_LINE
                if (snp_call_flag == 1) if (dnaat == snp_call_dnaat) snp_call_Ins_cnt++;
#endif

                *(dnacnts+dnaat) += 1;
                *(dnamms+dnaat) += 1;
            }
            dnaat -= cig_ops[j].len; 
            continue;
        }
        else if (cig_ops[j].cmd == 'S') // soft clip
        {
            for (k=0 ; k<cig_ops[j].len ; k++)
            {
                xcolor = black;
                d = x1 + (pxwidth * (i+coff));
                x = (int)d;
                d = x1 + (pxwidth * (i+1+coff));
                x1 = (int)(pxwidth * (double)(offset+i+coff)); 
                x2 = (int)(x1 + pxwidth);
                x = x1;
                xx = x2;
                i++;
                xcolor = black;
            }
            continue;
        }
        else if (cig_ops[j].cmd == 'H') // ignore ! hard clip
        {                    // ignore
        }

            /* SO it's draw nucleotide "matches" = 'M' ------------------*/
        for (k=0 ; k<cig_ops[j].len ; k++)
        {
            if (dnaat >= dnasize) break;
            ch = seq[i];
// fprintf(stderr,"in snp_call_aldetails(), dnaat=%d snd=%d '%c'\n", dnaat,snp_call_dnaat,ch); fflush(stderr); 
            if (fwdflag == 0) xcolor = grays[2]; else  xcolor = grays[3];
            if ((dnaat >= 0) && (dnaat < dnasize))
            {
                if ((basecolors_flag == 1) || ( (toupper(ch)) != toupper(*(dnaspace+dnaat)) ))
                {
                    if      (ch == 'A') xcolor = dna_a; else if (ch == 'C') xcolor = dna_c;
                    else if (ch == 'G') xcolor = dna_g; else if (ch == 'T') xcolor = dna_t;
                }
                d = x1 + (pxwidth * (i+coff));
                x = (int)d;
                d = x1 + (pxwidth * (i+1+coff));
                x1 = (int)(pxwidth * (double)(offset+i+coff)); 
                x2 = (int)(x1 + pxwidth);
                x = x1;
                xx = x2;
                if (
                     (x>=0)  && (x  < iw) &&
                     (xx>=0) && (xx < iw) &&
                     (y2>=0) && (y2 < ih) &&
                     (y1>=0) && (y1 < ih) 
                   )
                {
#if 1
           if (qual_flag == 1)
           {                // override if qual flag on
// quality values co-incident with sequence 
               char qch;
                if (quality) qch = *(quality+i);
                else         qch = 0;
#if 0
/*
The character '!' represents the lowest quality while '~' is the highest. Here are the quality value characters in left-to-right increasing order of quality (ASCII):
  Vdoublequote
 !.#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz
 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
 0         1         2         3         4         5         6         7         8         9
*/
#endif
            int idx = qch;
// sprintf(m,"idx = %d from %c",idx,qch); jdebug(m);
// grays[0] = white
// grays[9] = black

            if (idx > 65)      xcolor = grays[1];
            else if (idx > 61) xcolor = grays[2];
            else if (idx > 57) xcolor = grays[3];
            else if (idx > 53) xcolor = grays[4];
            else if (idx > 49) xcolor = grays[5];
            else if (idx > 45) xcolor = grays[6];
            else if (idx > 41) xcolor = grays[7];
            else if (idx > 37) xcolor = grays[8];
            else               xcolor = grays[9];
            }
#endif
// fprintf(stderr,"in snp_call_aldetails pxw=%f offset=%d, before ImagFillRectangle %d %d %d %d",pxwidth,offset, x,y2,xx,y1); fflush(stderr);
                //     ImageFilledRectangle(im,x,y2,xx,y1,xcolor);
                }

                if (dnacnts) *(dnacnts+dnaat) += 1;
                if ((dnaspace) && (dnamms))
                {
                    if ( (toupper(ch)) != toupper(*(dnaspace+dnaat)) ) { *(dnamms+dnaat) += 1; }
                }
#if CMD_LINE
            if ((snp_call_flag == 1) && (dnaat == snp_call_dnaat))
            {
                switch (toupper(ch) ) {
                    case 'A': snp_call_A_cnt++; break;
                    case 'C': snp_call_C_cnt++; break;
                    case 'G': snp_call_G_cnt++; break;
                    case 'T': snp_call_T_cnt++; break;
                    case 'N':                   break;
                    default: fprintf(stderr, "ERROR, invalid toupper in snp_call_aldetail for snp_call\n"); 
                }
                if (dnaspace) snp_call_referece = toupper(*(dnaspace+dnaat));
            }
#endif
            }
            i++;
            dnaat += 1;
        }
    }
    return;
}


// WHY ?
static void paint_align_box(int diff, int offset, int len, int wcode, int notch_mm_where, int qual, 
        int fwdflag, int bowtie_flag,int split_flag, int whichmm1, int whichmm2, int optical_dupe_flag )
{
    int jump = 1;
    int xcolor;
    int frcolor;
    int i;
    int x1,x2,y1,y2;
    int mmx1, mmx2;
    int mmcolor;
    int mm_min = -1; 
    int mm_max = -1; 
    int mm2_min = -1; 
    int mm2_max = -1; 
    int mm_minLQ = -1; 
    int mm_maxLQ = -1; 
    int keepgoing = 1;

 
    x1 = x2 = y1 = y2 = 0;
    frcolor = purple;

    if (ppp_firsttime == 1) 
    {
        ppp = (unsigned char *)malloc(((iw*ih) * sizeof(unsigned char)) + iw); // ppp (with some padding)
        if (ppp == (void *)0) 
        {
            printf("ERROR- can't malloc space for image.  (bad news, we're out of memory) in paint_align_box()\n"); 
            exit(0); 
        }
        memset(ppp,0,(iw*ih) * sizeof(unsigned char));
        ppp_firsttime = 0;
    }
    jump = 1;
    while (jump < (ih-60) ) 
    {
        x1 = (int)(pxwidth * (double)offset); 
        x2 = x1 + (int)(pxwidth * (double)len); 
        y1 = ih - 50 - jump;
        y2 = y1 - 3;
if ((y2<=0)|| (y1<=0)) { y2 = 1; y1=4; break; }
        if (x1<0) x1 = 0;
        if (x2<0) x2 = 0;
        if (y1<0) y1 = 0;
        if (y2<0) y2 = 0;
        for (keepgoing = 0, i=x1 ; i<=(x2+1) ; i++) 
        {
             if (*(ppp+(y2*iw)+i) == 1) //  if (*(ppp[i][y2] == 1)
             {   
                 jump += 5;
                 keepgoing = 1;
                 break;
             }
        }
// printf("i = %d , y1 = %d , jump = %d , keepgoing = %d ",i,y1,jump,keepgoing);fflush(stdout); 
        if (keepgoing == 0) break;
    }

  
    for (i=x1;i<=(x2+1);i++) 
    {
        if ( ((y2+iw)+i) < (iw*ih) ) *(ppp+(y2*iw)+i) = (char)1;        // ppp[i][y2] = 1;
    }
    if (wcode == 0)      xcolor = gray;        // pm1 
    else if (wcode == 1) xcolor = gray;        // pm1 
    else if (wcode == 2) xcolor = lightpink;   // pm++
    else if (wcode == 3) xcolor = green;       // mm1
    else if (wcode == 4) xcolor = lightgreen;  // mm++
    else                xcolor = red;
    if (bowtie_flag == 1) 
    {
       if (split_flag == 1) 
           xcolor = lightbrown;
       else
           xcolor = lightblue;
    }
    if (optical_dupe_flag == 1) 
    {
        xcolor = lightpink;
    }
    
#if 0 // USE_GD
    gdImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
#endif
    ImageFilledRectangle(im,x1,y2,x2,y1,xcolor);

                      // make a notch at start (fwd or rev) for read 
    if (fwdflag == 1)
#if 0 // USE_GD
        gdImageFilledRectangle(im,x1  ,y2+1,x1+1,y1-1,lightpink);
#endif
        ImageFilledRectangle(im,x1  ,y2+1,x1+1,y1-1,lightpink);

    else
#if 0 // USE_GD
        gdImageFilledRectangle(im,x2-1,y2+1,x2,  y1-1,lightpink);
#endif
        ImageFilledRectangle(im,x2-1,y2+1,x2,  y1-1,lightpink);

    if (bowtie_flag) 
    {
        if (bi1>=0)
        {
            notch_mm_where = bi1;
            mmx1 = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mmx2 = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
            mm_min = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mm_max = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
            mmcolor = set_which_mm(whichmm1); 
#if 0 // USE_GD
            gdImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);
#endif
            ImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);

        }
        if (bi2>=0)
        {
            notch_mm_where = bi2;
            mmx1 = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mmx2 = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
        mm2_min = (int)(pxwidth * (double)(offset+notch_mm_where)); 
        mm2_max = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
// nope! if (qual < QUALCUTOFF ) mmcolor = pink; else mmcolor = orange;
            mmcolor = set_which_mm(whichmm2); 
#if 0 // USE_GD
            gdImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);
#endif
            ImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);

        }
        for (i=x1 ;i<=x2 ; i++)
        {
            if ((wcode == 1) || (wcode == 3)) // single hit
                *(line+i) = *(line+i) + 1;
            *(line2+i) = *(line2+i) + 1;
            if ((i >= mm_min) && (i <= mm_max))
                *(mms+i) = *(mms+i) + 1;
            if ((i >= mm2_min) && (i <= mm2_max))
                *(mms+i) = *(mms+i) + 1;
        }
    }
    else
    {
        if (notch_mm_where >= 0) 
        {
            mmx1 = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mmx2 = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
#if 0
            if (qual < QUALCUTOFF ) mmcolor = pink; else mmcolor = orange;
#endif
            mmcolor = set_which_mm(whichmm1); 
#if 0 // USE_GD
            gdImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);
#endif
            ImageFilledRectangle(im,mmx1,y2,mmx2,y1,mmcolor);
    
            mm_min = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mm_max = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
        }

        if ( ((qual > 0)) && (qual < QUALCUTOFF ) && (notch_mm_where >= 0) )
        {
            mm_minLQ = (int)(pxwidth * (double)(offset+notch_mm_where)); 
            mm_maxLQ = (int)(pxwidth * (double)(offset+notch_mm_where+1)); 
        }
        for (i=x1 ;i<=x2 ; i++)
        {
            if ((wcode == 1) || (wcode == 3)) // single hit
                *(line+i) = *(line+i) + 1;
            *(line2+i) = *(line2+i) + 1;
            if ((i >= mm_min) && (i <= mm_max))
            {
                *(mms+i) = *(mms+i) + 1;
            }
        }
    }
}



int gotten_genes_index;
#define MAXGOTGENES 50
char gotgenes[MAXGOTGENES][50];
int firstgene_flag  = 1;

static void paint_gene_annot(char khr[], unsigned int s, unsigned int e) 
{
char m[512];
    size_t read_status = 0;
    unsigned int start = (long int)0;
    unsigned int end = (long int)0;
    unsigned int prev_start = (long int)0;
    unsigned int prev_end = (long int)0;
    long int spot;
    int did_gene_name_flag = 0;
    char prev[512];
    double local_pxwidth;
    int frcolor;
    int j = 0;
    int x1,x2,y1,y2;
    struct flattype *z;
    struct flattype f;
    unsigned int jj;

    if (firstgene_flag == 1)
    {
        memset(gotgenes,0,sizeof(gotgenes));
        firstgene_flag = 0;
        gotten_genes_index = 0;
    }
    local_pxwidth = (double)iw/(double)(e - s);
    x1 = x2 = 0;
    y1 = 1;
    y2 = 8;
    frcolor = purple;
    strcpy(prev,"");

    fix_up_support_file_paths();
sprintf(m,"in paint_gene_annot() blds = [%s], khr=%s",blds ,khr);  jdebug(m); 
sprintf(m,"still in paint_gene_annot() blds = [%s] refflat_d_n=[%s] refflat_o_fn=[%s]",blds ,refflat_d_fn,refflat_o_fn);  jdebug(m); 
    if (setup_refflat(refflat_d_fn,refflat_o_fn) < 0)              // 0 == good
    {
sprintf(m,"ERROR: after setup_refflat() - in paint_gene_annot ");  jdebug(m); 
        return;
    }
sprintf(m,"after setup_refflat() ");  jdebug(m); 

    memset(&f,0,sizeof(f));
    strcpy(f.chrom,khr);
    f.txStart = s;
    f.txEnd = s+100;

sprintf(m,"in paint_gene_annot() , before binary_search_refflat_file() [%s]:%d.%d",khr,s,s+100);  jdebug(m); 
    flat_search_count = 0;
    z = binary_search_refflat_file(&f,0,refflat_fixed_hi);
sprintf(m,"after binary_search_refflat_file() z=%p",z);  jdebug(m); 
// NO !!! if (!z) return; WE JUST WANT TO GET CLOSE !!!


// not really expecting a match, but need last looked at location in the binary refflat for bearings
    spot = last_refflat_fseek_place - (SIZE_REFFLAT_REC*50);
    if (spot < 0) spot = 0;

sprintf(m,"in paint_gene_annot() after binary_search_refflat_file khr=%s  spot=%ld ",khr,spot); jdebug(m);
    fseek(fp_refflat_d,last_refflat_fseek_place,SEEK_SET);
    j = 0;
    while (j<300)
    {
        // NO ! structs messed up sizes in different compilers . stat = fread(&temp_refflat_space,1,SIZE_REFFLAT_REC,fp_refflat_d);

        if (fread (&temp_refflat_space.geneName,MAX_GENENAME,1,fp_refflat_d) != 1) break;
        if (fread (&temp_refflat_space.name,MAX_GENENAME,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.chrom,MAX_GENENAME,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.txStart,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.txEnd,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.cdsStart,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.cdsEnd,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.exonCount,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.exonStarts,4,1,fp_refflat_d) != 1 ) break; // really offsets
        if (fread (&temp_refflat_space.exonEnds,4,1,fp_refflat_d) != 1 ) break;
        if (fread (&temp_refflat_space.strand,4,1,fp_refflat_d) != 1 ) break;

// sprintf(m,"j=%d wantchrloc=%s gotchr=%s gene?=%s gotloc=%d s=%d e=%d",j,khr,temp_refflat_space.chrom,temp_refflat_space.geneName,temp_refflat_space.txStart,s,e); jdebug(m);

        if (( (s >= temp_refflat_space.txStart ) && (s <= temp_refflat_space.txEnd )) ||   // start within 
            ( (e >= temp_refflat_space.txStart ) && (e <= temp_refflat_space.txEnd )) ||   // end within 
            ( (s <= temp_refflat_space.txStart ) && (e > temp_refflat_space.txEnd )) )      // spans
        {

/*
 sprintf(m,"in paint_gene_annot() here 3 j=%d s=%d e=%d txs=%d txe=%d cnt=%d",
    j,s,e,temp_refflat_space.txStart,temp_refflat_space.txEnd,temp_refflat_space.exonCount); jdebug(m);
*/

// sprintf(m,"in paint_gene_annot() -- fseek %u ",temp_refflat_space.exonStarts); jdebug(m);

            fseek(fp_refflat_o, (long int)(temp_refflat_space.exonStarts), SEEK_SET);    // really an offset into an "overflow" file
            for (jj=0 ; ((jj<temp_refflat_space.exonCount) && (jj < 2000)) ; jj++)  // kickout after 2000 to prevent lock up if read data is corrupted
            {
// sprintf(m,"in paint_gene_annot() -- looping %d ",jj);  jdebug(m);
                end = start = (unsigned int) 0;
                read_status = fread(&start,sizeof(unsigned int),1,fp_refflat_o); if (read_status != 1) break;
                read_status = fread(&end,sizeof(unsigned int),1,fp_refflat_o); if (read_status != 1) break;

                x1 = (int)(local_pxwidth * (double)(start -s)); 
                x2 = (int)(local_pxwidth * (double)(end -s)); 
                if (x2==x1) x2 = x1 + 1;
// sprintf(m,"in paint_gene_annot() -- looping %d , x1 = %d , x2 = %d ,iw = %d ",jj,x1,x2,iw);  jdebug(m);
                if ((x1<0) && (x2 < 0)) continue;
                if ((x1>=iw) && (x2 >= iw)) continue;
                if (x1 < 0) x1 = 0;
                if (x1>=iw) x1 = iw-1;;
                if (x2 < 0) x2 = 0;
                if (x2>=iw) x2 = iw-1;;

#if 0 // USE_GD
                gdImageFilledRectangle(im,x1,y1,x2,y2,black);
#endif
                ImageFilledRectangle(im,x1,y1,x2,y2,black);

                if (strcmp(prev,temp_refflat_space.geneName))
                {
// sprintf(m,"in paint_gene_annot() -- still here 4 [%s]",temp_refflat_space.geneName);  jdebug(m);
                    unsigned char *ucptr = (unsigned char *)&temp_refflat_space.geneName[0];
#if 0 // USE_GD
                    int brect[8];
#ifdef WIN32
                    gdImageString(im,gdFontGetLarge(),x1,22,ucptr,blue);
#else
                    gdImageStringFT (im, brect, blue, "" , 10, 0, x1,22 /*ih-30*/,(char *)ucptr);
#endif
#endif

                    ImageString(im,x1,22,ucptr,blue);

                      if (gotten_genes_index < MAXGOTGENES)
                      {
                          strcpy(gotgenes[gotten_genes_index++],temp_refflat_space.geneName);
                      }
                      strcpy(prev,temp_refflat_space.geneName);
                      did_gene_name_flag = 1;
                }
                if (jj)
                {
                    int x3 = (int)(local_pxwidth * (double)(prev_start-s));
                    if (x3 < 0) x3 = 0;
                    if (x3 >= iw) x3 = iw;
#if 0 // USE_GD
                    gdImageFilledRectangle(im,x3,y1+2,x2,y1+3,black);
#endif
                    ImageFilledRectangle(im,x3,y1+2,x2,y1+3,black);
                }
                prev_end = end;
                prev_start = start;

            }
        }
        j++;
    }

sprintf(m,"in paint_gene_annot() END "); jdebug(m);

    return;
}



static char *twobit(char khr[], int s, int e,char blds_arg[] ) // before sure to free return value
{
    int status;
    char *ptr;
    char fn[FILENAME_MAX];
char m[5120];

    status = 0;

sprintf(m,"in twobit() blds_arg=[%s] before setup_2bit_file",blds_arg); jdebug(m);
    setup_2bit_file(fn,blds_arg); 
    ptr = rd2bit(fn,khr,s,e,&status);
    if (status < 0) 
    {
sprintf(m,"ERROR: twobit() fn=%s %s %d %d ptr=%p, status=%d ",fn,khr,s,e,ptr,status); jdebug(m);
    }
sprintf(m,"after twobit() fn=%s %s %d %d ptr=%p, status=%d ",fn,khr,s,e,ptr,status); jdebug(m);
    return ptr;
}


   /* draw reference geneome region as colorful "rainbow" on bottom of image */
static void setup_reference_dna(char khr[], int s, int e,int paint_it_flag)
{
    double local_pxwidth;
    int yoffset = 18;
    int xcolor;
    int frcolor;
    int i;
    int x1,x2,y1,y2;
    int lastx1, lastx2;
    int slice_len_in_bases;
char m[512];

    if (dnaspace) free(dnaspace);
// use 2bit
    dnaspace = twobit(khr, s , e+1,blds); // free return if return != 0
    if (!dnaspace) 
    {
sprintf(m,"ERROR: in setup_reference_dna() dnaspace from twobit() is null "); jdebug(m);
        return;
    }
    if (paint_it_flag == 0) return;
 
    slice_len_in_bases = e - s;
    local_pxwidth = (double)iw/(double)(slice_len_in_bases);
    x1 = x2 = y1 = y2 = 0;
    frcolor = purple;

    lastx1 = lastx2 = 0;
    y1 = ih-3-yoffset;
    y2 = ih-7-yoffset;
sprintf(m,"in rainbow(), got dnaspace  %p len=%d y1=%d y2=%d ih=%d off=%d",dnaspace,slice_len_in_bases,y1,y2,ih,yoffset); jdebug(m);
    for (i=0 ; i<slice_len_in_bases ; i++)
    {
        x1 = (int)(local_pxwidth * (double)(i)); 
        x2 = (int)(local_pxwidth * (double)(i+1)); 
        if (x1==x2) x2++;
        if ((lastx1 == x1) && (lastx2 == x2)) continue;
        lastx2 = x2;
        if (x1>=iw) break;
        if (x2>=iw) break;
 
// sprintf(m,"dnaspace(%d) %c %d ",i,*(dnaspace+i),*(dnaspace+i)); jdebug(m);
        switch (*(dnaspace+i))
        {
            case 'A': xcolor = dna_a; break;
            case 'C': xcolor = dna_c; break;
            case 'G': xcolor = dna_g; break;
            case 'T': xcolor = dna_t; break;
            default : xcolor = pink;
        }
        ImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
        lastx1 = x1;
        lastx2 = x2;
    }
    return;
}


static void draw_dnacnts_and_dnamms(char khr[], int s, int e, int offset)
{
    double d;
    double local_pxwidth;
    int xcolor;
    int i;
    int x1,x2,y1,y2;
    int len;
// char m[512];

// sprintf(m,"in draw_dnacnts_and_dnamms() %s %d %d",khr,s,e); jdebug(m);
    if (!dnacnts) 
    {
        if (dnamms) free(dnamms);
        dnamms = (unsigned short int *)0;
    }
    if (!dnamms) return;
 
    len = e - s;
    local_pxwidth = (double)iw /(double)(len);
    x1 = x2 = y1 = y2 = 0;

    y1 = ih-3-offset;
// sprintf(m,"in draw_dnacnts_and_dnamms() - len = %d ,pxwidth=%10.3f",len,local_pxwidth); jdebug(m);

// THREE PASSES ... reds first, then blues ... then reds
    xcolor=red;
#define MMHEIGHT 22.0
#define MMCUTOFF 8
#define CHKMM if (((*(dnacnts+i) >= MMCUTOFF)  && (*(dnamms+i) >= 2)) || ((*(dnacnts+i) >= (MMCUTOFF+4))  && (*(dnamms+i) >= 1)))
    for (i=0 ; i<len ; i++)
    {
        CHKMM
        {
//sprintf(m,"in draw_dnacnts_and_dnamms() HIT - i=%d of %d mms=%d cnts=%d  %10.5f",
//                      i,len,*(dnamms+i),*(dnacnts+i),((double)*(dnamms+i) / (double)*(dnacnts+i)));
//jdebug(m);
            d = ((double)*(dnamms+i) / (double)*(dnacnts+i));
            if (d <= 0.6666) continue;
            y2 = y1 - (int)(d * MMHEIGHT);
            x1 = (int)(local_pxwidth * (double)(i)); 
            x2 = (int)(local_pxwidth * (double)(i+1)); 
            if (x1==x2) x2++;
            if (x1>=iw) break;
            if (x2>=iw) break;
#if 0 // USE_GD
            gdImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
#endif
            ImageFilledRectangle(im,x1,y2,x2,y1,xcolor);

        }
    }
    xcolor = blue;
    for (i=0 ; i<len ; i++)
    {
        CHKMM
        {
//sprintf(m,"in draw_dnacnts_and_dnamms() HIT - i=%d of %d mms=%d cnts=%d  %10.5f",
//                      i,len,*(dnamms+i),*(dnacnts+i),((double)*(dnamms+i) / (double)*(dnacnts+i)));
//jdebug(m);
            double d;
            d = ((double)*(dnamms+i) / (double)*(dnacnts+i));
            if ((d>0.6666) || (d < 0.3333)) continue;
            y2 = y1 - (int)(d * MMHEIGHT);
            x1 = (int)(local_pxwidth * (double)(i)); 
            x2 = (int)(local_pxwidth * (double)(i+1)); 
            if (x1==x2) x2++;
            if (x1>=iw) break;
            if (x2>=iw) break;
#if 0 // USE_GD
            gdImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
#endif
            ImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
        }
    }
    xcolor=red;
    for (i=0 ; i<len ; i++)
    {
        CHKMM
        {
//sprintf(m,"in draw_dnacnts_and_dnamms() HIT - i=%d of %d mms=%d cnts=%d  %10.5f",
//                      i,len,*(dnamms+i),*(dnacnts+i),((double)*(dnamms+i) / (double)*(dnacnts+i)));
//jdebug(m);
            d = ((double)*(dnamms+i) / (double)*(dnacnts+i));
            if (d>0.33333)  continue;
            y2 = y1 - (int)(d * MMHEIGHT);
            x1 = (int)(local_pxwidth * (double)(i)); 
            x2 = (int)(local_pxwidth * (double)(i+1)); 
            if (x1==x2) x2++;
            if (x1>=iw) break;
            if (x2>=iw) break;
#if 0 // USE_GD
            gdImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
#endif
            ImageFilledRectangle(im,x1,y2,x2,y1,xcolor);
        }
    }
// sprintf(m,"end draw_dnacnts_and_dnamms() i=%d of %d",i,len); jdebug(m);

// dnacnts and dnamms this uses lots of memory, so delete asap
    if (dnamms) free(dnamms);
    dnamms = (unsigned short int *)0;
    if (dnacnts) free(dnacnts);
    dnacnts = (unsigned short int *)0;

    return;
}



#if 0
static void underline(int diff, int offset, int len, int downamt)
{
    int kolor;
    int x1,x2,y1;

 
    y1 = ih - 5;
    x1 = (int)(pxwidth * (double)offset); 
    x2 = x1 + (int)(pxwidth * (double)len); 

    if (downamt == 0)  kolor = lightblue;
    else               kolor = yellow;
#if 0 // USE_GD
    gdImageFilledRectangle(im,x1,y1+downamt,x2,y1+downamt+1,kolor);
#endif
    ImageFilledRectangle(im,x1,y1+downamt,x2,y1+downamt+1,kolor);

    return;
}
#endif


#if 0
static int fit_this_line(int start, int stop, int * p_pos, int i, int space_arg)
{
        if (p_pos[2*i] == 0 && p_pos[2*i+1] == 0)
                return 1;
        if (start > (p_pos[2*i+1] + space_arg))
                return 1;
        if (stop < (p_pos[2*i] - space_arg))
                return 1;
        return 0;
}
#endif
 

#if 0
static int find_f_pos(int start, int stop, int * p_pos, int space_arg, int num)
{
        int i =0;
 
        for(i =0; i<num; ++i)
        {
                if (fit_this_line(start, stop, p_pos, i, space_arg))
                {
                        if (p_pos[2*i] == 0 && p_pos[2*i+1] == 0)
                        {
                                p_pos[2*i] = start;
                                p_pos[2*i+1] = stop;
                        }
                        else
                        {
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
                                p_pos[2*i] = min(start, p_pos[2*i]);
                                p_pos[2*i +1] = max(stop, p_pos[2*i+1]);
                        }
                        return i;
                }
        }
        return -1;

}
#endif


#define H 10

#ifdef WEB_SERVER
void print_some_javascript_for_popups(void)
{
    printf("<script type=\"text/javascript\"> \n");
    printf("var numcheckboxes;\n");

    printf("function popSubmit(form,cnt)\n");
    printf("{\n");
    printf("    var i;\n");
    printf("    var ischkd;\n");
    printf("    var fldchkd;\n");
    printf("    var dat;\n");
    printf("    var url;\n");
    printf("    var srch;\n");

    printf("    dat = \"\";\n");
    printf("    url = \"\";\n");
    printf("    srch = \"\";\n");
    printf("    url = \"http://%s/%s/%s?\";\n",URL,CGIDIRIFNECESSARY,progname);
    printf("    for (i=0;i<=numcheckboxes-1;i++)\n");
    printf("    {\n");
    printf("        fldchkd = \"form.k\" + i + \".checked\";\n");
    printf("        ischkd =eval(fldchkd);\n");
    printf("        if ( (ischkd) && (ischkd != \"undefined\") )\n");
    printf("        {\n");
    printf("            val = \"form.k\" + i + \".getAttribute(\" + \"\\\"VALUE\\\"\" + \")\" ;\n");
    printf("            val = (eval(val));\n");
    printf("            fldname  = \"\" + i; \n");
    printf("            dat = dat + \"&\" + fldname + \"=\" + 1;\n");
    printf("        }\n");
    printf("        else\n");
    printf("        {\n");
    printf("            val = \"form.k\" + i + \".getAttribute(\" + \"\\\"VALUE\\\"\" + \")\" ;\n");
    printf("            val = (eval(val));\n");
    printf("            fldname  = \"\" + i; \n");
    printf("            dat = dat + \"&\" + fldname + \"=\" + 0;\n");
    printf("        }\n");
    printf("    }\n");
    printf("    srch = srch + \"&\" + \"fl1=\" + \"%s\" \n",fl1);
    printf("    srch = srch + \"&\" + \"fl2=\" + \"%s\" \n",fl2);

//     if (strlen(typ)) printf("    srch = srch + \"&\" + \"typ=\" + \"%s\" \n",typ);
    printf("    srch = srch + \"&xyz=\" + Math.random(); \n");
printf("    alert(url+srch+dat); \n");
printf("    alert(\"hello\"); \n");
// printf("    alert(document.getElementById(\"form1\").TEXTAREA) \n");
printf("    alert(document.forms[0].elements[0]); \n");
    
    printf("    window.location.search = srch\n");
    printf("    window.location.href = url + srch + dat;\n");
    printf("    return false;\n");
    printf("}\n"); 
printf("    \n"); 


printf("function SubmitPICK(obj)\n"); 
printf("{   \n"); 
printf("    var u;   \n"); 
printf("    \n"); 
printf("    u = \"../cgi-bin/alview?pickq=\";   \n"); 
printf("    u = u + encodeURIComponent(obj.pickq.value);   \n"); 
printf("    window.location = u;   \n"); 
printf("}   \n"); 

printf("function SubmitX(obj)    \n"); 
printf("{   \n"); 
printf("    var u=\"\";   \n"); 
printf("    \n"); 
printf("    u = \"../cgi-bin/alview?filenum=%d&\";   \n",filez_id); 
printf("    u = u + \"position=\" + encodeURIComponent(obj.position.value);   \n"); 
printf("    u = u + \"&iw=\" + encodeURIComponent(obj.iw.value);   \n"); 
printf("    u = u + \"&ih=\" + encodeURIComponent(obj.ih.value);   \n"); 
printf("    if (obj.spliceonly.checked) { u = u + \"&spliceonly=on\"; } \n"); 
printf("    if (obj.qual.checked) { u = u + \"&qual=on\"; } \n"); 

// note    printf("        fldchkd = \"form.k\" + i + \".checked\";\n");
// printf("    u = u + \"&nso=\" + encodeURIComponent(obj.nso.value);   \n"); 

// printf("    if (obj.geneannot.checked) { u = u + \"&geneannot=on\"; } \n");  always on 
printf("    if (obj.basecolors.checked) { u = u + \"&basecolors=on\"; } \n"); 
#if 0
printf("    if (obj.altonly.checked) { u = u + \"&altonly=on\"; } \n"); 
printf("    if (obj.uniq.checked) { u = u + \"&uniq=on\"; } \n"); 
printf("    if (obj.mmspot.checked) { u = u + \"&mmspot=on\"; } \n"); 
printf("    if (obj.exp.checked) { u = u + \"&exp=on\"; } \n"); 

// set tds_val
printf("    if (obj.tds[0].checked) { u = u + \"&tds=1\"; } \n"); 
printf("    else if (obj.tds[1].checked) { u = u + \"&tds=2\"; } \n"); 
printf("    else if (obj.tds[2].checked) { u = u + \"&tds=3\"; } \n"); 
printf("    else if (obj.tds[3].checked) { u = u + \"&tds=4\"; } \n"); 
printf("    else if (obj.tds[4].checked) { u = u + \"&tds=5\"; } \n"); 
printf("    else if (obj.tds[5].checked) { u = u + \"&tds=6\"; } \n"); 
printf("    else if (obj.tds[6].checked) { u = u + \"&tds=7\"; } \n"); 
printf("    else if (obj.tds[7].checked) { u = u + \"&tds=8\"; } \n"); 
printf("    else if (obj.tds[8].checked) { u = u + \"&tds=9\"; } \n"); 
printf("    else if (obj.tds[9].checked) { u = u + \"&tds=10\"; } \n"); 
printf("    else if (obj.tds[10].checked) { u = u + \"&tds=11\"; } \n"); 
#endif

//  printf("alert(\"yyy in SubmitX end = \"+u);\n");
printf("    window.location = u;   \n"); 
printf("}   \n"); 
printf("    \n"); 

printf("    \n"); 

printf("function SubmitY(obj,p)    \n"); 
printf("{   \n"); 
printf("    var u;   \n"); 
printf("    \n"); 
// printf("alert(\"yyy in SubmitY \");\n");
printf("    u = \"../cgi-bin/alview?filenum=%d&\";   \n",filez_id); 
printf("    u = u + \"position=\" + encodeURIComponent(p);   \n"); 
printf("    u = u + \"&iw=\" + encodeURIComponent(obj.iw.value);   \n"); 
printf("    u = u + \"&ih=\" + encodeURIComponent(obj.ih.value);   \n"); 

#if 0
printf("    if (obj.altonly.checked) { u = u + \"&altonly=on\"; } \n"); 
#endif

printf("    if (obj.basecolors.checked) { u = u + \"&basecolors=on\"; } \n"); 
printf("    if (obj.uniq.checked) { u = u + \"&uniq=on\"; } \n"); 
printf("    if (obj.qual.checked) { u = u + \"&qual=on\"; } \n"); 

#if 0
printf("    if (obj.mmspot.checked) { u = u + \"&mmspot=on\"; } \n"); 
printf("    if (obj.exp.checked) { u = u + \"&exp=on\"; } \n"); 
printf("    if (obj.nso.checked) { u = u + \"&nso=on\"; } \n"); 
// set tds_val
printf("    if (obj.tds[0].checked) { u = u + \"&tds=1\"; } \n"); 
printf("    else if (obj.tds[1].checked) { u = u + \"&tds=2\"; } \n"); 
printf("    else if (obj.tds[2].checked) { u = u + \"&tds=3\"; } \n"); 
printf("    else if (obj.tds[3].checked) { u = u + \"&tds=4\"; } \n"); 
printf("    else if (obj.tds[4].checked) { u = u + \"&tds=5\"; } \n"); 
printf("    else if (obj.tds[5].checked) { u = u + \"&tds=6\"; } \n"); 
printf("    else if (obj.tds[6].checked) { u = u + \"&tds=7\"; } \n"); 
printf("    else if (obj.tds[7].checked) { u = u + \"&tds=8\"; } \n"); 
printf("    else if (obj.tds[8].checked) { u = u + \"&tds=9\"; } \n"); 
printf("    else if (obj.tds[9].checked) { u = u + \"&tds=10\"; } \n"); 
printf("    else if (obj.tds[10].checked) { u = u + \"&tds=11\"; } \n"); 
#endif
printf("    window.location = u;   \n"); 
printf("}   \n"); 
printf("    \n"); 


    printf("function popOn(form,cnt)\n");
    printf("{\n");
printf("    for(i=0;i<form.length;i++) \n");
printf("    {\n");
printf("        c = form.elements[i].checked = true;\n");
printf("    }\n");
    printf("}\n");

    printf("function popOff(form,cnt)\n");
    printf("{\n");
printf("    for(i=0;i<form.length;i++) \n");
printf("    {\n");
printf("        c = form.elements[i].checked = false;\n");
printf("    }\n");
    printf("}\n");


#if 0
    printf("function popOffX(form,cnt)\n");
    printf("{\n");
printf("    for(i=0;i<form.length;i++) \n");
printf("    {\n");
printf("        c = form[i].elements;\n");
printf("        for(j=0;j<c.length;j++) \n");
printf("        {\n");
printf("            if(c[j].getAttribute(\"type\") == \"checkbox\") \n");
printf("            {\n");
printf("                    c[j].setAttribute(\"checked\",\"true\");\n");
printf("            }\n");
printf("        }\n");
printf("    }\n");
    printf("}\n");
#endif

/*
    printf("var changecallshadowbox = document.getElementById(\"changecallshadow\"); \n");
    printf("changecallshadowbox.style.display = \"none\"; \n");
    printf("changecallshadowbox.style.position = \"absolute\"; \n");
*/
printf("function submitenter(this_form_arg,e) \n"); 
printf("{ \n"); 
printf("var keycode; \n"); 
printf("if (window.event) keycode = window.event.keyCode; \n"); 
printf("else if (e) keycode = e.which; \n"); 
printf("else return true; \n"); 
printf(" \n"); 
printf("if (keycode == 13) \n"); 
printf("   { \n"); 
// printf("   myfield.form.submit(); \n"); 
printf("SubmitX(this_form_arg); \n"); 
printf("   return false; \n"); 
printf("   } \n"); 
printf("else \n"); 
printf("   return true; \n"); 
printf("} \n"); 

printf("function picksubmitenter(this_form_arg,e) \n"); 
printf("{ \n"); 
printf("var keycode; \n"); 
printf("if (window.event) keycode = window.event.keyCode; \n"); 
printf("else if (e) keycode = e.which; \n"); 
printf("else return true; \n"); 
printf(" \n"); 
printf("if (keycode == 13) \n"); 
printf("   { \n"); 
// printf("   myfield.form.submit(); \n"); 
printf("SubmitX(this_form_arg); \n"); 
printf("   return false; \n"); 
printf("   } \n"); 
printf("else \n"); 
printf("   return true; \n"); 
printf("} \n"); 

    printf("</script> \n");
    printf("\n"); 

    return;
}
#endif


static void submitbutton(int flag,char *label, char chr[], int start, int end)
{
   if (flag) 
       printf("<span class=\"button\"><INPUT TYPE=BUTTON OnClick=\"SubmitX(this.form); return false;\" VALUE=\"%s\"></span>\n",label);
   else
       printf("<span class=\"button\"><INPUT TYPE=BUTTON OnClick=\"SubmitY(this.form,'%s:%d-%d'); return false;\" VALUE=\"%s\"></span>\n",chr,start,end,label);
}



#if 0
static int has_cov_file(char *s,char fnarg[])
{
    int j;
    int got;

    fnarg[0] = (char)0;
    for (j=0 ; covfilez[j] ; j++)
    {
        got = 0;
        if (strstr(s,covfilez[j]))
        {
            strcpy(fnarg,covfilez[j]);
            return 1;
         }
    }

    return 0;
}
#endif


#if WEB_SERVER
static int get_full_path_for_short_file_name(char sfn[]) 
{
    int i;

    for (i=0 ; num_filez ; i++) 
    {
       if (strstr(filez[i].fullpath,sfn) )
       {
           return filez[i].file_id;
       }
    }
    return -1;
}

static int get_tds_for_lfn(char lfn[]) 
{
    int i;
    for (i=0 ; i < num_filez ; i++) 
    {
       if (strcmp(filez[i].fullpath,lfn) == 0)
       {
           return filez[i].file_id;
       }
    }
    return -1;
}


int find_filez_id(char *frce)
{
    int i;

    for (i=0 ; i< num_filez; i++) 
    {
       if (strcmp(filez[i].fullpath,frce) == 0)
       {
           return filez[i].file_id;
       }
    }
    return -1;
}


int findtcga(char tcga[])
{
#define RIM(A,B) if (strcmp(A,tcga)== 0) return find_filez_id(B);

RIM("TCGA-09-0365.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-09-0365-10A-01W-0372-09_capture.bam" )  
RIM("TCGA-09-0365.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-09-0365-01A-02W-0372-09_capture.bam" )  
RIM("TCGA-09-0369.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-09-0369-10C-01W-0372-09_capture.bam" )  
RIM("TCGA-09-0369.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-09-0369-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-10-0931.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-10-0931-11A-01W-0420-08.bam" )  
RIM("TCGA-10-0931.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-10-0931-01A-01W-0420-08.bam" )  
RIM("TCGA-10-0934.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-10-0934-11A-01W-0420-08.bam" )  
RIM("TCGA-10-0934.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-10-0934-01A-02W-0420-08.bam" )  
RIM("TCGA-13-0724.N" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0724-10B-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0724.T" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0724-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0726.N" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0726-10B-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0726.T" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0726-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0755.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0755-10A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0755.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0755-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0760.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0760-10A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0760.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0760-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0768.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0768-10A-01W-0371-08.bam" )  
RIM("TCGA-13-0768.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0768-01A-01W-0371-08.bam" )  
RIM("TCGA-13-0795.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0795-10A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0795.T" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0795-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0800.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0800-10A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0800.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0800-01A-01W-0372-09_capture.bam" )  
RIM("TCGA-13-0803.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0803-10A-01W-0371-08.bam" )  
RIM("TCGA-13-0803.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0803-01A-01W-0371-08.bam" )  
RIM("TCGA-13-0890.N" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0890-10A-01D-0421-09_whole.bam" )  
RIM("TCGA-13-0890.T" , "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-13-0890-01A-01D-0421-09_whole.bam" )  
RIM("TCGA-13-0904.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0904-10A-01W-0420-08.bam" )  
RIM("TCGA-13-0904.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0904-01A-02W-0420-08.bam" )  
RIM("TCGA-13-0913.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0913-10A-01W-0420-08.bam" )  
RIM("TCGA-13-0913.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0913-01A-01W-0420-08.bam" )  
RIM("TCGA-23-1116.N" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-23-1116-10A-01W-0486-08.bam" )  
RIM("TCGA-23-1116.T" , "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-23-1116-01A-01W-0486-08.bam" )  
RIM("TCGA-10-0927.N", "/tcga_next_gen/TCGABCM/BAMFile/TCGA-10-0927-11A-01D-0445-10_SOLiD.bam" )  
RIM("TCGA-10-0927.T", "/tcga_next_gen/TCGABCM/BAMFile/TCGA-10-0927-01A-02D-0445-10_SOLiD.bam" )  
RIM("TCGA-13-0720.N", "/tcga_next_gen/TCGABCM/BAMFile/TCGA-13-0720-10B-01D-0445-10_SOLiD.bam" )  
RIM("TCGA-13-0720.T", "/tcga_next_gen/TCGABCM/BAMFile/TCGA-13-0720-01A-01D-0445-10_SOLiD.bam" )  
RIM("TCGA-13-0723.N", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0723-10B-01W-0372-09_whole.bam" )  
RIM("TCGA-13-0723.T", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0723-01A-02W-0372-09_whole.bam" )  
// RIM("TCGA-13-0751.N", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0751-10A-01W-0371-08.bam" )  
// RIM("TCGA-13-0751.T", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0751-01A-01W-0371-08.bam" )  
RIM("TCGA-13-0751.N", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0751-10A-01D-0446-08.bam" )  
RIM("TCGA-13-0751.T",  "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0751-01A-01D-0446-08.bam")  
RIM("TCGA-13-0890.N", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0890-10A-01W-0421-09_whole.bam" )  
RIM("TCGA-13-0890.T", "/tcga_next_gen/NG_buckets/bucket2/bam/TCGA-13-0890-01A-01W-0421-09_whole.bam" )  
RIM("TCGA-24-0980.N", "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-24-0980-10A-01D-0515-09_whole.bam" )  
RIM("TCGA-24-0980.T", "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-24-0980-01A-01D-0515-09_whole.bam" )  
RIM("TCGA-24-1103.N", "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-24-1103-10A-01D-0515-09_whole.bam" )  
RIM("TCGA-24-1103.T", "/tcga_next_gen/NG_buckets/bucket1/bam/TCGA-24-1103-01A-01D-0515-09_whole.bam" )  

return -1;
}

int findforcefile(char frce[])
{
    int i;
    for (i=0 ; i<num_filez ; i++) 
    {
       if (strcmp(filez[i].fullpath,frce) == 0)
           return filez[i].file_id;
    }
    return -1;
}

int find_fileid_for_bamname(char frce[])
{
    int i;
    char m[1024];

    for (i=0 ; i<num_filez ; i++) 
    {
       if (strstr(filez[i].fullpath,frce))
       {
sprintf(m,"FOUND at %d cmp %s, %s in find_fileznum_for_bamname<br>", i, filez[i].fullpath ,frce);  jdebug(m); 
           return filez[i].file_id; // return id, not index
       }
    }
    return -1;
}


void do_filez_form(void)
{
    int i;
    char m[512];

sprintf(m,"in do_filez_form() start , num_filez = %d",num_filez);  jdebug(m); 
    printf("<br>\n");
    for (i=0 ; i<num_filez ; i++) 
    {
#if 0
sprintf(m,"in do_filez_form() i=%d , num_filez = %d",i,num_filez);  jdebug(m); 
sprintf(m,"in do_filez_form() i=%d , num_filez = %d, name=[%s] ",i,num_filez,filez[i].fullpath);  jdebug(m); 
#endif
       if (filez[i].fullpath) 
           printf("<a href=\"../cgi-bin/alview?position=chr1:100-500&iw=1000&ih=350&filenum=%d\">%s</a>", filez[i].file_id,filez[i].fullpath);
       else
           printf("%s", filez[i].info); // comment or information line from file
/* char fn[512]; if (has_cov_file(filez[i].fullpath,fn) ) printf("&nbsp;&nbsp <small>*** %s</small>\n",fn); */
       printf("<br>\n"); 
    }
sprintf(m,"in do_filez_form() end , i = %d of %d ",i,num_filez);  jdebug(m); 

    return;

}
// if WEB_SERVER
#endif


void do_form(char chr[], int start, int end, int s, int e)
{
    int j;
    int tmp1;
    int flushtr;

    tmp1 =  0;
    flushtr = 0;

    printf("<br>\n");
    printf("<form name=\"form1\" action=\"../cgi-bin/%s\">\n",progname);
    printf("Position: ");
#if 0
    printf("<input type=\"text\" name=\"position\" id=\"position\" size=\"40\" value=\"%s:%d-%d\"> \n",chr,start,end);
#else
    // support carriage return 
    printf("<input type=\"text\" name=\"position\" id=\"position\" size=\"40\" value=\"%s:%d-%d\" onKeyPress=\"return submitenter(this.form,event)\"> \n",chr,start,end);

#endif
    if (start == end) 
    {
        start--;
        if (start < 1)  start = 1;
    }
    if (start >= end) 
    {
        end = start + 1;
    }
    printf("Width:<input type=\"text\" name=\"iw\" id=\"iw\" size=\"6\" value=\"%d\"> \n",iw);
    printf("Height:<input type=\"text\" name=\"ih\" id=\"ih\" size=\"6\" value=\"%d\"> \n",ih);

    if (spliceonly_flag == 1) 
        printf("&nbsp;&nbsp; Splice Only On:<input type=\"checkbox\" name=\"spliceonly\" checked>");
    else
        printf("&nbsp;&nbsp; Splice Only On:<input type=\"checkbox\" name=\"spliceonly\" >");

#if NONREFSPLICE
// old stuff
    if (altonly_flag == 1) 
        printf("&nbsp;&nbsp; Nonref-splice Only:<input type=\"checkbox\" name=\"altonly\" checked>");
    else
        printf("&nbsp;&nbsp; Nonref-splice Only:<input type=\"checkbox\" name=\"altonly\" >");
#endif

    if (basecolors_flag == 1) 
        printf("&nbsp;&nbsp; BaseColor:<input type=\"checkbox\" name=\"basecolors\" checked>");
    else
        printf("&nbsp;&nbsp; BaseColor:<input type=\"checkbox\" name=\"basecolors\" >");

    if (uniq_flag == 1) 
        printf("&nbsp;&nbsp; Uniq:<input type=\"checkbox\" name=\"uniq\" checked>");
    else
        printf("&nbsp;&nbsp; Uniq:<input type=\"checkbox\" name=\"uniq\" >");

    if (qual_flag == 1) 
        printf("&nbsp;&nbsp; Quality:<input type=\"checkbox\" name=\"qual\" checked>");
    else
        printf("&nbsp;&nbsp; Quality:<input type=\"checkbox\" name=\"qual\" >");

#if 0
    if (extrapolate_flag) 
        printf("&nbsp;&nbsp; Extrapolate:<input type=\"checkbox\" name=\"exp\" checked>");
    else
        printf("&nbsp;&nbsp; Extrapolate:<input type=\"checkbox\" name=\"exp\" >");

    if (nso_flag) 
        printf("&nbsp;&nbsp; Non Specfics On:<input type=\"checkbox\" name=\"nso\" checked>");
    else
        printf("&nbsp;&nbsp; Non Specfics On:<input type=\"checkbox\" name=\"nso\" >");

    if (mmspot_flag) 
        printf("&nbsp;&nbsp; MMmark:<input type=\"checkbox\" name=\"mmspot\" checked>");
    else
        printf("&nbsp;&nbsp; MMmark:<input type=\"checkbox\" name=\"mmspot\" >");
#endif

    printf("<br>");

    printf("&nbsp;\n");
//printf("<span class=\"submit\"><INPUT TYPE=BUTTON OnClick=\"SubmitX(this.form,'%s:%d-%d'); return false;\" VALUE=\"Submit\"></span>\n",chr,start,end);

    j = (end - start) / 2;

    submitbutton(1,(char *)"Submit",chr,start,end);
printf("<font size=-1>%u-%u , mismatch </font>\n",s,e); 
printf("<font color=\"%s\">A</font>&nbsp;\n",dna_a_s);  // dna_a
printf("<font color=\"%s\">C</font>&nbsp;\n",dna_c_s); 
printf("<font color=\"%s\">G</font>&nbsp;\n",dna_g_s); 
printf("<font color=\"%s\">T</font>&nbsp;\n",dna_t_s); 
printf("<font color=\"%s\">Ins</font>&nbsp;\n",dna_I_s); 
printf("<font color=\"%s\">Del</font>&nbsp;\n",dna_D_s); 
       printf( "<a href=\"../cgi-bin/alview?menu=1\">mainmenu</a>\n");
    if (opticaldupecount)
        printf("NOTE: %d optical/pcr dupes stripped",opticaldupecount); 
    if (badqualcount)
        printf("NOTE: %d bad quals stripped",badqualcount); 
    if (nocigarcount)
        printf("NOTE: %d bad cigars stripped",nocigarcount); 
    printf(" ins=%d del=%d tot=%ld blds=%s",icnt,dcnt,totalcnt,blds);
printf("<br>\n");
    submitbutton(0,(char *)"<Page",chr,start-(end-start),end-(end-start));
    submitbutton(0,(char *)"Page>",chr,start+(end-start),end+(end-start));
    submitbutton(0,(char *)"lefthalf",chr,start,start+j);
    submitbutton(0,(char *)"righthalf",chr,end-j,end);
    submitbutton(0,(char *)"ZoomIN",chr,start+(j/2),end-(j/2));

    tmp1 = (int)(start-((end-start)*1.5));
    if (tmp1 < 0) tmp1 = 1; 
     
    submitbutton(0,(char *)"ZoomOUT",chr,tmp1,(int)(end+((end-start)*1.5)));

    submitbutton(0,(char *)"base",chr,(start+((end-start)/2))-(iw/2),(start+((end-start)/2)-(iw/2))+iw);

printf("<br>\n");
    submitbutton(0,(char *)"<10",chr,start-10,end-10);
    submitbutton(0,(char *)"10>",chr,start+10,end+10);
    submitbutton(0,(char *)"<100",chr,start-100,end-100);
    submitbutton(0,(char *)"100>",chr,start+100,end+100);
    submitbutton(0,(char *)"<1000",chr,start-1000,end-1000);
    submitbutton(0,(char *)"1000>",chr,start+1000,end+1000);
    submitbutton(0,(char *)"<10,000",chr,start-10000,end-10000);
    submitbutton(0,(char *)"10,000>",chr,start+10000,end+10000);
printf("<br>\n");
    submitbutton(0,(char *)"<100,000",chr,start-100000,end-100000);
    submitbutton(0,(char *)"100,000>",chr,start+100000,end+100000);
    submitbutton(0,(char *)"<1,000,000",chr,start-1000000,end-1000000);
    submitbutton(0,(char *)"1,000,000>",chr,start+1000000,end+1000000);
    submitbutton(0,(char *)"<10,000,000",chr,start-10000000,end-10000000);
    submitbutton(0,(char *)"10,000,000>",chr,start+10000000,end+10000000);

/*
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">RightPage</a> | \n",
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Left100</a> | \n",
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Right100</a> | \n",chr,start+100,end+100);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Left500</a> | \n",chr,start-500,end-500);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Right500</a> | \n",chr,start+500,end+500);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Left1000</a> | \n",chr,start-1000,end-1000);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Right1000</a> | \n",chr,start+1000,end+1000);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Left10000</a> | \n",chr,start-10000,end-10000);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">Right10000</a> | \n",chr,start+10000,end+10000);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">ZoomOut 2X</a> | \n",chr,start-((end-start)*2),end+((end-start)*2));
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d\">ZoomIn 1/2</a> \n",chr,start+j,end-j);
printf("<br>\n");
*/

//printf("<span class=\"submit\"><INPUT type=\"submit\" name=\"Submit\" value=\"Submit\"></span>\n");
    printf("</form>\n");

/*
printf("<script type=\"text/javascript\"> \n");
printf("numcheckboxes=%d;\n",numspec);
printf("</script>\n");
*/
}


void print_some_css(void)
{
    printf("<style type=\"text/css\"> \n");
// style from http://css-tricks.com/examples/ButtonMaker/# 
printf(".button { \n");
printf("   border-top: 1px solid #96d1f8; \n"); 
printf("   background: #65a9d7; \n"); 
printf("   background: -webkit-gradient(linear, left top, left bottom, from(#3e779d), to(#65a9d7)); \n"); 
printf("   background: -webkit-linear-gradient(top, #3e779d, #65a9d7); \n"); 
printf("   background: -moz-linear-gradient(top, #3e779d, #65a9d7); \n"); 
printf("   background: -ms-linear-gradient(top, #3e779d, #65a9d7); \n"); 
printf("   background: -o-linear-gradient(top, #3e779d, #65a9d7); \n"); 
printf("   padding: 2px 4px; \n"); 
printf("   -webkit-border-radius: 8px; \n"); 
printf("   -moz-border-radius: 8px; \n"); 
printf("   border-radius: 8px; \n"); 
printf("   -webkit-box-shadow: rgba(0,0,0,1) 0 1px 0; \n"); 
printf("   -moz-box-shadow: rgba(0,0,0,1) 0 1px 0; \n"); 
printf("   box-shadow: rgba(0,0,0,1) 0 1px 0; \n"); 
printf("   text-shadow: rgba(0,0,0,.4) 0 1px 0; \n"); 
printf("   color: white; \n"); 
printf("   font-size: 13px; \n"); 
printf("   font-family: Helvetica, Arial, Sans-Serif; \n"); 
printf("   text-decoration: none; \n"); 
printf("   vertical-align: middle; \n"); 
printf("   } \n"); 
printf(".button:hover { \n"); 
printf("   border-top-color: #28597a; \n"); 
printf("   background: #28597a; \n"); 
printf("   color: #ccc; \n"); 
printf("   } \n"); 
printf(".button:active { \n"); 
printf("   border-top-color: #9fceed; \n"); 
printf("   background: #9fceed; \n"); 
printf("   } \n"); 
printf("    </style> \n");
}


void get_out_of_this_mess(char msg[])
{
        printf("<head>\n"); 
        printf("</head>\n"); 
        printf("<body bgcolor=\"#%s\">\n",FAVEBACKGROUNDCOLOR); 
        printf("%s <br>\n",msg);
        printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview\">For ALVIEW Main Page, click here.</a></font>&nbsp;&nbsp;\n"); 

#if PICK
        printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?pickmenu=1\">For ALVIEW Pick Page, click here.</a></font><br>\n"); 
#endif

        zap_tmpalview_files();
        printf("</body>\n"); 
        printf("</html>\n"); 
        return;
}


#define MAXKEEPS 30000

struct keep_these_type
{
   char *name;
   char *hugo;
};

struct keep_these_type keep_these[MAXKEEPS];
int numkeeps=0;


void comma_gin(char *s)
{
    char tmps[MAXBUFF+10];
    char *z;


// sprintf(m,"in comma gin 1 [%s]",s);  jdebug(m); 

    while (1)
    {
        z = strstr(s,",");
        if (!z) break;
        strcpy(tmps,z+1);
        strcpy(z,tmps);
    }
}


#if 0
static int chr_sizes_hg18[] =
{
/* chr1 */        
    247249704 , 242951134 , 199501812 , 191273048 , 180857851 , 170899977 , 158821409 , 146274811 ,
    140273237 , 135374722 , 134452369 , 132349519 , 114142965 , 106368570 , 100338900 , 88827239 , 78774727 ,
    76117138 , 63811636 , 62435949 , 46944308 , 49691417 , 154913739 , 57772939 , 16556 , -1
};
#endif


int parse_position(const char argposition[],char chr[],int *start,int *end)
{ 
    int build;
    int i = 0;
    char tmps[1024];
    char t[MAXBUFF];
    char tmps1[MAXBUFF];
    char tmps2[MAXBUFF];
char m[1024];

    if (strcmp(blds,"hg18") == 0) build = 18;
    else if (strcmp(blds,"hg19") == 0) build = 19;
    else 
    {
        sprintf(m,"ERROR: invalid build in parse_position() %s",blds);
        jdebug(m);
        return -1;
    }

// wth ? if ((i>0)&&(i<=24))

sprintf(m,"start parse_position()  %s, bld=%d",argposition,build);  jdebug(m); 
    tmps[0] = tmps1[0] = tmps2[0] = t[0] = (char)0;
    strcpy(t,argposition);
    un_escape((unsigned char *)t);

    strcpy(tmps,t);
    for (i=0 ; tmps[i] ; i++)
    {
        if (tmps[i] == ':') { tmps[i] = ' '; }
        if (tmps[i] == '-') { tmps[i] = ' '; }
    }

    sscanf(tmps,"%s %s %s",chr,tmps1,tmps2);
sprintf(m,"in parse_position before comma_gin %s ",tmps1);  jdebug(m); 
    comma_gin(tmps1);
sprintf(m,"in parse_position after  comma_gin %s ",tmps1);  jdebug(m); 

    *start = atoi(tmps1); 
    comma_gin(tmps2);
    *end = atoi(tmps2); 

    i = atoi(&tmps[3]);
    i--;
    if (strcmp("chrX",chr) == 0) i = 22;
    else if (strcmp("chrY",chr) == 0) i = 23;
    else if (strcmp("chrM",chr) == 0) i = 24;
sprintf(m,"in parse_position i= %d , start=%d end=%d",i,*start,*end);  jdebug(m); 
    if ((i>0)&&(i<=24))
    {
        if (build == 18) 
        {
            if (*end>chrsizes_hg18[i]) *end = chrsizes_hg18[i];
            if (*start>chrsizes_hg18[i]-2) *start = chrsizes_hg18[i]-2;
        }
        else if (build == 19) 
        {
sprintf(m,"in parse_position i= %d , start=%d end=%d bld=%d",i,*start,*end,build);  jdebug(m); 
            if (*end>chrsizes_hg19[i]) *end = chrsizes_hg19[i];
            if (*start>chrsizes_hg19[i]-2) *start = chrsizes_hg19[i]-2;
sprintf(m,"in parse_position FIXED? i= %d , start=%d end=%d bld=%d",i,*start,*end,build);  jdebug(m); 
        }
    }
sprintf(m,"end parse_position %s %d %d",chr,*start,*end);  jdebug(m); 
    return 0;
}


int atypical_go_front_end(void)
{
    int i = 0;
    long contentlength;
    char *endptr = (char *)0;
    char *ptr2contentlength = (char *)0;
    char *ptr2querystring = (char *)0;
    char dbmsg[MAXBUFF];


           /* make some assumtions until proven otherwise ... */
    WINDOWS = 1;   /* where do you want to go today ??? */
    MSIE = 1;   /* microsoft internet explorer */
    MOZILLA = 0;  /* netscape or mozilla : bow to the beast !!! */

    ptr2querystring = (char *)0;
    ptr2contentlength = getenv("CONTENT_LENGTH");
    if (ptr2contentlength == (char *)0)
    {
        ptr2querystring  = getenv("QUERY_STRING");
        global_QUERY_STRING_ptr = ptr2querystring;
        if (ptr2querystring == (char *)0)
        {
             jdebug("no CONTENT_LENGTH in go_front_end");

             printf("Content-type: text/html\n\n");
             did_content_type = 1;
             printf("<html><body><H1>Invalid Input: No \"Content Length\" and No Query String for FE program<br>\n");
//             printf("Click <a href=\"%s/goev.html\"> here </a> to start or restart</H1></body></html>\n", FEURL);

             debug_cgi_environment_variables_or_whatever_they_are();

             fflush(stdout);
             return -1;
        }
    }

    if (ptr2contentlength)
    {
        contentlength = strtol(ptr2contentlength, &endptr, 10);
        if ( contentlength > MAXBUFF - 10)
        {
            bark();
        }
        fread(buff, contentlength, 1, stdin);

#if 1 /* debug stuff */
printf("Content-type: text/html\n\n");
did_content_type = 1;
printf("<html><body><br>\n");
printf("hello.<br>\n");
debug_cgi_environment_variables_or_whatever_they_are();
printf("<br>");
printf("buff len=%ld content = [%s]<br>\n",contentlength, buff);
printf("</body></html>\n");
return -1;
#endif
    }
else {
  // sprintf(dbmsg,"ptr2contentlength = NULL "); jdebug(dbmsg); 
     }

    if (ptr2querystring)
    {
        memset(buff,0,sizeof(buff));
        i = strlen(ptr2querystring);
        if (i > (int)(sizeof(buff) - 2) )
            i = (int)(sizeof(buff) - 2); 

        strncpy(buff,ptr2querystring,i);

    }
    else { sprintf(dbmsg,"ptr2querystring = NULL\n"); jdebug(dbmsg); }

    get_some_environment_variables(); /* will also set some global variables */

//     jdebug(dbmsg);

    return 0;
}


static int globalreadscount = 0;

static int g_min_mapQ = 0, g_flag_on = 0, g_flag_off = 0;
static char *g_library, *g_rg;

static inline int __g_skip_aln(const bam_header_t *h, const bam1_t *b)
{
        if (b->core.qual < g_min_mapQ || ((b->core.flag & g_flag_on) != g_flag_on) || (b->core.flag & g_flag_off))
                return 1;
        if (g_rg) {
                uint8_t *s = bam_aux_get(b, "RG");
                if (s && strcmp(g_rg, (char*)(s + 1)) == 0) return 0;
        }
        if (g_library) {
                const char *p = bam_get_library((bam_header_t*)h, b);
                return (p && strcmp(p, g_library) == 0)? 0 : 1;
        }
        return 0;
}

#define MAXSEQ 10240
char seq[MAXSEQ];

static int check_uniq_tag(const bam_header_t *header, const bam1_t *b) // search SAM tags
{
    uint8_t *s;
    s = bam1_aux(b);
    char *z;

#if 0
        if (strstr((char *)s,"XT:A:U")) return 1;
        return 0;
#endif
        if (!s) return 0;

        z = strstr((char *)s,"XTAU" ); 
        if (!z) return 0;
        return 1;
        while (s < b->data + b->data_len) 
        {
                uint8_t type, key[2];
                key[0] = s[0]; key[1] = s[1];
                s += 2; type = *s; ++s;
#if 0
                ksprintf(&str, "\t%c%c:", key[0], key[1]);
if (type == 'A') { ksprintf(&str, "A:%c", *s); ++s; }
                else if (type == 'C') { ksprintf(&str, "i:%u", *s); ++s; }
                else if (type == 'c') { ksprintf(&str, "i:%d", *s); ++s; }
                else if (type == 'S') { ksprintf(&str, "i:%u", *(uint16_t*)s); s += 2; }
                else if (type == 's') { ksprintf(&str, "i:%d", *(int16_t*)s); s += 2; }
                else if (type == 'I') { ksprintf(&str, "i:%u", *(uint32_t*)s); s += 4; }
                else if (type == 'i') { ksprintf(&str, "i:%d", *(int32_t*)s); s += 4; }
                else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
                else if (type == 'd') { ksprintf(&str, "d:%lg", *(double*)s); s += 8; }
                else if (type == 'Z' || type == 'H') { ksprintf(&str, "%c:", type); while (*s) kputc(*s++, &str); ++s; }
#endif
                if ((key[0] == 'X') && (key[1] == 'T') && (type=='A'))
                {
                   if (strstr((char *)s,"XT:A:U")) return 1;
                }
        }
     return 0;
}




#if 0
// old - out custom sam tag
static int check_XR_tag(const bam_header_t *header, const bam1_t *b)
{
 char m[512];
	// kstring_t str;
	uint8_t *s;

        s = bam_aux_get(b, "XR");
        if ((s = bam_aux_get(b, "XR")) != (void *)0)
        {
           if (*s=='A') 
           {
               if (*(s+1) == 'G') return 0;
               if (*(s+1) == 'N') return 1;
               if (*(s+1) == 'A') return 2;
               if (*(s+1) == 'E') return 3;
           }
        }
else
{
jdebug("faile to get XR tag");
}
        return -1;


#if 0
        s = bam1_aux(b);
        if (s == (void *)0) 
            return -1;
        while (s < b->data + b->data_len) 
        {
                uint8_t type, key[2];
                key[0] = s[0]; key[1] = s[1];
                s += 2; type = *s; ++s;
                if ((key[0] == 'X') && (key[1] == 'R')) 
                {
                    if (*(char *)s == 'G') return 0;
                    if (*(char *)s == 'N') return 1;
                    if (*(char *)s == 'A') return 2;
                    if (*(char *)s == 'E') return 3;
                }
        }
      return 0;
#endif
}
#endif


// old? static int kknt = 0;
const char *rpf_bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"; // a local copy, samtools library on windows is not exporting data for some reason

// rpf_format2() is called often
// rpf_format2 is called from view_func(), rpf_format2 call aldetails() 
static char *rpf_format2(const bam_header_t *header, const bam1_t *b)
{
    char cigar[512];
    int8_t *qbuff = (int8_t *)0;
    int splice_source;
    int i = 0;
    const bam1_core_t *c = &b->core;

#if 0
char m[1024];
        memset(&mystr,0,sizeof(mystr));
	mystr.l = mystr.m = 0; mystr.s = 0;

	ksprintf(&mystr, "%s\t%d\t", bam1_qname(b), c->flag);
	if (c->tid < 0) kputs("*\t", &mystr);
	else ksprintf(&mystr, "%s\t", header->target_name[c->tid]);
	ksprintf(&mystr, "%d\t%d\t", c->pos + 1, c->qual);
	if (c->n_cigar == 0) kputc('*', &mystr);
	else{
		for (i = 0; i < c->n_cigar; ++i)
			ksprintf(&mystr, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	} 
	kputc('\t', &mystr);
	if (c->mtid < 0) kputs("*\t", &mystr);
	else if (c->mtid == c->tid) kputs("=\t", &mystr);
	else ksprintf(&mystr, "%s\t", header->target_name[c->mtid]);

	ksprintf(&mystr, "%d\t%d\t", c->mpos + 1, c->isize);
	for (i = 0; i < c->l_qseq; ++i) kputc(rpf_bam_nt16_rev_table[bam1_seqi(s, i)], &mystr);
	kputc('\t', &mystr);
	if (t[0] == 0xff) kputc('*', &mystr);
	else for (i = 0; i < c->l_qseq; ++i) kputc(t[i] + 33, &mystr);
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
		ksprintf(&mystr, "\t%c%c:", key[0], key[1]);
		if (type == 'A') { ksprintf(&mystr, "A:%c", *s); ++s; }
		else if (type == 'C') { ksprintf(&mystr, "i:%u", *s); ++s; }
		else if (type == 'c') { ksprintf(&mystr, "i:%d", *s); ++s; }
		else if (type == 'S') { ksprintf(&mystr, "i:%u", *(uint16_t*)s); s += 2; }
		else if (type == 's') { ksprintf(&mystr, "i:%d", *(int16_t*)s); s += 2; }
		else if (type == 'I') { ksprintf(&mystr, "i:%u", *(uint32_t*)s); s += 4; }
		else if (type == 'i') { ksprintf(&mystr, "i:%d", *(int32_t*)s); s += 4; }
		else if (type == 'f') { ksprintf(&mystr, "f:%g", *(float*)s); s += 4; }
		else if (type == 'd') { ksprintf(&mystr, "d:%lg", *(double*)s); s += 8; }
		else if (type == 'Z' || type == 'H') { ksprintf(&mystr, "%c:", type); while (*s) kputc(*s++, &mystr); ++s; }
	}
// example @SOLEXA9_26:5:75:1325:688 chr10 123384948 50
/*
        i = c->pos+1 - start;
        if (i<0) i = 0;
        for (j=0;  (j<(c->l_qseq))  ; j++,i++)
        {  
            if (i>=maxspace) break;
             *(space+i) += 1;
        }
*/

#define BAM_OFDEC          0
	if (c->n_cigar == 0) 
        {
        }
	else {
                int i;
		for (i = 0; i < c->n_cigar; ++i)
		{
			sprintf(m, "cig:%s:%d%c", header->target_name[c->tid],
                               bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
		}
	}


// unused code ^^^^ 
#endif


#define bam1_opticaldupe(b) (((b)->core.flag&BAM_FDUP) != 0) 
#define bam1_badqual(b) (((b)->core.flag&BAM_FQCFAIL) != 0) 

       // #define BAM_FQCFAIL      512
        if (bam1_opticaldupe(b) == 1)
            opticaldupecount++;
        else if (bam1_badqual(b) == 1) badqualcount++;
        else 
        {
            if (c->n_cigar == 0)  // number of cigar operations
            {
               strcpy(cigar,"*");
               nocigarcount++;
            }
            else 
            {
                cigar[0] = (char)0;
                for (i = 0; i < c->n_cigar; i++)
                {
                    char tmps[512];
                    sprintf(tmps,"%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
                    strcat(cigar,tmps);
                }
                if (c->l_qseq >= 2048) // l_qseq = length of query sequence read
/// WHY IS THIS HERE ????? 
                {
                    paint_align_box(gdiff,c->pos+1-gs ,c->l_qseq,/*wcode*/0,/*notch*/-1,/*qual*/-1,/*fwd*/bam1_strand(b),/*bow*/0,0,0,0,
                             bam1_opticaldupe(b));
                }
                else
                {
 	            uint8_t *s = bam1_seq(b); // , *t = bam1_qual(b);
	            for (i = 0; i < c->l_qseq; i++) 
	            {
                        seq[i] = rpf_bam_nt16_rev_table[bam1_seqi(s, i)];
    	            }
                    seq[i] = (char)0;
#if 0 
// debug 
for (i = 0; i < c->n_cigar; ++i)
                  ksprintf(&str, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);

#endif

                  splice_source = 0;
                  if (uniq_flag)
                  {
                     if (check_uniq_tag(header,b) == 0) 
                         return (char *)0;
                  }
#if 0
                  if ((spliceonly_flag == 1) ||( altonly_flag == 1))
                  {
                     splice_source = check_XR_tag(header,b); // our custom tag
                  }
                  if (altonly_flag == 1) 
                  {
                        if (splice_source < 2) return (char *)0;
                  }
#endif
                  totalcnt++;
// do nucleotide base quality for the read 
                if (qual_flag == 1) 
                {
                    uint8_t *seq_q;
                    int qlen = b->core.l_qseq;

                    qbuff = (int8_t *)malloc(qlen+1);
                    memset(qbuff,0,qlen+1);
                    seq_q = bam1_qual(b);
                    for (i = 0; i < qlen; ++i)
                            qbuff[i] = 33 + seq_q[i];
                    if (b->core.flag & 16) { // reverse
                            for (i = 0; i < qlen>>1; ++i) {
                                    int8_t t = qbuff[qlen - 1 - i];
                                    qbuff[qlen - 1 - i] = qbuff[i];
                                    qbuff[i] = t;
                            }
                    }
                }

#ifdef CMD_LINE
                if (snp_call_flag == 1) 
                    snp_call_aldetails(gdiff,c->pos+1-gs ,c->l_qseq,
                                /*fwd*/bam1_strand(b), bam1_opticaldupe(b),seq,(char *)qbuff,c->l_qseq,cigar,splice_source);
                else
                    aldetails(gdiff,c->pos+1-gs ,c->l_qseq,
                                /*fwd*/bam1_strand(b), bam1_opticaldupe(b),seq,(char *)qbuff,c->l_qseq,cigar,splice_source);
#else
                aldetails(gdiff,c->pos+1-gs ,c->l_qseq,
                                /*fwd*/bam1_strand(b), bam1_opticaldupe(b),seq,(char *)qbuff,c->l_qseq,cigar,splice_source);
#endif
                }
            }
        }
        globalreadscount++;
	return (char *)0; // mystr.s;
}


// p_haspair,P_mappedproperpair, u_unmapped, U_mateunmapped, r_strand, R_matestrand, 1_pair1st, 2_pair2nd,s_notprime,
// f_failvendor,d_pcrdupe
// extern char *bam_flag2char_table;  windows samtools library not exporting data 
const char *rpf_bam_flag2char_table = "pPuUrR12sfd\0\0\0\0\0";


char *rpf_format_sam(const bam_header_t *header, const bam1_t *b, int of)
{
	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	int i;
	const bam1_core_t *c = &b->core;
	kstring_t str;
	str.l = str.m = 0; str.s = 0;

	ksprintf(&str, "%s\t", bam1_qname(b));

	if (of == BAM_OFDEC) ksprintf(&str, "%d\t", c->flag);
	else if (of == BAM_OFHEX) ksprintf(&str, "0x%x\t", c->flag);
	else { // BAM_OFSTR
		for (i = 0; i < 16; ++i)
			if ((c->flag & 1<<i) && rpf_bam_flag2char_table[i])
				kputc(rpf_bam_flag2char_table[i], &str);
		kputc('\t', &str);
	}
	if (c->tid < 0) kputs("*\t", &str);
	else ksprintf(&str, "%s\t", header->target_name[c->tid]);
	ksprintf(&str, "%d\t%d\t", c->pos + 1, c->qual);
	if (c->n_cigar == 0) kputc('*', &str);
	else {
		for (i = 0; i < c->n_cigar; ++i)
			ksprintf(&str, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	}
	kputc('\t', &str);
	if (c->mtid < 0) kputs("*\t", &str);
	else if (c->mtid == c->tid) kputs("=\t", &str);
	else ksprintf(&str, "%s\t", header->target_name[c->mtid]);
	ksprintf(&str, "%d\t%d\t", c->mpos + 1, c->isize);
	if (c->l_qseq) {
		for (i = 0; i < c->l_qseq; ++i) kputc(rpf_bam_nt16_rev_table[bam1_seqi(s, i)], &str);
		kputc('\t', &str);
		if (t[0] == 0xff) kputc('*', &str);
		else for (i = 0; i < c->l_qseq; ++i) kputc(t[i] + 33, &str);
	} else ksprintf(&str, "*\t*");
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
		ksprintf(&str, "\t%c%c:", key[0], key[1]);
		if (type == 'A') { ksprintf(&str, "A:%c", *s); ++s; }
		else if (type == 'C') { ksprintf(&str, "i:%u", *s); ++s; }
		else if (type == 'c') { ksprintf(&str, "i:%d", *s); ++s; }
		else if (type == 'S') { ksprintf(&str, "i:%u", *(uint16_t*)s); s += 2; }
		else if (type == 's') { ksprintf(&str, "i:%d", *(int16_t*)s); s += 2; }
		else if (type == 'I') { ksprintf(&str, "i:%u", *(uint32_t*)s); s += 4; }
		else if (type == 'i') { ksprintf(&str, "i:%d", *(int32_t*)s); s += 4; }
		else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'd') { ksprintf(&str, "d:%lg", *(double*)s); s += 8; }
		else if (type == 'Z' || type == 'H') { ksprintf(&str, "%c:", type); while (*s) kputc(*s++, &str); ++s; }
	}
printf("%s\n",str.s); 
	return str.s;
}  

int maxcut = 5;

void fetch(char fn[],long int pos,char whatthere[],int len)
{
    register int i;
    register int k;
    int headlen;
    long spot;
    FILE *fp;
    char s[1024];

// printf("fetch : %s %ld %d\n",fn,pos,len);
    whatthere[0] = (char)0;
    fp = fopen(fn,"rb");
    if (fp == (void *)0) { fprintf(stderr,"%s not open error.  Invalid filename. \n",fn); exit(0); }
    fgets(s,1022,fp);  // eat head line
    headlen = strlen(s) ;
    spot = headlen +  pos + (pos/50);
// printf("spot = %ld headlen = %d pos = %ld  len=%d\n",spot,headlen,pos,len);
    fseek(fp,spot,SEEK_SET);
    i  = 0;
    while (i < len)
    {
        k = fgetc(fp);
        if ((char)k == '\r') continue;
        if ((char)k == '\n') continue;
        whatthere[i++] = (char)k;
    }
    whatthere[i+1] = (char)0;
// printf("whatthere %d %s\n",i,whatthere);
    fclose(fp);
    return;
}

#if 0
#define MAX_DNA 1000002
static char whatthere[MAX_DNA ];

static int fetchwrap2(char *chr,int s, int e)
{
    int k;
    char fn[512];
    if ((e-s)>MAX_DNA-2)
    {
        printf("ERROR- too big %s %d %d , cant be bigger than %d \n",chr,s,e,MAX_DNA);
        return -1;
    }
    sprintf(fn,"/h1/finneyr/amd64/hg18/%s.fa",chr);
    whatthere[0] = (char)0;
    fetch(fn,s,whatthere,e-s);
/*
    for (k=0;whatthere[k] ;k++)
    {
        printf("%c",whatthere[k]);
        if ((k%50) == 49) printf("\n");
    }
    printf("\n");
*/
    k = strlen(whatthere); 
    return k;
}
#endif

char compli(char ch)
{
    if (ch == 'A') return 'T';
    else if (ch == 'T') return 'A';
    else if (ch == 'C') return 'G';
    else if (ch == 'G') return 'C';
    else if (ch == 'N') return 'N';
    else if (ch == 'a') return 't';
    else if (ch == 't') return 'a';
    else if (ch == 'c') return 'g';
    else if (ch == 'g') return 'c';
/*
    else if (ch == 'a') return 't';
    else if (ch == 't') return 'a';
    else if (ch == 'c') return 'g';
    else if (ch == 'g') return 'c';
*/
    return '.';
}

void reverse_and_compliment(char s[], char outs[], int len)
{
    register int i,j;

    j = 0;
    i = len - 1;
    while (i >= 0) 
        outs[j++] = compli(s[i--]);
    return;
}

#if GOTOH
// alignment code swiped from internet, non-profit use is OK ....
// *****   KEEP GOTOH set to 0 if you are a commercial outfit. *******   
// It should be easy to re-engineer if you realy need it.
// this is somewhat hacked by Richard Finney for alview

/*-------------------------------------------------------------

Copyright (C) 2000 Peter Clote. 
All Rights Reserved.

Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.


THE AUTHOR AND PUBLISHER MAKE NO REPRESENTATIONS OR
WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
THIS SOFTWARE OR ITS DERIVATIVES.

-------------------------------------------------------------*/


/*************************************************

	Program: gotoh.c

	Peter Clote, 24 April 1998

Program for sequence alignment, assuming that sequences
s,t both begin and end with the same character, which can
artificially be ensured by giving dummy start and stop characters
(such as #, but not $ or * because the shell interprets these).
This is the final, corrected version, which rectifies the
initialization errors in the presentation of Nei and Waterman by using
P(0,j) = Q(i,0) = INF.

The traceback path is determined by choosing that direction, for
which D(i-1,j-1), D(i-1,j) and D(i,j-1) is minimum, i.e.
for which 

                    _
                   |
                   | D(i-1,j-1) + GAMMA(i,j)  (diagonal)
D(i,j) =  MIN of   | P(i,j)                   (up)
                   | Q(i,j)                   (left)
                   |_

is a minimum.

*************************************************/


// #include <stdio.h>
// #include <string.h>
// #include <stdlib.h>

#define SPACE ' '
#define MIN2(x,y)   ( (x) <= (y) ? (x) : (y) )
#define MAX2(x,y)   ( (x) <= (y) ? (y) : (x) )
#define MIN(x,y,z)  MIN2(MIN2(x,y),z)
#define MISMATCH  3     /* cost of mismatched letters */
#define D(i,j)      d[(i)*(m+1) + (j)]
#define P(i,j)      p[(i)*(m+1) + (j)]
#define Q(i,j)      q[(i)*(m+1) + (j)]
#define INF 100

#define U  2
#define V  2
#define W(i)    U*(i)+V       /* cost of  gap of length i, i>0 */
#define S(i)    s[(i)-1]
#define T(i)    t[(i)-1]
	/* D is integer matrix with n+1 rows and m+1 cols,
	   where n = length of s, m = length of t */
#define GAMMA(i,j) ( S(i) == T(j) ? 0 : MISMATCH )


/*----------------------------------------------------------
Triplet is a triple (d,p,q) of pointers to the arrays
D,P,Q used in Gotoh's algorithm.
-----------------------------------------------------------*/

struct pointer_triplet
{
  int *first;
  int *second;
  int *third;
};

typedef struct pointer_triplet Triple;




void error ( char *s );     /* prototype of error function */

void error ( char *s )
{
	fprintf(stderr,"%s\n",s);
	exit(1);
}


Triple dist ( char *s, char *t )
{
	int i,j,k,n = strlen(s), m = strlen(t);
	Triple val;
	int *d, *p, *q;

	if ( s[0] != t[0] || s[n-1] != t[m-1] )
		error((char *)"Strings don't have same beginning or ending character.");

	/******* allocate arrays *******************/
	d = (int *) calloc((n+1)*(m+1), sizeof(int));
	p = (int *) calloc((n+1)*(m+1), sizeof(int));
	q = (int *) calloc((n+1)*(m+1), sizeof(int));

	/******** initialize border values of matrices */
	D(0,0) = 0;
        P(0,0) = Q(0,0) = 0;
	for (i=1;i<=n;i++) D(i,0) = W(i);
	for (j=1;j<=m;j++) D(0,j) = W(j);
	for (j=1;j<=m;j++)  P(0,j) = INF;
	for (i=1;i<=n;i++)  Q(i,0) = INF;

	/******** compute interior values using recurrence relation.
	 This loop is really of the form
				for k = 2 to n+m
					for i,j satisfying i+j = k
						.....
	*************************************************************/


	for (k=2;k<=n+m;k++)
		for (i=1;i<=MIN2(n,k-1);i++)
			{
			j = k-i;
			if ( j > m ) continue;     /* index out of bounds */
			P(i,j) = MIN2( D(i-1,j) + W(1), P(i-1,j) + U );
			Q(i,j) = MIN2( D(i,j-1) + W(1), Q(i,j-1) + U );
			D(i,j) = MIN( D(i-1,j-1) + GAMMA(i,j), P(i,j), Q(i,j) );
			}
	val.first = d;
	val.second = p;
	val.third = q;
	return(val);
}



char *backtrack ( Triple val, char *s, char *t, int n, int m )
{
	/* compute the directions, starting from the LAST
	   character of sequence s,t in order to determine
	   the best alignment */

	int i=n,j=m;
// rpf	char c;
        char *r, *dir;
	int *d, *p, *q;


	d = val.first;
	p = val.second;
	q = val.third;
	r = (char *) calloc(n+m+1,sizeof(char));
	dir = r+n+m;    /* dir (direction) points to last cell of r */
	*dir = '\0';                /* last cell in r is null character */

	while ( i != 1 && j != 1 )
		{
		if ( D(i,j) == P(i,j)  )
			{
			i--;
			*(--dir) = 'u';   /* up */
			}
		else if ( D(i,j) == D(i-1,j-1) + GAMMA(i,j) )
			{
			i--; j--;
			*(--dir) = 'd';   /* diagonal */
			}
		else if ( D(i,j) == Q(i,j) )
			{
			j--;
			*(--dir) = 'l';   /* left */
			}
		}
	while ( i != 1 )
		{
		i--;
		*(--dir) = 'u';   /* up */
		}
	while ( j != 1 )
		{
		j--;
		*(--dir) = 'l';   /* left */
		}
	return(dir);
}


void display ( Triple val, char *s, char *t )
{
	/* display the distance matrix d(i,j) between strings
	   s,t where
			0<=i<=n = strlen(s)
	   and
			0<=j<=m = strlen(t) */

	int i,j, n = strlen(s), m = strlen(t);
	int *d, *p, *q;

	d = val.first;
	p = val.second;
	q = val.third;
	printf("Distance Matrix for alignment of %s with %s.\n\n",s,t);

	/****************** print letters of string t ***********/
	printf("%5c",SPACE);
	printf("%5c",SPACE);
	for (j=1;j<=m;j++) printf("%5c",T(j));
	printf("\n");


	/********************* print D(0,0),...,D(0,m) ************/
        printf("D matrix\n\n");
	printf("%5c",SPACE);
	for (j=0;j<=m;j++)
		printf("%5d",D(0,j));
	printf("\n");

	for (i=1;i<=n;i++)
		{
		printf("%5c",S(i));
		for (j=0;j<=m;j++)
			printf("%5d",D(i,j));
		printf("\n");
		}
	/********************* print P(0,0),...,P(0,m) ************/
        printf("P matrix\n\n");
	printf("%5c",SPACE);
	for (j=0;j<=m;j++)
		printf("%5d",P(0,j));
	printf("\n");


	for (i=1;i<=n;i++)
		{
		printf("%5c",S(i));
		for (j=0;j<=m;j++)
			printf("%5d",P(i,j));
		printf("\n\n\n");
		}
	/********************* print Q(0,0),...,Q(0,m) ************/
        printf("Q matrix\n\n");
	printf("%5c",SPACE);
	for (j=0;j<=m;j++)
		printf("%5d",Q(0,j));
	printf("\n");


	for (i=1;i<=n;i++)
		{
		printf("%5c",S(i));
		for (j=0;j<=m;j++)
			printf("%5d",Q(i,j));
		printf("\n");
		}
}


void dump_without_the_Z(char *name,char *note,char *s0,char *t0)
{
    register char *z;
    register char *z2;
    z = s0;

    printf("%s - %s\n",name,note); 
    for (z = s0;*z;z++)
    {
        if (*z == 'Z') continue;
        else if (*z >= 123) printf("*");
        else if (*z <= 32 ) printf("*");
        else                printf("%c",*z);
    }
    printf("\n"); 
    for (z = s0, z2=t0 ; (*z && *z2) ; z++,z2++)
    {
        if (*z == 'Z') continue;
        if (*z == *z2) printf("|");
        else           printf(" ");
    }
    printf("\n"); 

    for (z = t0;*z;z++)
    {
        if (*z == 'Z') continue;
        else if (*z >= 123) printf("*");
        else if (*z <= 32 ) printf("*");
        else                printf("%c",*z);
    }
    printf("\n\n"); 
    return;
}


void alignment ( char *direction, char *s, char *t ,char *name, char *note)
{
	/** Computes strings s0,t0 of length MAX(len(s),len(t)),
		where s0 [resp t0] is s [resp t] possibly with
		some interspersed '-' marks **/

	int n = strlen(s), m = strlen(t), 
	max = n+m;

	int i,j,k;
	char *s0, *t0;

	s0 = (char *) calloc(max,sizeof(char));
	t0 = (char *) calloc(max,sizeof(char));
	for (k=0;k<max;k++) s0[k] = t0[k] = '\0';  
		 /* initialize to null char */


	s0[0] = s[0];
	t0[0] = t[0];   /* We assume first letter of s,t same */
	i=1; j=1;       /* indices running through s, t */

	for (k=0;k<max-1;k++)   
		/* direction has at most max-1 components */
		{
		if ( direction[k] == 'd' )
			{
			s0[k+1] = s[i++];
			t0[k+1] = t[j++];
			}
		else if ( direction[k] == 'u' )
			{
			s0[k+1] = s[i++];
			t0[k+1] = '-';
			}
		else if ( direction[k] == 'l' )
			{
			s0[k+1] = '-';
			t0[k+1] = t[j++];
			}
		}

        dump_without_the_Z(name,note,s0,t0);
	// printf("%s\n%s\n",s0,t0);
}




static void drive_gotoh(char *name,char *s1, char *s2, char *note)
{
char tmps1[2048]; 
char tmps2[2048]; 
	void display(Triple val, char *, char *);
	void alignment(char*, char*, char*, char *, char *);
	char *backtrack(Triple, char *, char *, int, int );
	char *s, *t;
	Triple val;	
	size_t n,m;

sprintf(tmps1,"Z%sZ",s1); 
sprintf(tmps2,"Z%sZ",s2); 
        s = &tmps1[0]; 
        t = &tmps2[0]; 

	n = strlen(s);
	m = strlen(t);
	val = dist(s,t);
	alignment(backtrack(val,s,t,(int)n,(int)m),s,t,name,note);   

        if (val.first) free(val.first); 
        if (val.second) free(val.second); 
        if (val.third) free(val.third); 

#if 0
	display(val,s,t); 
	printf("%s\n",backtrack(val,s,t,n,m)); 
#endif

}

void chkhit3(char khr[], int loc, char seq[],int len, char *name)
{
    char t[512];
    int i;
    int k;
    int lo;
    int hi;

    lo = loc;
    hi = loc+len;
    fetchwrap2(khr,lo,hi);
    k = strlen(whatthere); 

    for  (i=0;i<k;i++)
         whatthere[i] = toupper(whatthere[i]);

    reverse_and_compliment(seq,t,len);
    drive_gotoh(name,seq, whatthere,(char *)"original");
//    drive_gotoh(name,t  , whatthere,"revcomp");

    return;
}
#endif

static char *rpf_format_fasta(const bam_header_t *header, const bam1_t *b)
{
    int i = 0;
    char cigar[512];
    const bam1_core_t *c = &b->core;

 // char m[512];
//        int j,i;
//        kstring_t mystr;

#if 0
        memset(&mystr,0,sizeof(mystr));
	mystr.l = mystr.m = 0; mystr.s = 0;

	ksprintf(&mystr, "%s\t%d\t", bam1_qname(b), c->flag);
	if (c->tid < 0) kputs("*\t", &mystr);
	else ksprintf(&mystr, "%s\t", header->target_name[c->tid]);
	ksprintf(&mystr, "%d\t%d\t", c->pos + 1, c->qual);
	if (c->n_cigar == 0) kputc('*', &mystr);
	else{
		for (i = 0; i < c->n_cigar; ++i)
			ksprintf(&mystr, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	} 
	kputc('\t', &mystr);
	if (c->mtid < 0) kputs("*\t", &mystr);
	else if (c->mtid == c->tid) kputs("=\t", &mystr);
	else ksprintf(&mystr, "%s\t", header->target_name[c->mtid]);

	ksprintf(&mystr, "%d\t%d\t", c->mpos + 1, c->isize);
	for (i = 0; i < c->l_qseq; ++i) kputc(rpf_bam_nt16_rev_table[bam1_seqi(s, i)], &mystr);
	kputc('\t', &mystr);
	if (t[0] == 0xff) kputc('*', &mystr);
	else for (i = 0; i < c->l_qseq; ++i) kputc(t[i] + 33, &mystr);
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
		ksprintf(&mystr, "\t%c%c:", key[0], key[1]);
		if (type == 'A') { ksprintf(&mystr, "A:%c", *s); ++s; }
		else if (type == 'C') { ksprintf(&mystr, "i:%u", *s); ++s; }
		else if (type == 'c') { ksprintf(&mystr, "i:%d", *s); ++s; }
		else if (type == 'S') { ksprintf(&mystr, "i:%u", *(uint16_t*)s); s += 2; }
		else if (type == 's') { ksprintf(&mystr, "i:%d", *(int16_t*)s); s += 2; }
		else if (type == 'I') { ksprintf(&mystr, "i:%u", *(uint32_t*)s); s += 4; }
		else if (type == 'i') { ksprintf(&mystr, "i:%d", *(int32_t*)s); s += 4; }
		else if (type == 'f') { ksprintf(&mystr, "f:%g", *(float*)s); s += 4; }
		else if (type == 'd') { ksprintf(&mystr, "d:%lg", *(double*)s); s += 8; }
		else if (type == 'Z' || type == 'H') { ksprintf(&mystr, "%c:", type); while (*s) kputc(*s++, &mystr); ++s; }
	}
// example @SOLEXA9_26:5:75:1325:688 chr10 123384948 50
/*
        i = c->pos+1 - start;
        if (i<0) i = 0;
        for (j=0;  (j<(c->l_qseq))  ; j++,i++)
        {  
            if (i>=maxspace) break;
             *(space+i) += 1;
        }
*/


#define BAM_OFDEC          0
	if (c->n_cigar == 0) 
        {
        }
	else {
                int i;
		for (i = 0; i < c->n_cigar; ++i)
		{
			sprintf(m, "cig:%s:%d%c", header->target_name[c->tid],
                               bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
		}
	}
#endif

#define bam1_opticaldupe(b) (((b)->core.flag&BAM_FDUP) != 0) 
#define bam1_badqual(b) (((b)->core.flag&BAM_FQCFAIL) != 0) 
#define bam1_unmapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)

        if (bam1_unmapped(b) == 1)
        {
            unmapcnt++;
// fprintf(stderr,"READ unmapped %ld of %ld  pos=%d percent=%f\n",unmapcnt,rdcnt,c->pos,(double)unmapcnt/(double)rdcnt);
            return 0;
        }

       // #define BAM_FQCFAIL      512
        if (bam1_opticaldupe(b) == 1)
            opticaldupecount++;
        else if (bam1_badqual(b) == 1) badqualcount++;
        else 
            if (c->n_cigar == 0) 
            {
               strcpy(cigar,"*");
               nocigarcount++;
            }
            else 
            {
                cigar[0] = (char)0;
                for (i = 0; i < c->n_cigar; i++)
                {
                    char tmps[512];
                    sprintf(tmps,"%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
                    strcat(cigar,tmps);
                }
                if (c->l_qseq >= 2048)
                {
                    paint_align_box(gdiff,c->pos+1-gs ,c->l_qseq,/*wcode*/0,/*notch*/-1,/*qual*/-1,/*fwd*/bam1_strand(b),/*bow*/0,0,0,0,
                             bam1_opticaldupe(b));
                }
                else
                {
 	            uint8_t *s = bam1_seq(b); // , *t = bam1_qual(b);
	            for (i = 0; i < c->l_qseq; i++) 
	            {
                        seq[i] = rpf_bam_nt16_rev_table[bam1_seqi(s, i)];
    	            }
                    seq[i] = (char)0;

                    if (global_code == 1) // fasta 
                    {
                        if (global_fastafastq_option ==  0) // fasta only 
                        {
/*
/h1/finneyr/samtools-0.1.7a/bam_stat.c:                 if ((c)->flag & BAM_FREAD1) ++(s)->n_read1;                                     \
/h1/finneyr/samtools-0.1.7a/bam_stat.c:                 if ((c)->flag & BAM_FREAD2) ++(s)->n_read2;  
*/
                            printf(">%s",bam1_qname(b));
                            if (((c)->flag & BAM_FREAD1)) printf(".1");
                            else if (((c)->flag & BAM_FREAD2)) printf(".2");
                            printf("\n");
                            printf("%s\n",seq);
                        }
                        else  // fastq
                        {
                         char  qual[512] ;
	                 uint8_t *t = bam1_qual(b);
                         i = 0;
                         if (c->l_qseq) 
                         {
	                     for (i = 0; i < c->l_qseq; i++) 
	                     {
                                qual[i] = t[i] + 33;
    	                     }
                         }
                         qual[i] = (char)0;

                            printf("@%s\n",bam1_qname(b));
                            printf("%s\n",seq);
                            printf("+\n");
                            printf("%s\n",qual);
                        }
                    }
#if GOTOH
                    else // if (global_code == 2)  gotoh 
                    {
                        chkhit3(global_khr, c->pos, seq,i,bam1_qname(b));
                    }
#endif
                }
            }

//	if (mystr.s) free(mystr.s);
        globalreadscount++;

	return (char *)0; // mystr.s;
}


#if 0
static int rpfwrite(samfile_t *fp, const bam1_t *b)
{
     rpf_format2(fp->header, b);
     return 0;
}
#endif


// callback function for bam_fetch() ** THIS GETS CALLED A LOT 
static int view_func(const bam1_t *b, void *data)
{
    if (!__g_skip_aln(((samfile_t*)data)->header, b))
        rpf_format2(((samfile_t*)data)->header, b);
    return 0;
}


static int rpfwrite_fasta(samfile_t *fp, const bam1_t *b)
{
     rpf_format_fasta(fp->header, b);
     return 0;
}

static int view_func_fasta(const bam1_t *b, void *data)
{
    if (!__g_skip_aln(((samfile_t*)data)->header, b))
        rpfwrite_fasta((samfile_t*)data, b);
    return 0;
}

static int rpfwrite_sam(samfile_t *fp, const bam1_t *b)
{
     rpf_format_sam(fp->header, b,samoption); // last arg ("of") = 2 means "string"
// of (output format) defines ....
// #define BAM_OFDEC          0
// #define BAM_OFHEX          1
// #define BAM_OFSTR          2
     return 0;
}
static int view_func_sam(const bam1_t *b, void *data)
{
jdebug("in view_func_sam() "); 
    if (!__g_skip_aln(((samfile_t*)data)->header, b))
        rpfwrite_sam((samfile_t*)data, b);
jdebug("end view_func_sam() "); 
    return 0;
}


char junkfile[512];


void bam_zap_old_files(void)
{
    int aflag = 0;
    lflag = 1;
    char m[512];
    DIR *d;
    struct dirent *e;
    time_t now;
    time_t timer;
    time_t ptimer;        // past
    char filename2zap[2048];
    char directory[2048];


    strcpy(directory, "/tmp" );

    timer = time(NULL);         // returns the time since 00:00:00 GMT, Jan. 1, 1970, measured in seconds.
    ptimer = timer - (60 );           // sec min hours -- cut off is 4 hours
    
    d = opendir(directory);
    if (! d) {
        jdebug("cant open /tmp dir in bamstuff.c"); 
        return;
    }

    now = time(NULL);
    
    while (1) 
    {
        e = readdir(d);
        if (! e) {
            break;
        }
        if (aflag || e->d_name[0] != '.') {
            if (lflag) 
            {
                struct stat s;
                if (stat(e->d_name, &s) == -1) {
                    fprintf(stderr, "%s: %s: ", progname, e->d_name);
                    perror(NULL);
                    continue;
                }
                           //123456789012
                if (strncmp("okay2delete.",e->d_name ,11) == 0)
                {
                    sprintf(filename2zap,"%s/%s",directory,e->d_name);
                    if (s.st_mtime < ptimer) 
                    {
                        sprintf(m,"zapping %s ",filename2zap);jdebug(m);
#ifdef WIN32
                        _unlink(filename2zap);
#else
                        unlink(filename2zap);
#endif
                    }
                }
            }
        }
    }
    closedir(d);
}


char region[512];
char region2[512];
char region3[512];
char region4[512]; // try MT for chrmM ?

void setregions(char chr[],int s,int e)
{
    int jj;
    int nus; // new s
    int nue; // new e
char m[512];

sprintf(m,"in setregions() , chr=[%s], s=%d e=%d ",chr,s,e); jdebug(m);
    nus = s-150; // look ahead
    if (nus<=0) nus=1;
    nue = e;
    sprintf(region,"%s:%d-%d",chr,nus,nue);
    jj =  atoi(&chr[3]);
    if (jj < 1)
    {
        sprintf(region2,"%d:%d-%d",jj,nus,nue);  // chr
        if (strcmp("chrM",chr) == 0) { jj = 25; sprintf(region3,"M:%d-%d",nus,nue); }
        else if (strcmp("chrX",chr) == 0) { jj = 23; sprintf(region3,"X:%d-%d",nus,nue); }
        else if (strcmp("chrY",chr) == 0) { jj = 24; sprintf(region3,"Y:%d-%d",nus,nue); }
        return;
    }
    sprintf(region2,"%d:%d-%d",jj,nus,nue);  // chr
    sprintf(region3,"%d:%d-%d",jj,nus,nue);  // chr
    return;
}

unsigned short int *gg1 = (unsigned short int *)0;
size_t cov_len;
int global_cov_flag = 0;


#if 0 // old 
int fopen_and_read_cov_OLD(char dbgmsg[], char fullfn[],unsigned short int *z, long int spot, int len) // multiplies by 2 !
{
    int err;
    char m[512];
    FILE *fp;
    int status;

    fp = fopen(fullfn,"rb");
    if (fp == (FILE *)0)
    {
         err = errno;
sprintf(m,"ERROR: cant open fullfn=[%s] in fopen_and_read_cov_OLD() %s, errno = %d",fullfn,dbgmsg,err); jdebug(m);
         memset(z,0,len*2);
         return -1;
    }
    status = fseek(fp,spot,SEEK_SET);
    if (status != 0)
    {
sprintf(m,"ERROR: cant seek %s in fopen_and_read_cov_OLD()",fullfn); jdebug(m);
         return -1;
    }
    status  = fread(z,2,len,fp);
    if (status != len)
    {
sprintf(m,"ERROR: cant read %s in fopen_and_read_cov_OLD()",fullfn); jdebug(m);
         return -1;
    }
// sprintf(m,"read %d * 2 from %s <br>\n",len,fullfn); jdebug(m);

    global_cov_flag = 1;
    fclose(fp);
    return 0;
}


void start_do_cov(char bamfn[], unsigned int circ_genome_s, unsigned int circ_genome_e)
{
    int px,py;
    double local_pxwidth;
    int status;
    long spot;
    char m[512];
    char tmps[512];

    if (has_cov_file(bamfn,covfn))
    {
        // strcpy(khr,chromName); circ_genome_s = seqStart;  get_circ_genome_addr(khr,&cir_genome_s);
        // strcpy(khr,chromName); circ_genome_3 = seqEnd;  get_circ_genome_addr(khr,&cir_genome_e);
        cov_len = (size_t)(circ_genome_e - circ_genome_s); 
        gg1 = (int *)malloc(cov_len*2); if (gg1 == (void *)0)  // gg1
        if (!gg1)
        { 
             sprintf(m,"cant malloc %ld bytes\n",cov_len);
             jdebug(m);
             return; 
        }

        spot = (long)circ_genome_s*(long)2;

        sprintf(tmps,"/tcga_next_gen/TCGAWIBR/rich/%s.cov",covfn);
        status = fopen_and_read_cov_OLD("from start_do_cov",tmps,gg1,spot,cov_len);
        if (status == -1) return;

        int y1;
        int i;
        int x;

        local_pxwidth = (double)iw/(double)(cov_len);
        px = py = -1;
        for (i=0 ; i<cov_len ; i++)
        {
            x = (int)(local_pxwidth * (double)i);
            y1 = ih - *(gg1+i);
            if (y1<0) y1 = 0;

            if (px == -1) gdImageLine (im, x,y1,x,y1+1,black);
            else  gdImageLine (im, px,py,x,y1,black);
            px = x;
            py = y1;
        }
    }
    return;
}
#endif


void set_junk_file(char junkfile[]) 
{
#ifdef WIN32
    strcpy(junkfile,"nul");
#else
    strcpy(junkfile,"/dev/null");
#endif
}


int global_bamerr = 0;

int dobam_fasta(char fn[],int khroffset,int s, int e,char chr[],int kstart, int kend, int code )
{
    bam_index_t *idx = 0;
    int searcht;
    int numrecs = 0;
    char m[512];
    int tid, beg, end;
    int is_header = 0, is_header_only = 0, is_bamin = 1, ret = 0, is_uncompressed = 0, is_bamout = 0;
    samfile_t *in = 0;
    samfile_t *out = 0;
    char in_mode[24], out_mode[24]; 
    char *fn_list = 0;
    char bai_fn[1024];  // BAm Index
    char *fn_out = 0;

    numrecs = searcht = 0;
    globalreadscount = 0;
    global_code = code;

    if (ge <= gs) ge = gs+1;
    setregions(chr,gs,ge);

sprintf(m,"in dobam_fasta() START func, region=%s",region); jdebug(m);

    set_junk_file(junkfile);    // set as null (/dev/null) file, OS specific, "nul" for windows
    fn_out = &junkfile[0];

    bai_fn[0] = (char)0;

             /* parse command-line options */
    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");

    if (is_uncompressed) is_bamout = 1;
    if (is_header_only) is_header = 1;
    if (is_bamout) strcat(out_mode, "b");
    if (is_bamin) strcat(in_mode, "b");
    if (is_header) strcat(out_mode, "h");
    if (is_uncompressed) strcat(out_mode, "u");
    strcpy(in_mode, "r"); strcpy(out_mode, "w");
    fn_list = fn;

    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");

    // open file handlers
    sprintf(bai_fn,"%s.bai",fn_list);
    if ((in = samopen(fn, in_mode, fn_list)) == 0) 
    {
        sprintf(m, "[dobam_fasta] fail to open sam/bam file for reading : %s",fn);
        jdebug(m);
        global_bamerr = 1;
        // cant write output here  printf(m,"<font color=\"red\">ERROR: failed to open %s <br>\n",fn);
        goto view_end;
    }

    if ((out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0)  // why do I even need "out" ?
     {
        sprintf(m, "[dobam_fasta] fail to open file for writing. fn=%s\n",fn_out);
        jdebug(m); 
        goto view_end;
    }
    if (is_header_only) goto view_end; // no need to print alignments
sprintf(m,"opened %s %s %s",fn_list,bai_fn,fn_out);  jdebug(m); 
        if (is_bamin) idx = bam_index_load(fn); // load BAM index - load 

        if (idx == 0) { // index is unavailable
            sprintf(m, "[main_samview] random alignment retrieval only works for indexed BAM files.\n");
            jdebug(m);
            ret = 1;
            goto view_end;
        }
sprintf(m,"in bamtest. after bam_index_load , index is %p, region = [%s]",idx,region); jdebug(m); 

        bam_parse_region(in->header, region , &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
        if (tid < 0)
        {
//sprintf(m,"in bamtest after bam_parse_region1 BAD , tid=%d from %s",tid,region); jdebug(m); 
            bam_parse_region(in->header, region2 , &tid, &beg, &end); // parse a region in the format like `2:100-200'
sprintf(m,"in bamtest after bam_parse_region2 ???, tid=%d from %s",tid,region2); jdebug(m); 
        }
        else
        {
sprintf(m,"in bamtest after bam_parse_region1 GOOD , tid=%d from %s",tid,region); jdebug(m); 
        }

        if (tid < 0) { // reference name is not found
            sprintf(m, "[dobam_fasta] fail to get the reference name. will exit. tid=%d [%s]",tid,fn);
            jdebug(m); 
            goto view_end;
        }

sprintf(m, "in dobam_fasta() before bamfetch %d %d %d VFF ",tid,beg,end); jdebug(m); 
        bam_fetch(in->x.bam, idx, tid, beg, end, out, view_func_fasta); // fetch alignments - "view_func_fasta" is callback 
sprintf(m,"in dobam_fasta() after bamfetch "); jdebug(m); 

        bam_index_destroy(idx); // destroy the BAM index

view_end:
        // close files, free and return

 sprintf(m,"in dobam_fasta, before samclose "); jdebug(m); 
    if (in) samclose(in);
    in = (samfile_t *)0;
    if (out) samclose(out);
    out = (samfile_t *)0;
/*
 sprintf(m,"in dobam_fasta, before samclose out %p",out); jdebug(m); 
    samclose(out);
 sprintf(m,"in dobam_fasta, before unlink"); jdebug(m); 
*/
#if 0 // not needed now that we use nul of /dev/null 
#ifdef WIN32
    _unlink(junkfile);
#else
    unlink(junkfile);
#endif
#endif 

#if WEB_SERVER
    bam_zap_old_files(); // zap old "okay2delete*" files in /tmp  just in case
#endif

    ret = globalreadscount;
    return ret;
}


int dobam_sam(char fn[],int khroffset,int s, int e,char chr[],int kstart, int kend, int code )
{
    bam_index_t *idx = 0;
    int searcht;
    int numrecs = 0;
    char m[512];
    int tid, beg, end;
    int is_header = 0, is_header_only = 0, is_bamin = 1, ret = 0, is_uncompressed = 0, is_bamout = 0;
    samfile_t *in = 0;
    samfile_t *out = 0;
    char in_mode[24], out_mode[24]; 
    char *fn_list = 0;
    char bai_fn[128];  // BAm Index

    char *fn_out = 0;
    numrecs = searcht = 0;
    globalreadscount = 0;
    global_code = code;

    if (ge <= gs) ge = gs+1;
    setregions(chr,gs,ge);

sprintf(m,"in dobam_sam() START func, region=%s",region); jdebug(m);

    set_junk_file(junkfile); 
    fn_out = &junkfile[0];

    bai_fn[0] = (char)0;
             /* parse command-line options */
    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");

    if (is_uncompressed) is_bamout = 1;
    if (is_header_only) is_header = 1;
    if (is_bamout) strcat(out_mode, "b");
    if (is_bamin) strcat(in_mode, "b");
    if (is_header) strcat(out_mode, "h");
    if (is_uncompressed) strcat(out_mode, "u");
    strcpy(in_mode, "r"); strcpy(out_mode, "w");
    fn_list = fn;

    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");

    // open file handlers
    sprintf(bai_fn,"%s.bai",fn_list);
    if ((in = samopen(fn, in_mode, fn_list)) == 0) 
    {
        sprintf(m, "[dobam_sam] fail to open file for reading : %s",fn);
        jdebug(m);
        global_bamerr = 2;
        // cant write output here  printf(m,"<font color=\"red\">ERROR: failed to open %s <br>\n",fn);
        goto view_end;
    }

    if ((out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) 
     {
        sprintf(m, "[dobam_sam] fail to open file for writing. fn=%s \n",fn_out);
        jdebug(m); 
        goto view_end;
    }
    if (is_header_only) goto view_end; // no need to print alignments
sprintf(m,"opened %s %s %s",fn_list,bai_fn,fn_out);  jdebug(m); 


// sprintf(m,"in bamtest . convert part "); jdebug(m); 
        if (is_bamin) idx = bam_index_load(fn); // load BAM index - load 

        if (idx == 0) { // index is unavailable
            sprintf(m, "[main_samview] random alignment retrieval only works for indexed BAM files.\n");
            jdebug(m);
            ret = 1;
            goto view_end;
        }
// sprintf(m,"in bamtest. after bam_index_load , index is %p, region = [%s]",idx,region); jdebug(m); 

        bam_parse_region(in->header, region , &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
        if (tid < 0)
        {
//sprintf(m,"in bamtest after bam_parse_region1 BAD , tid=%d from %s",tid,region); jdebug(m); 
            bam_parse_region(in->header, region2 , &tid, &beg, &end); // parse a region in the format like `2:100-200'
sprintf(m,"in bamtest after bam_parse_region2 ???, tid=%d from %s",tid,region2); jdebug(m); 
        }
        else
        {
sprintf(m,"in bamtest after bam_parse_region1 GOOD , tid=%d from %s",tid,region); jdebug(m); 
        }

        if (tid < 0) { // reference name is not found
            sprintf(m, "[dobam_sam] fail to get the reference name. will exit. tid=%d [%s]",tid,fn);
            jdebug(m); 
            goto view_end;
        }

sprintf(m,"in dobam_sam() before bamfetch %d %d %d",tid,beg,end); jdebug(m); 

/*******            HERE                              VVVVV at view_func_xxxxx()  big function call here ****/
        bam_fetch(in->x.bam, idx, tid, beg, end, out, view_func_sam); // fetch alignments - "view_func_sam" is callback 

        bam_index_destroy(idx); // destroy the BAM index

view_end:
        // close files, free and return

 sprintf(m,"in dobam_sam , before samclose in dobam_sam"); jdebug(m); 
    if (in) samclose(in);
    in = (samfile_t *)0;
    if (out) samclose(out);
    out = (samfile_t *)0;

 sprintf(m,"in dobam_sam , before unlink dobam_sam %s",junkfile); jdebug(m); 
/* don't need to delete NULL file
    if (junkfile[0]) 
    {
#ifdef WIN32
    _unlink(junkfile);
#else
    unlink(junkfile);
#endif
    }
*/

#if WEB_SERVER
    bam_zap_old_files(); // zap old "okay2delete*" files in /tmp  just in case
#endif

    ret = globalreadscount;
sprintf(m,"end dobam_sam %d",ret); jdebug(m); 
    return ret;
}


int dobam(char fn[],int khroffset,int s, int e,char chr[])
{
    bam_index_t *idx = 0;
    int errornumber; // holds errno
    int searcht;
    int numrecs = 0;
    char m[5120];
    int tid, beg, end;
    int is_header = 0, is_header_only = 0, is_bamin = 1, ret = 0, is_uncompressed = 0, is_bamout = 0;
    samfile_t *in = (samfile_t *)0;
    samfile_t *out = (samfile_t *)0;
    char in_mode[240], out_mode[24]; 
    char *fn_list = (char *)0;
    char bai_fn[1280];  // BAM Index ".bai" file
    char *fn_out = 0;

// debug fprintf(stderr,"in dobam(file=%s) %s %d %d",fn,chr,s,e);  fflush(stderr); 

    numrecs = searcht = 0;
    globalreadscount = 0;
    memset(line,0,5000);
    memset(line2,0,5000);
    memset(mms,0,5000);
    memset(mmsLQ,0,5000);
sprintf(m,"in dobam 0.1, before setregions %d %d",gs,ge);  jdebug(m); 

    gs = s;
    ge = e;
    if (ge<=gs) ge = gs+1;
    setregions(chr,gs,ge);
// sprintf(m,"in dobam 0.2");  jdebug(m); 

    pxwidth = (double)iw / (double)gdiff;

sprintf(m,"in dobam() START func, region=%s, pxw=%f (%d/%d)",region,pxwidth,iw,gdiff); jdebug(m);

    set_junk_file(junkfile); 
    fn_out = &junkfile[0];

    bai_fn[0] = (char)0;
             /* parse command-line options */
    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");

    if (is_uncompressed) is_bamout = 1;
    if (is_header_only) is_header = 1;
    if (is_bamout) strcat(out_mode, "b");
    if (is_bamin) strcat(in_mode, "b");
    if (is_header) strcat(out_mode, "h");
    if (is_uncompressed) strcat(out_mode, "u");
    strcpy(in_mode, "r"); strcpy(out_mode, "w");
    fn_list = fn;

    strcpy(in_mode, "rb"); strcpy(out_mode, "wb");
// sprintf(m,"in dobam 2.0 ");  jdebug(m); 
    // open file handlers
    sprintf(bai_fn,"%s.bai",fn_list);
    if ((in = samopen(fn, in_mode, fn_list)) == 0) 
    {
        sprintf(m, "[dobam] fail to open file for reading : [%s]",fn); jdebug(m);
        global_bamerr = 3;
sprintf(m,"fail to open [%s] ",fn); jdebug(m);
        // cant write output here  printf(m,"<font color=\"red\">ERROR: failed to open %s <br>\n",fn);
        goto view_end;
    }

    if ((out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) 
    {
        errornumber = errno;
        sprintf(m, "[dobam] fail to open file for writing. fn=%s, errno=%d ***\n",fn_out,errornumber);
        jdebug(m); 
        goto view_end;
    }
    if (is_header_only) goto view_end; // no need to print alignments
sprintf(m,"opened %s %s %s",fn_list,bai_fn,fn_out);  jdebug(m); 

// sprintf(m,"in bamtest . convert part "); jdebug(m); 
        if (is_bamin) idx = bam_index_load(fn); // load BAM index - load 

        if (idx == 0) { // index is unavailable
            sprintf(m, "[main_samview] random alignment retrieval only works for indexed BAM files.\n");
            jdebug(m);
            ret = 1;
            goto view_end;
        }
// sprintf(m,"in bamtest. after bam_index_load , index is %p, region = [%s]",idx,region); jdebug(m); 

        bam_parse_region(in->header, region , &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
        if (tid < 0)
        {
//sprintf(m,"in bamtest after bam_parse_region1 BAD , tid=%d from %s",tid,region); jdebug(m); 
            bam_parse_region(in->header, region2 , &tid, &beg, &end); // parse a region in the format like `2:100-200'
            if (tid < 0)
            {
                bam_parse_region(in->header, region3 , &tid, &beg, &end); //parse a region in the format like `2:100-200'
            }
sprintf(m,"in bamtest after bam_parse_region2 ???, tid=%d from %s",tid,region2); jdebug(m); 
        }
        else
        {
sprintf(m,"in bamtest after bam_parse_region1 GOOD , tid=%d from %s",tid,region); jdebug(m); 
        }

        if (tid < 0) { // reference name is not found
            sprintf(m, "[dobam] fail to get the reference name. will exit. tid=%d [%s]",tid,fn);
            jdebug(m); 
            goto view_end;
        }

// fprintf(stderr,"in do_bam before bam-fetch viewfunc \n");   fflush(stderr);  
// sprintf(m, "in dobam() before bamfetch tid=%d beg=%d end=%d out=%p vf=%p VF",tid,beg,end,out,view_func); jdebug(m); 
        bam_fetch(in->x.bam, idx, tid, beg, end, out, view_func); // fetch alignments - "view_func" is callback 
// sprintf(m, "in dobam() after view_func() globalreadscount=%d",globalreadscount); jdebug(m); 

        bam_index_destroy(idx); // destroy the BAM index

view_end:
        // close files, free and return

sprintf(m,"in dobam , after samclose in dobam() (before unlink()?) "); jdebug(m); 
    if (in) samclose(in);
    in = (samfile_t *)0;
    if (out) samclose(out);
    out = (samfile_t *)0;

#if 0
// not needed even though nul or /dev/null ??
#ifdef WIN32
    _unlink(junkfile);
#else
    unlink(junkfile);
#endif
#endif

// sprintf(m,"rpf: in bam view end, reads = %d maxspace = %d",globalreadscount,maxspace); jdebug(m);
#if WEB_SERVER
    bam_zap_old_files(); // zap old "okay2delete*" files in /tmp  just in case
#endif
    ret = globalreadscount;
 sprintf(m,"end dobam ret(globalreadscount)=%d ",ret); jdebug(m); 
    return ret;
}


int do_fasta( int start ,int end, int image_type, char shortfilename[],char chr[],int khrstart, int khrend)
{
    int ki;
    int status;
    int e = 0;
    int s = 0;
    char m[2048];

sprintf(m,"in do_fasta chr=%s start=%d end=%d",chr,start,end); jdebug(m); 
    s = start;
    e = end;
    status = get_chr_lo_hi(chr,&s,&e,&ki);
    if (status == -1)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d in do_fasta()",chr,start,end); 
        jdebug(m);
        return -1;
    }

    dobam_fasta(fn_bam,ki,s,e,chr,start,end,1);

sprintf(m,"end do_fasta (fn_bam=%s)",fn_bam);  jdebug(m); 
    return 0;
}


int do_sam( int start ,int end, int image_type, char shortfilename[],char chr[],int khrstart, int khrend)
{
    int ki;
    int status;
    int e = 0;
    int s = 0;
    char m[2048];

sprintf(m,"in do_sam chr=%s start=%d end=%d",chr,start,end); jdebug(m); 

    s = start;
    e = end;
    status = get_chr_lo_hi(chr,&s,&e,&ki);
    if (status == -1)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d in do_sam()",chr,start,end); 
        jdebug(m);
        return -1;
    }

    dobam_sam(fn_bam,ki,s,e,chr,start,end,1);

sprintf(m,"end do_sam (fn_bam=%s)",fn_bam);  jdebug(m); 
    return 0;
}


int do_gotoh( int start ,int end, int image_type, char shortfilename[],char chr[],int khrstart, int khrend)
{
    int ki;
    int status;
    int e = 0;
    int s = 0;
    char m[2048];

sprintf(m,"in do_gotoh chr=%s start=%d end=%d",chr,start,end); jdebug(m); 

    s = start;
    e = end;
    status = get_chr_lo_hi(chr,&s,&e,&ki);
    if (status == -1)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d in go_gotoh()",chr,start,end); 
        jdebug(m);
        return -1;
    }

    dobam_fasta(fn_bam,ki,s,e,chr,start,end,2);

sprintf(m,"end do_gotoh (fn_bam=%s)",fn_bam);  jdebug(m); 
    return 0;
}


void draw_mapability_lines_at_bottom(int khroffset,int s, int e,char chr[],int kstart, int kend, int readlen, 
                  int yoffset, int kolor1,  int kolor2)
{
    int x1,x2;
    int stat;
    int runlen;
    unsigned char uc;
    int i;
    long l;
    int diff;
    FILE *fp;
    char m[1024];
    char fn[1024]; // filename

sprintf(m,"in draw_mapability_lines_at_bottom() 000 , khroffset=%u s=%d e=%d chr=%s kstart=%d kend=%d",khroffset,s,e,chr,kstart,kend); jdebug(m);
    diff = e - s;
    pxwidth = (double)iw/(double)diff;

    sprintf(fn,"/app1/cgwb-shortread/COV/%s.all.%d",chr,readlen);
    fp = fopen(fn,"rb");
    if (fp == (FILE *)0) 
    {
        sprintf(m,"ERROR: in draw_mapability_lines_at_bottom() , can't open [%s]",fn); jdebug(m);
        return;
    }
    l = (long)kstart;
    if (fseek(fp,l,SEEK_SET) != 0)  // listart by earlier function calls
    {
        sprintf(m,"ERROR in draw_mapability_lines_at_bottom(), fseek(%p) to %ld failed ",fp,l);  jdebug(m); 
        return;
    }
    i = 0;
    runlen = kend - kstart;
    while (1)
    {
        if (i>=runlen) break;
        uc = 0xff;
        stat = fread(&uc,1,1,fp); 
        if (stat != 1) break;
        x1 = (int)(pxwidth * (double)i); 
        x2 = (int)(pxwidth * (double)(i+1)); 
        if (uc & (0x1 << 1))
        {
#if 0 // USE_GD
           gdImageFilledRectangle(im,x1,ih-3-(yoffset*4),x2,ih-2-(yoffset*4),kolor1);
#endif
           ImageFilledRectangle(im,x1,ih-3-(yoffset*4),x2,ih-2-(yoffset*4),kolor1);
        }
        if (uc & (0x1 << 0))
        {
#if 0 // USE_GD
           gdImageFilledRectangle(im,x1,ih-5-(yoffset*4),x2,ih-4-(yoffset*4),kolor2);
#endif
           ImageFilledRectangle(im,x1,ih-5-(yoffset*4),x2,ih-4-(yoffset*4),kolor2);

        }
// sprintf(m,"%d got %x from freading draw_mapability_lines_at_bottom() ",i,uc);  jdebug(m); 
        i++;
    }
    fclose(fp); 
sprintf(m,"end draw_mapability_lines_at_bottom() 999 "); jdebug(m);

    return;
}


void draw_mapability_lines_at_bottom2(int khroffset,int s, int e,char chr[],int kstart, int kend, int readlen, 
                  int yoffset)
{
    int y1,y2;
    int x1,x2;
    int stat;
    int runlen;
    unsigned char uc;
    int i;
    long l;
    int diff;
    FILE *fp;
    char m[1024];
    char fn[1024]; // filename
    int kolor = blue;

    diff = e - s;
    pxwidth = (double)iw/(double)diff;

    sprintf(fn,"/app1/cgwb-shortread/COV/%s.all.%d",chr,readlen);
    fp = fopen(fn,"rb");
    if (fp == (FILE *)0) 
    {
        sprintf(m,"ERROR: in draw_mapability_lines_at_bottom2(), can't open [%s]",fn); jdebug(m);
        return;
    }
    l = (long)kstart;
    if (fseek(fp,l,SEEK_SET) != 0)  // listart by earlier function calls
    {
        sprintf(m,"ERROR in draw_mapability_lines_at_bottom2(), fseek(%p) to %ld failed ",fp,l);  jdebug(m); 
        return;
    }
    i = 0;
    runlen = kend - kstart;
    while (1)
    {
        if (i>=runlen) break;
        uc = 0xff;
        stat = fread(&uc,1,1,fp); 
        if (stat != 1) break;

        x1 = (int)(pxwidth * (double)i); 
        x2 = (int)(pxwidth * (double)(i+1)); 
        if (uc == 1) /* ns */    kolor = white;
        else if (uc == 3) /* ns+s*/  kolor = lightblue;
        else if (uc == 2) /* s */   kolor = blue;
        else                         kolor = pink;
        y1 = ih-(yoffset*4)-2; 
        y2 = y1 + 3; 
#if 0 // USE_GD
        gdImageFilledRectangle(im,x1,y1,x2,y2,kolor);
#endif
        ImageFilledRectangle(im,x1,y1,x2,y2,kolor);

// sprintf(m,"%d got %x from freading draw_mapability_lines_at_bottom2() ",i,uc);  jdebug(m); 
        i++;
    }
    fclose(fp); 
    return;
}


void fix_num_to_chrstyle(char chr[],char chrchr[]) 
{
    if (strcmp(chr,"1")== 0) strcpy(chrchr,"chr1"); 
    else if (strcmp(chr,"2")== 0) strcpy(chrchr,"chr2"); 
    else if (strcmp(chr,"3")== 0) strcpy(chrchr,"chr3"); 
    else if (strcmp(chr,"4")== 0) strcpy(chrchr,"chr4"); 
    else if (strcmp(chr,"5")== 0) strcpy(chrchr,"chr5"); 
    else if (strcmp(chr,"6")== 0) strcpy(chrchr,"chr6"); 
    else if (strcmp(chr,"7")== 0) strcpy(chrchr,"chr7"); 
    else if (strcmp(chr,"8")== 0) strcpy(chrchr,"chr8"); 
    else if (strcmp(chr,"9")== 0) strcpy(chrchr,"chr9"); 
    else if (strcmp(chr,"10")== 0) strcpy(chrchr,"chr10"); 
    else if (strcmp(chr,"11")== 0) strcpy(chrchr,"chk11"); 
    else if (strcmp(chr,"12")== 0) strcpy(chrchr,"chr12"); 
    else if (strcmp(chr,"13")== 0) strcpy(chrchr,"chr13"); 
    else if (strcmp(chr,"14")== 0) strcpy(chrchr,"chr14"); 
    else if (strcmp(chr,"15")== 0) strcpy(chrchr,"chr15"); 
    else if (strcmp(chr,"16")== 0) strcpy(chrchr,"chr16"); 
    else if (strcmp(chr,"17")== 0) strcpy(chrchr,"chr17"); 
    else if (strcmp(chr,"18")== 0) strcpy(chrchr,"chr18"); 
    else if (strcmp(chr,"19")== 0) strcpy(chrchr,"chr19"); 
    else if (strcmp(chr,"20")== 0) strcpy(chrchr,"chr20"); 
    else if (strcmp(chr,"21")== 0) strcpy(chrchr,"chr21"); 
    else if (strcmp(chr,"22")== 0) strcpy(chrchr,"chr22"); 
    else if (strcmp(chr,"X")== 0) strcpy(chrchr,"chrX"); 
    else if (strcmp(chr,"Y")== 0) strcpy(chrchr,"chrY"); 
    else if (strcmp(chr,"M")== 0) strcpy(chrchr,"chrM"); 
    else if (strcmp(chr,"MT")== 0) strcpy(chrchr,"chrM"); 
}


int imgen( char shortfilename[],char chr[],int khrstart, int khrend,int image_type)
{
    int s,e;
    int len;
    int chrom_index;
    int status;
    char chrchr[1024];  
    char fn[1024];
    char m[2048];


    strcpy(chrchr,chr);
    fix_num_to_chrstyle(chr,chrchr); 

#ifdef CMD_LINE
strcpy(fn,shortfilename);
#endif

    s = e = len = status = 0;
    totalcnt = 0;
    len = khrstart - khrend;
    if (len > MAX_CHROMSIZE)
    {
sprintf(m,"ERROR: in imgen() size is too big "); jdebug(m);
        return -1;
    }

    pxwidth = (double)iw /(double)(len);
    gs = khrstart;
    ge = khrend;
    gdiff = ge - gs;
    ppp_firsttime = 1; 

sprintf(m,"start imgen() iw=%d ih=%d pxw=%f file=%s",iw,ih,pxwidth,shortfilename); jdebug(m);
#if 0 // USE_GD
    im = gdImageCreate(iw,ih);
    if (im == (void *)0)
    {
        sprintf(m,"ERROR: in imgen, NO im, gdImageCreate failed");  jdebug(m);
        return -1;
    }
#endif
    im = &image;
    im->width = iw;
    im->height = ih;
    im->data = (unsigned char *)malloc(iw*ih*3);
    if (im->data == (void *)0)
    {
        sprintf(m,"ERROR: in imgen, NO im, image create failed");  jdebug(m);
        return -1;
    }

    inited_colors_flag = 0; 
    init_colors(); //  inited_colors_flag  may already be set  

// sprintf(m,"in imgen() 2 %d %d ",iw,ih); jdebug(m);
    status = get_chr_lo_hi(chr,&s,&e,&chrom_index);

// sprintf(m,"in imgen after get_chr_lo_hi(), status= %d: chr=%s gs:%d ge:%d [chrom_index=%d] %s",status,chr,gs,ge,chrom_index,blds); jdebug(m);

    if (status == -1)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d %s in imgen()",chr,gs,ge,blds); 
        jdebug(m);
        return -1;
    }
#if 0 // USE_GD
    gdImageFilledRectangle(im,0,0,iw,ih,black);
    gdImageFilledRectangle(im,1,1,iw-2,ih-2,verylightgray);
#endif
    ImageFilledRectangle(im,0,0,iw,ih,black);
    ImageFilledRectangle(im,1,1,iw-2,ih-2,verylightgray);

    setup_reference_dna(chrchr,gs,ge,1); // importantly, get 2 bit data 
    setup_dnacnts_and_dnamms(ge-gs+1);
    dobam(fn_bam,chrom_index,gs,ge,chr);

    draw_dnacnts_and_dnamms(chr,gs,ge,24);
    if (geneannot_flag == 1) 
    {
// sprintf(m,"in imgen() before paint_gene_annot() ");  jdebug(m); 
        paint_gene_annot(chrchr,gs,ge);
// sprintf(m,"in imgen() after paint_gene_annot() ");  jdebug(m); 
    }

//    old start_do_cov(fn_bam, s,e);

    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,36,1);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,50,2);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,75,3);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,100,4);
// sprintf(m,"in imgen() after draw_mapability_lines_at_bottom2's ");  jdebug(m); 

// put notch in middle
#if 0 // USE_GD
    gdImageFilledRectangle(im, (iw/2)-2, ih-4, (iw/2)+2, ih-1, red);
#endif
    ImageFilledRectangle(im, (iw/2)-2, ih-4, (iw/2)+2, ih-1, red);

sprintf(m,"in imgen() in imgen() after gdImageFilledRectangle");  jdebug(m); 

#if 0 // QT_GUI
    strcpy(shortfilename,"tmp.png");
#ifdef WIN32
    sprintf(fn,"%s\\%s",cwd,shortfilename);   // NOTE NOT URL STYLE, DOS STYLE
    _unlink(fn);
#else
    sprintf(fn,"%s/%s",cwd,shortfilename);
    unlink(fn);
#endif

#endif


#if WEB_SERVER
#ifdef WIN32
    sprintf(shortfilename,"tmpalview.%d.png",rand());
#else
    sprintf(shortfilename,"tmpalview.%d.%d.png",getpid(),rand());
#endif

    sprintf(fn,"%s/%s",HTDOCSDIR,shortfilename); // must put in place apache/nginx can see it
sprintf(m,"WEB_SERVER: in imgen before savefile fn = [%s] ",fn); jdebug(m);
sprintf(m,"WEB_SERVER: HTDOCSDIR = [%s] ",HTDOCSDIR); jdebug(m);
#endif

sprintf(m,"in imgen before savefile fn = [%s] ",fn); jdebug(m);
    savefile(fn,im,GIF_IMAGE); // **********************************

sprintf(m,"in imgen() end imgen (fn=%s), returning status = 0 ",fn);  jdebug(m); 
    return 0;
}


int find_build(char fn_bam_arg[])
{
   char *z;
   int ret = 19; // hg18 is default 
   samfile_t space;
   samfile_t *inf = 0;
char m[512];


sprintf(m,"in find_build [%s] " , fn_bam_arg);  jdebug(m); 
   strcpy(blds,"hg19"); // default 
   if (fn_bam_arg[0] == (char)0) return 18;
sprintf(m,"in find_build 2" );  jdebug(m); 

   memset(&space,0,sizeof(space));
   inf = &space;

#if SAMTOOLS1
#else
// "type" field is not there in samtools 1.0 which is new as of September 2014
#define TYPE_BAM 1
   inf->type |= TYPE_BAM;
// 
#endif
   inf->x.bam = bam_open(fn_bam_arg, "r");
   if (!inf->x.bam) return -1;

   inf->header = bam_header_read(inf->x.bam);
// sprintf(m,"in find_build 3 [%s]",inf->header );  jdebug(m); 

   if (inf->header) 
   {
        z = inf->header->text;
        if (strstr(z,"NCBI-Build-36")) { ret = 18;  }
        if (strstr(z,"hg18")) { ret = 18;  }
        if (strstr(z,"NCBI-human-build36")) { ret = 18;  }
        if (strstr(z,"AS:HG18")) { ret = 18;  }
        if (strstr(z,"NCBI36")) { ret = 18;  }
        if (strstr(z,"LN:247249719")) { ret = 18;  }

        if (strstr(z,"LN:249250621")) ret = 19;
        if (strstr(z,"hg19")) { ret = 19;  }
        if (strstr(z,"GRCh37")) { ret = 19;  }
        if (strstr(z,"LB:HG19")) { ret = 19;  }
        if (strstr(z,"hg19")) { ret = 19;  }
        if (strstr(z,"GRCh37")) { ret = 19;  }
        if (strstr(z,"LB:HG19")) { ret = 19;  }

// new 2012 november , added hack support for variables sized output from novoalign
        if (strstr(inf->header->text,"SamTranscriptomeParser")) ret = 19;
        bam_header_destroy(inf->header);
   }
sprintf(m,"in find_build 4 ");  jdebug(m); 

   if (inf->x.bam)bam_close(inf->x.bam);
   if (ret == 19) strcpy(blds,"hg19"); 
   else if (ret == 18) strcpy(blds,"hg18"); 
   return ret;
}


unsigned char *imgen_mem(char fn[], char chr[],int khrstart, int khrend, int h, int w,int *ret_status) 
{
    int mallocsize = 0;
    int stati;
    int s,e;
    int len;
    int chrom_index;
    int status;
    char chrchr[1024];  
    char m[2048];


sprintf(m,"in imgen_mem 0"); jdebug(m); 
 sprintf(m,"in imgen_mem 0.1 s=%d",khrstart); jdebug(m); 
 sprintf(m,"in imgen_mem 0.2 e=%d",khrend); jdebug(m); 
 sprintf(m,"in imgen_mem 0.3 h=%d",h); jdebug(m); 
 sprintf(m,"in imgen_mem 0.4 w=%d",w); jdebug(m); 
 sprintf(m,"in imgen_mem 0.5 chr=%s",chr); jdebug(m); 
 sprintf(m,"in imgen_mem 0.6 fn=%p",fn); jdebug(m); 
 sprintf(m,"in imgen_mem 0.7 fn=%s",fn); jdebug(m); 

    if (fn[0] == (char)0)
    {
        *ret_status = -86;
        return (unsigned char *)0;
    }

    if (&fn_bam[0] != &fn[0]) 
        strcpy(fn_bam,fn);
    stati = find_build(fn_bam); // 18=hg18, 19=hg19, else -1 is error , sets global string variable "blds"

 sprintf(m,"in imgen_mem 0.6 stati=%d",stati); jdebug(m); 

    iw = w;
    ih = h;
    strcpy(chrchr,chr);
    fix_num_to_chrstyle(chr,chrchr); 

// sprintf(m,"in imgen_mem 1"); jdebug(m); 
    s = e = len = status = 0;
    totalcnt = 0;
    len = khrend - khrstart;
    if (len > MAX_CHROMSIZE)
    {
        sprintf(m,"ERROR: in imgen_mem() size is too big "); jdebug(m);
        *ret_status = -1;
        return (unsigned char *)0;
    }

    if (len <= 0)
    {
        sprintf(m,"ERROR: in imgen_mem() size is too small "); jdebug(m);
        *ret_status = -1;
        return (unsigned char *)0;
    }

// sprintf(m,"in imgen_mem 2"); jdebug(m); 
    pxwidth = (double)iw /(double)(len);
    gs = khrstart;
    ge = khrend;
    gdiff = ge - gs;
    ppp_firsttime = 1; 
// sprintf(m,"in imgen_mem 3"); jdebug(m); 

sprintf(m,"start imgen_mem() iw=%d ih=%d pxw=%f file=%s",iw,ih,pxwidth,fn_bam); jdebug(m);

#if 0 // USE_GD
    im = gdImageCreate(iw,ih);
    if (im == (void *)0)
    {
        sprintf(m,"ERROR: in imgen_mem, NO im, gdImageCreate failed");  jdebug(m);
        *ret_status = -2;
        return (unsigned char *)0;
    }
#endif
    im = &image;
sprintf(m,"in imgen_mem() 1 w=%d h=%d ",iw,ih); jdebug(m);
    mallocsize = (iw*ih*3);
    im->data = (unsigned char *)malloc(mallocsize);
    im->width = iw;
    im->height = ih;
sprintf(m,"in imgen_mem() 2 %d %d ,malloced %d",iw,ih,mallocsize); jdebug(m);

    if (im->data == (void *)0)
    {
        sprintf(m,"ERROR: in imgen_mem, NO im, image create failed");  jdebug(m);
        *ret_status = -2;
        return (unsigned char *)0;
    }
    memset(im->data,0,mallocsize);

sprintf(m,"in imgen_mem() 3 %d %d ",iw,ih); jdebug(m);
    init_colors();

sprintf(m,"in imgen_mem() 4 %d %d ",iw,ih); jdebug(m);
    status = get_chr_lo_hi(chr,&s,&e,&chrom_index);

sprintf(m,"in imgen_mem after get_chr_lo_hi(), status= %d: chr=%s gs:%d ge:%d [chrom_index=%d] %s",status,chr,gs,ge,chrom_index,blds); 
jdebug(m);

    if (status == -1)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d %s in imgen_mem()",chr,gs,ge,blds); 
        jdebug(m);
        *ret_status = -3;
        return (unsigned char *)0;
    }
#if 0 // USE_GD
    gdImageFilledRectangle(im,0,0,iw,ih,black);
    gdImageFilledRectangle(im,1,1,iw-2,ih-2,verylightgray);
#endif
    ImageFilledRectangle(im,0,0,iw,ih,black);
    ImageFilledRectangle(im,1,1,iw-2,ih-2,verylightgray);

    setup_reference_dna(chrchr,gs,ge,1); // importantly, get 2 bit data 
    setup_dnacnts_and_dnamms(ge-gs+1);
    dobam(fn_bam,chrom_index,gs,ge,chr);

    draw_dnacnts_and_dnamms(chr,gs,ge,24);

    if (geneannot_flag == 1) 
    {
        paint_gene_annot(chrchr,gs,ge);
    }

//    old start_do_cov(fn_bam, s,e);

    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,36,1);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,50,2);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,75,3);
    draw_mapability_lines_at_bottom2(chrom_index,gs,ge,chr,gs,ge,100,4);
sprintf(m,"in imgen() after draw_mapability_lines_at_bottom2's ");  jdebug(m); 

// put notch in middle
#if 0 // USE_GD
    gdImageFilledRectangle(im, (iw/2)-2, ih-4, (iw/2)+2, ih-1, red);
#endif
    ImageFilledRectangle(im, (iw/2)-2, ih-4, (iw/2)+2, ih-1, red);

sprintf(m,"in imgen_mem(after gdImageFilledRectangle, did notch, near end");  jdebug(m); 

#if 0
nope !!!!
    freedom_for_memory();
#endif

sprintf(m,"end imgen_mem() , will return %p",im->data);  jdebug(m); 
    *ret_status = 0; // good 
    return im->data;

}



int load_file_into_memory(char *fn,long int expected_size, char *p)
{
    long int li;
    FILE *fp;

fprintf(stderr,"Loading %ld bytes from %s\n",expected_size,fn); fflush(stderr);
    fp = fopen(fn,"rb");
    if (fp == (FILE *)0)
    {
        fprintf(stderr,"Can not OPEN file %s \n",fn);
        return -1;
    }
    if (fseek(fp,0,SEEK_END) != 0)
    {
        fprintf(stderr,"Can not fseek in file %s \n",fn);
        return -1;
    }
    li = ftell(fp);
    if (li != expected_size)
    {
        fclose(fp);
        fprintf(stderr,"ERROR: reading file %s expected_size is not equal to size (%ld != %ld)\n",
                         fn,expected_size,li);
        return -1;
    }
    rewind(fp);
    fread(p,li,1,fp);
    
/*
    li = 0L;
    while (1)
    {
        c = getc(fp);
        if (c == -1) break;
        *(p+li++) =  c;
    }
*/
    fclose(fp);
fprintf(stderr,"Loaded %ld bytes from %s\n",li,fn); fflush(stderr);
    return 0;
}




void reverse(char s[], char outs[])
{
    register int i,j;

    i = (int)strlen(s);
    i--;

    j = 0;
    while (i >= 0) 
    {
        outs[j] = s[i];
        j++;
        i--;
    }
    outs[j] = (char)0;
    strcpy(s,outs);
    return;
}

void compliment(char s[], char outs[])
{
    register int i,j;

    i = (int)strlen(s);
    i--;

    j = 0;
    while (i >= 0) 
    {
        outs[j] = s[i];
        j++;
        i--;
    }
    outs[j] = (char)0;
    return;
}



char revmap[256];


#define MAXSOLEXAS 6335899
// 5935899

struct solexas_type
{
   char *name;
   char *seq;
   int  pm;
   int  mm;
   int  mm2;
};

struct solexas_type solexas[MAXSOLEXAS];
int numsolexas=0;

int cmp_solexas(const void *arg1, const void *arg2)
{
    struct solexas_type *z1;
    struct solexas_type *z2;

    z1 = (struct solexas_type  *)arg1;
    z2 = (struct solexas_type  *)arg2;

    return strcmp(z1->name,z2->name);
}










char jobstring[MAXBUFF];
int alv_start,alv_end;



#define MAXCHRS 26

struct range_type
{
    unsigned int lo;
    unsigned int hi;
};

int numchrranges = 25;

struct range_type chr_ranges[MAXCHRS] =
{
    { 0 , 247249703 } ,   /* chr1  */
    { 247249703 , 490200836 } ,   /* chr2  */
    { 490200836 , 689702647 } ,   /* chr3  */
    { 689702647 , 880975694 } ,   /* chr4  */
    { 880975694 , 1061833544 } ,   /* chr5  */
    { 1061833544 , 1232733520 } ,   /* chr6  */
    { 1232733520 , 1391554928 } ,   /* chr7  */
    { 1391554928 , 1537829738 } ,   /* chr8  */
    { 1537829738 , 1678102974 } ,   /* chr9  */
    { 1678102974 , 1813477695 } ,   /* chr10  */
    { 1813477695 , 1947930063 } ,   /* chr11  */
    { 1947930063 , 2080279581 } ,   /* chr12  */
    { 2080279581 , 2194422545 } ,   /* chr13  */
    { 2194422545 , 2300791114 } ,   /* chr14  */
    { 2300791114 , 2401130013 } ,   /* chr15  */
    { 2401130013 , 2489957251 } ,   /* chr16  */
    { 2489957251 , 2568731977 } ,   /* chr17  */
    { 2568731977 , 2644849114 } ,   /* chr18  */
    { 2644849114 , 2708660749 } ,   /* chr19  */
    { 2708660749 , 2771096697 } ,   /* chr20  */
    { 2771096697 , 2818041004 } ,   /* chr21  */
    { 2818041004 , 2867732420 } ,   /* chr22  */
    { 2867732420 , 3022646158 } ,   /* chr23  */
    { 3022646158 , 3080419096 } ,   /* chr24  */
    { 3080419096 , 3080435651 } ,   /* chr25  */
    { 0 , 0 }     /* one overflew the chr-chr nest  */
};


char *chromnames[MAXCHRS];
int numchromoffsets = 0;


int bsearch_ranges( unsigned int key,  unsigned int *lo)
{
    size_t l, u, idx;

    l = 0;
    u = numchromoffsets;
    while (l < u) 
    {
        idx = (l + u) / 2;
        if ( (key>=chr_ranges[idx].lo) && (key<chr_ranges[idx].hi) )
        {
           *lo = chr_ranges[idx].lo;
           return (int)idx;
        }
        else if (key < chr_ranges[idx].lo)
           u = idx;
        else 
            l = idx + 1;
    } 
    return -1;
} 


void do_gotgenes()
{
    int i;
    char tmps[512];
    char q[51200];

    char m[512];

sprintf(m,"in do_gotgenes "); jdebug(m); 
    q[0] = (char)0;
#if KEGGINATOR
    if (gotten_genes_index)
    {
        printf("<a href=\"../cgi-bin/kg?txtarea=");
        for (i=0;i<gotten_genes_index;i++)
        {
           if (i) printf("%%0D%%0A");
           printf("%s",gotgenes[i]);
        }
        printf("&Submit=Submit\">Kegginator! </a> | \n");
    }
#endif
    if (gotten_genes_index > 0)
        printf("NCBI: ");
    for (i=0 ; i<gotten_genes_index ; i++)
    {
        if (i > 0)
        {
             sprintf(tmps,"%%20OR%%20%s",gotgenes[i]);
             strcat(q,tmps);
        }
        else
        {
             sprintf(tmps,"%%28%%20%s",gotgenes[i]); // "( "
             strcat(q,tmps);
        }
        printf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s%%20AND%%20Human%%5BOrganism%%5D\" target=\"_blank\">%s</a> | \n",gotgenes[i],gotgenes[i]);
    }

    if (gotten_genes_index)
    {
        sprintf(tmps,"%%29%%20AND%%20Human%%5BOrganism%%5D");  // ") AND Human[Organism]"
        strcat(q,tmps);
        printf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s\" target=\"_blank\">all</a> | \n",q);
    }
      
sprintf(m,"in do_gotgenes end - gotten_genes_index = %d ",gotten_genes_index); jdebug(m); 
    printf("<br>");
}

void do_gotgenes_trawler()
{
    int i;

    if (gotten_genes_index > 0)
        printf("TRAWLER: ");
    for (i=0 ; i<gotten_genes_index ; i++)
    {
// https://cgwb.nci.nih.gov/cgi-bin/trawler?q=TP53
        printf("<a href=\"https://cgwb.nci.nih.gov/cgi-bin/trawler?q=%s\" target=\"_blank\">%s</a> | \n",gotgenes[i],gotgenes[i]);
    }
    printf("<br>");
}


void fix_doubleslash_to_singleslash(char fixme[])
{
    int i;
     char t[2010];

     for (i=0;fixme[i];i++)
     {
         if ((fixme[i] == '/') && (fixme[i+1] == '/')) 
         {
             strcpy(t,&fixme[i+1]); 
             strcpy(&fixme[i],t);
         }
     }
    return;
}


int isT(char *s)
{
    char *z;
    z = strstr(s,"TCGA");
    if (z) 
    {
        if (strncmp(z+13,"01",2) == 0) return 1;
        if (strncmp(z+13,"02",2) == 0) return 1;
    }
    return 0;
}

int isN(char *s)
{
    char *z;
    z = strstr(s,"TCGA");
    if (z) 
    {
        if (strncmp(z+13,"11",2) == 0) return 1;
        if (strncmp(z+13,"10",2) == 0) return 1;
    }
    return 0;
}


char *getbasename(char s[])
{
    size_t i;

    i = strlen(s);
    if (i <= 0) return (char *)0;
    while (i>0)
    {
        if (s[i] == '/') return &s[i+1];
        i--;
    }
    return &s[i];
}

;


int gettissue(char t[], char puthere[])
{
strcpy(puthere,"Unknown");
#if 0
    int i;
    for (i=0;ovs[i];i++) { if (strncmp(ovs[i],t,12) == 0) { strcpy(puthere,"Ovary"); return 1; } }
    for (i=0;gbms[i];i++) { if (strncmp(ovs[i],t,12) == 0) { strcpy(puthere,"Brain"); return 1; } }
    strcpy(puthere,"other");
#endif
    return 0;
}


#if WEB_SERVER
void do_describe_based_on_filename(char fn[])
{
    int fzidx = -1;
    char *z;
    char msg[500];
    char m[500];

jdebug("in do_describe_based_on_filename start"); 
jdebug(fn); 
    if (filez_id < 0) 
    {
sprintf(m,"ERROR in do_describe_based_on_filename , invalid filez_id"); jdebug(m);
        return;
    }

    fzidx = get_filez_index_for_id(filez_id);
    if (fzidx < 0) 
    {
sprintf(m,"ERROR in do_describe_based_on_filename , invalid fzidx"); jdebug(m);
        return;
    }
// sprintf(m,"in do_describe_based_on_filename, filez_id=%d",filez_id);  jdebug(m); 
// sprintf(m,"%s",fn);  jdebug(m); 

    printf("<small>");
    printf("%s",filez[fzidx].fullpath);
    if (strstr(filez[fzidx].fullpath,"TCGA"))
    {
        z = getbasename(fn);
        if (z) 
        {
            if (strstr(z,"TCGA"))
            {
                if (isT(z)) printf(" <font color=\"red\" >TUMOR</font>"); 
                if (isN(z)) printf(" <font color=\"blue\" >NORMAL</font>"); 
                if (gettissue(z,msg)) printf(" %s ",msg); 
            }
        }
jdebug("in do_describe_based_on_filename 4"); 
    }
    else if (strstr(filez[fzidx].fullpath,"TARGET"))
    {
jdebug("in do_describe_based_on_filename before getbasename"); 
        z = getbasename(fn);
        if (z) 
        {
            if (strstr(z,"TARGET"))
            {
sprintf(m,"rpf in do_describe_based_on_filename(), fn=%s",fn);  jdebug(m); 
sprintf(m,"rpf in do_describe_based_on_filename(), z=%s",z);  jdebug(m); 
            if (strncmp(z+17,"01",2) == 0) printf(" <font color=\"red\" >TUMOR</font>"); 
            else if (strncmp(z+17,"02",2) == 0) printf(" <font color=\"red\" >TUMOR METASTATIC</font>"); 
            else if (strncmp(z+17,"10",2) == 0) printf(" <font color=\"blue\" >NORMAL BLOOD</font>"); 
            else if (strncmp(z+17,"11",2) == 0) printf(" <font color=\"blue\" >NORMAL TISSUE</font>"); 
            if (gettissue(z,msg)) printf(" %s ",msg); 
            }
        }
    }

    printf("</small>");
    fflush(stdout); 
sprintf(m,"end do_describe_based_on_filename(), filez_id=%d",filez_id);  jdebug(m); 
    return;
}
#endif




void badbuild()
{
    printf("Content-type: text/html\n\n");
    did_content_type = 1;
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>CBIIT ALVIEW</title>\n"); 
    print_some_css();
    printf("</head>\n"); 
    printf("<body style=\"margin:10px;padding:0px;\">\n"); 
    printf("Can't determine build, file may not exist.<br>\n"); 
    printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview\">For ALVIEW Main Page, click here.</a></font><br>\n"); 
    printf("</body>");
    fflush(stdout);
    exit(0);
}

void badness(char *msg)
{
    printf("Content-type: text/html\n\n");
    did_content_type = 1;
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>CBIIT ALVIEW</title>\n"); 
    print_some_css();
    printf("</head>\n"); 
    printf("<body style=\"margin:10px;padding:0px;\">\n"); 
    printf("Bad situation <br>\n"); 
    printf("%s<br>\n",msg); 
    printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview\">For ALVIEW Main Page, click here.</a></font><br>\n"); 
    printf("</body>");
    fflush(stdout);
    exit(0);
}



struct wigtype
{
    int SeqCoverageID   ;
    char *SuperProjectName;
    char *ProjectName;
    char *BAMFileName;
    char *BAMFilePath;
    char *COVFileName;
    char *COVFilePath;
    char *BWFileName;
    char *BWFilePath;
    char *ProjectSampleName;
    char *DiseaseStatus;
    char *Tissue;
    char *genome_build;
    int ViewPriority;
    int SuperProjectID;
    int ProjectID;
    int SampleID;
    int sort_order;
    int TrackType;
    char *labelbefore;
};

#define MAXWIGS 14000
struct wigtype wigs[MAXWIGS];
int numwigs = 0;

// static int loaded_wigs_flag = 0;


int cmp_wigs(const void *arg1, const void *arg2)
{
    struct wigtype *z1;
    struct wigtype *z2;

    z1 = (struct wigtype *)arg1;
    z2 = (struct wigtype *)arg2;

    if (z1->ProjectID > z2->ProjectID) return 1;
    else if (z1->ProjectID < z2->ProjectID) return -1;

    if (z1->sort_order > z2->sort_order) return 1;
    else if (z1->sort_order < z2->sort_order) return -1;

    return 0;
}


void set_samples_string(char psn[], char ss[])  // args=wigs[i].ProjectSampleName,samples_string
{
   char t[512]; // temp string
   int i;

   ss[0] = (char)0;
   for (i=0;i<numwigs;i++)
   {
       if (wigs[i].ProjectSampleName)
       {
           if (strcmp(psn,wigs[i].ProjectSampleName) == 0)
           {
              sprintf(t,"%s",wigs[i].ProjectSampleName);
              if (ss[0]) strcat(ss,",");
              strcat(ss,t);
           }
       }
   }

}


#if USE_MYSQL
void load_wigs(void)     // load (misnamed?) wig "wiggle" tracks 
{
    int len;
    char *z;
    int rowcount;
    MYSQL_ROW row;
    MYSQL *conn = (MYSQL *)0;
    MYSQL_RES *result = (MYSQL_RES *)0;
    char q[4096];
    char m[2048];


#if 1
    sprintf(q,"select distinct sq.SeqCoverageID,sp.SuperProjectName,p.ProjectName,sq.BAMFileName,sq.BAMFilePath,sm.SampleName,sm.DiseaseStatus,sp.SuperProjectID,sq.ProjectID,sm.SampleID ,sq.CoverageFileName,sq.CoverageFilePath,sq.BWFileName,sq.BWFilePath,sm.Tissue, sq.ViewPriority , p.TrackType ,sq.genome_build from cgwb.Project as p , cgwb.SuperProject as sp, hg18.SeqCoverage as sq, hg18.Sample as sm where (sq.obsolete=0) and (sq.ViewPriority<=1) and sm.SampleID=sq.SampleID and sm.SampleID=sq.SampleID and p.ProjectID=sq.ProjectID and sp.SuperProjectID=p.SuperProjectID  order by SuperProjectName,ProjectName,SampleName,DiseaseStatus");


#else // test 
    sprintf(q,"select distinct sq.SeqCoverageID,sp.SuperProjectName,p.ProjectName,sq.BAMFileName,sq.BAMFilePath,sm.SampleName,sm.DiseaseStatus,sp.SuperProjectID,sq.ProjectID,sm.SampleID ,sq.CoverageFileName,sq.CoverageFilePath,sq.BWFileName,sq.BWFilePath,sm.Tissue, sq.ViewPriority , p.TrackType ,sq.genome_build from cgwb.Project as p , cgwb.SuperProject as sp, hg18.SeqCoverage as sq, hg18.Sample as sm where (sq.obsolete=0) and (sq.ViewPriority<=1) and sm.SampleID=sq.SampleID and sm.SampleID=sq.SampleID and p.ProjectID=sq.ProjectID and sp.SuperProjectID=p.SuperProjectID and sq.genome_build like '%s%%' order by SuperProjectName,ProjectName,SampleName,DiseaseStatus",DB);
#endif


    conn = mysql_init(conn);
    if (conn == (MYSQL *)0)
    {
        sprintf(m,"ERROR: CAN NOT sql init"); jdebug(m);
        printf("ERROR: CAN NOT sql init");
        return;
    }

    mysql_options(conn,MYSQL_OPT_COMPRESS,0);
    mysql_options(conn,MYSQL_READ_DEFAULT_GROUP,"odbc");
    if (!mysql_real_connect(conn,"localhost","root","ncigenome","hg18",0,NULL,0))
    {
        sprintf(m, "Failed to connect to database - load_wigs: Error: %s\n", mysql_error(conn)); jdebug(m); mysql_close(conn);
        return;
    }

    if (mysql_real_query(conn, q, strlen(q)) != 0)
    {
        sprintf(m, "ERROR: Failed to query database - load_wigs: Error: %s\n", mysql_error(conn)); jdebug(m); mysql_close(conn);
        return;
    }

    result = mysql_store_result(conn);
    if (result == (void *)0)
    {
        sprintf(m, "ERROR: no results \n"); jdebug(m);
        printf( "ERROR: no results <br>\n");
        mysql_close(conn);
        return;
    }

    numwigs = rowcount = 0;
    while ((row = mysql_fetch_row(result)) != (void *)0)
    {
        if (numwigs == MAXWIGS )
        {  
           mysql_close(conn); 
           sprintf(m,"BIG ERROR: overflow in load_wigs() %d",numwigs);  jdebug(m); 
           printf("pound define error code:4410 "); 
           return; 
        }
        wigs[numwigs].SeqCoverageID = atoi(row[0]);

        if (row[1]) 
        {
            len = strlen(row[1]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4411"); return; } 
            memset(z,0,len);
            strcpy(z,row[1]);
            wigs[numwigs].SuperProjectName = z;
        }
        else wigs[numwigs].SuperProjectName  = (char *)0;

        if (row[2]) 
        {
            len = strlen(row[2]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4412"); return; } 
            memset(z,0,len);
            strcpy(z,row[2]);
            wigs[numwigs].ProjectName = z;
        }
        else wigs[numwigs].ProjectName  = (char *)0;

        if (row[3]) 
        {
            len = strlen(row[3]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4413"); return; } 
            memset(z,0,len);
            strcpy(z,row[3]);
            wigs[numwigs].BAMFileName = z;
        }
        else wigs[numwigs].BAMFileName  = (char *)0;

        if (row[4])
        {
            len = strlen(row[4]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4414"); return; } 
            memset(z,0,len);
            strcpy(z,row[4]);
            wigs[numwigs].BAMFilePath = z;
        }
        else wigs[numwigs].BAMFilePath  = (char *)0;

        if (row[5])
        {
            len = strlen(row[5]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4415"); return; } 
            memset(z,0,len);
            strcpy(z,row[5]);
            wigs[numwigs].ProjectSampleName = z;
        }
        else wigs[numwigs].ProjectSampleName   = (char *)0;

        if (row[6])
        {
            len = strlen(row[6]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4416"); return; } 
            memset(z,0,len);
            strcpy(z,row[6]);
            wigs[numwigs].DiseaseStatus = z;
        }
        else wigs[numwigs].DiseaseStatus   = (char *)0;

        if (row[7])wigs[numwigs].SuperProjectID = atoi(row[7]);
        else       wigs[numwigs].SuperProjectID = -1;
        if (row[8])wigs[numwigs].ProjectID = atoi(row[8]);
        else       wigs[numwigs].ProjectID = -1;
        if (row[9])wigs[numwigs].SampleID = atoi(row[9]);
        else       wigs[numwigs].SampleID = -1;

        if (row[10]) 
        {
            len = strlen(row[10]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4415"); return; } 
            memset(z,0,len);
            strcpy(z,row[10]);
            wigs[numwigs].COVFileName = z;
        }
        else wigs[numwigs].COVFileName  = (char *)0;

        if (row[11])
        {
            len = strlen(row[11]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4416"); return; } 
            memset(z,0,len);
            strcpy(z,row[11]);
            wigs[numwigs].COVFilePath = z;
        }
        else wigs[numwigs].COVFilePath  = (char *)0;

        if (row[12]) 
        {
            len = strlen(row[12]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4415"); return; } 
            memset(z,0,len);
            strcpy(z,row[12]);
            wigs[numwigs].BWFileName = z;
        }
        else wigs[numwigs].BWFileName  = (char *)0;

        if (row[13])
        {
            len = strlen(row[13]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4416"); return; } 
            memset(z,0,len);
            strcpy(z,row[13]);
            wigs[numwigs].BWFilePath = z;
        }
        else wigs[numwigs].BWFilePath  = (char *)0;


        if (row[14])
        {
            len = strlen(row[14]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4417"); return; } 
            memset(z,0,len);
            strcpy(z,row[14]);
            wigs[numwigs].Tissue = z;
        }
        else wigs[numwigs].Tissue  = (char *)0;

        if (row[15]) wigs[numwigs].ViewPriority  = atoi(row[15]);
        else wigs[numwigs].ViewPriority  = 99;

        if (row[16]) wigs[numwigs].TrackType = atoi(row[16]);
        else         wigs[numwigs].TrackType = -1;

        if (row[17])
        {
            len = strlen(row[17]) + 2; 
            z = malloc(len);
            if (z == (char *)0) { mysql_close(conn);  printf("alloc error 4419"); return; } 
            memset(z,0,len);
            strcpy(z,row[17]);
            wigs[numwigs].genome_build = z;
        }
        else wigs[numwigs].genome_build  = (char *)0;

        numwigs++;
    }

    mysql_free_result(result);
    mysql_close(conn);

sprintf(m,"end load_wigs , numwigs = %d ",numwigs); jdebug(m);

    return;
}
#endif


#if 0
void dump_wigs(void) // debug routine
{
int k;
FILE *fp;
fp = fopen("/tmp/junk.me","w");
   for (k=0;k<numwigs;k++)
   {
        fprintf(fp,"%d %s\n",k,wigs[k].BAMFileName);
   }
fclose(fp);
    return;
}
#endif

#if WEB_SERVER
void pick_form()
{
    int tmp1;
    int flushtr;

    tmp1 =  0;
    flushtr = 0;

    print_some_javascript_for_popups();
    printf("<br>\n");
    printf("PICK Page<br>\n");
    printf("<form name=\"form1\" action=\"../cgi-bin/%s\">\n",progname);
    printf("Query: ");
    // support carriage return 
    printf("<input type=\"text\" name=\"pickq\" id=\"pickq\" size=\"40\" onKeyPress=\"return picksubmitenter(this.form,event)\"> \n");

    printf("<span class=\"submit\"><INPUT TYPE=BUTTON OnClick=\"SubmitPICK(this.form); return false;\" VALUE=\"%s\"></span>\n","SUBMIT");
    printf("<br>\n");
    printf(" Enter free form query. <br>\n");
    printf(" Token recognized DISEASE SAMPLE rsSNPID affySNPID HUGOgenename position. <br>\n");
    printf("Example : \"OV rs2453 \" - generate links to rs2453 for Ovaraian <br>\n");
    printf("Example : \"TCGA-B6-A0I6-01A TP53\" - generate links to tumor,TP53<br>\n");
    printf("Example : \"TCGA-B6-A0I6 PTEN\" - generate links to tumor and normal for sample,PTEN<br>\n");
    printf("Example : \"GBM rs4934282\" - generate links for ALL GBM to rs4934282<br>\n");
    printf("Example : \"UCEC chr10:50-90\" - generate links for ALL Uterine Cancer to chr10:50-90<br>\n");
    printf("Example : \"OV SNP_A-2257531\" - find affysnposition SNP_A-2257531 (rs1323341) for all ovarians<br>\n");
    printf("Example : \"GBM SNP_A-8316240\" - find affysnposition  for all gbm<br>\n");
    printf("Example : \"GBM SNP_A-8489775\" - find affysnposition  for all gbm<br>\n");
    printf("Example : \"GBM SNP_A-840057\" - find affysnposition   for all gbm<br>\n");
    printf("Example : \"GBM SNP_A-8587253\" - find affysnposition  for all gbm<br>\n");
    printf("Example: \"OV rs4934282\" C10orf116/AGAP11<br>\n");
    printf("Example: \"OV rs1323341\" ZDHHC21 <br>\n");
    printf("Example: \"OV rs1857623\" DNAH14 <br>\n");

    printf("</form>\n");
}
#endif


int isdisease(char *param)
{
    int i;
    for (i=0;diseases[i].ln;i++)
        if (strcmp(diseases[i].sn,param) == 0) return i;
    return -1;
}


#if USE_MYSQL
int ishugo(char *gene)
{
    int count = 1;
    int ret;
    char m[MAXBUFF];
    char q[MAXBUFF];
    MYSQL_ROW row;
    MYSQL_RES *result = (MYSQL_RES *)0;
    MYSQL *conn = (MYSQL *)0;

    ret = 0;
    sprintf(q,"SELECT geneName from refFlat where geneName='%s' or name='%s'",gene,gene);
// printf("chkhugo %s %s <br>\n",gene,q);
/*
Field	Type	Null	Key	Default	Extra
geneName	varchar(255)	NO	MUL		
name	varchar(255)	NO	MUL		
chrom	varchar(255)	NO	MUL		
strand	char(1)	NO			
txStart	int(10) unsigned	NO		0	
txEnd	int(10) unsigned	NO		0	
cdsStart	int(10) unsigned	NO		0	
cdsEnd	int(10) unsigned	NO		0	
exonCount	int(10) unsigned	NO		0	
exonStarts	longblob	NO			
exonEnds	longblob	NO			
*/
    conn = mysql_init(conn);
    if (conn == (MYSQL *)0)
    {
        sprintf(m, "ERROR: CAN NOT sql init");
        jdebug(m);
        return 0;
    }

    mysql_options(conn,MYSQL_OPT_COMPRESS,0);
    mysql_options(conn,MYSQL_READ_DEFAULT_GROUP,"odbc");
    if (!mysql_real_connect(conn,"localhost","root","ncigenome","hg19",0,NULL,0))
    {
        sprintf(m, "Failed to connect to database: Error: %s\n", mysql_error(conn));
        jdebug(m);
        mysql_close(conn);
        return 0;
    }

    if (mysql_real_query(conn, q, strlen(q)) != 0)
    {
        sprintf(m, "ERROR: Failed to query to database: q=\"%s\" mysqlerror: %s\n", q,mysql_error(conn)); jdebug(m);
        return 0;
    }
    result = mysql_store_result(conn);
    if (result == (void *)0)
    {
        sprintf(m, "ERROR: no results \n"); jdebug(m); printf( "ERROR: no results <br>\n");
        return 0;
    }
    count = 1;
    while ((row = mysql_fetch_row(result)) != (void *)0)
    {
        ret = 1;
        count++;
        ret = 1;
    }
    mysql_free_result(result);
    mysql_close(conn);
    return ret;
}

/*
hg19.snpArrayAffy6
Field   Type    Null    Key     Default Extra
bin     int(10) unsigned        NO              0
chrom   varchar(255)    NO      MUL
chromStart      int(10) unsigned        NO              0
chromEnd        int(10) unsigned        NO              0
name    varchar(255)    NO      MUL
score   int(10) unsigned        NO              0
strand  char(1) NO
observed        blob    NO
rsId    varchar(64)     YES             NULL
*/


int getrs_4affysnp(char *snp,char *rs)
{
    int count = 1;
    int ret;
    char m[MAXBUFF];
    char q[MAXBUFF];
    char locDB[MAXBUFF];
    MYSQL_ROW row;
    MYSQL_RES *result = (MYSQL_RES *)0;
    MYSQL *conn = (MYSQL *)0;

/*
chrom   varchar(31)     NO      MUL
chromStart      int(10) unsigned        NO
chromEnd        int(10) unsigned        NO
name    varchar(15)     NO      MUL
*/

    ret = 0;
    sprintf(q,"SELECT rsID from hg19.snpArrayAffy6 where name='%s'",snp);
    conn = mysql_init(conn);
    if (conn == (MYSQL *)0)
    {
        sprintf(m, "ERROR: CAN NOT sql init");
        jdebug(m);
        return 0;
    }

    mysql_options(conn,MYSQL_OPT_COMPRESS,0);
    mysql_options(conn,MYSQL_READ_DEFAULT_GROUP,"odbc");
    if (!mysql_real_connect(conn,"localhost","root","ncigenome",locDB,0,NULL,0))
    {
        sprintf(m, "Failed to connect to database: Error: %s\n", mysql_error(conn));
        jdebug(m);
        mysql_close(conn);
        return 0;
    }

    if (mysql_real_query(conn, q, strlen(q)) != 0)
    {
        sprintf(m, "ERROR: Failed to query to database: q=\"%s\" mysqlerror: %s\n", q,mysql_error(conn)); jdebug(m);
        return 0;
    }
    result = mysql_store_result(conn);
    if (result == (void *)0)
    {
        sprintf(m, "ERROR: no results \n"); jdebug(m); printf( "ERROR: no results <br>\n");
        return 0;
    }
    count = 1;
    while ((row = mysql_fetch_row(result)) != (void *)0)
    {
       strcpy(rs,row[0]);
       break;
    }
    mysql_free_result(result);
    mysql_close(conn);
    return 1;
}

int getrs_withbuild(char *rs, char *bld, char *pos)
{
    int count = 1;
    int ret;
    char m[MAXBUFF];
    char q[MAXBUFF];
    char locDB[MAXBUFF];
    MYSQL_ROW row;
    MYSQL_RES *result = (MYSQL_RES *)0;
    MYSQL *conn = (MYSQL *)0;

/*
chrom   varchar(31)     NO      MUL
chromStart      int(10) unsigned        NO
chromEnd        int(10) unsigned        NO
name    varchar(15)     NO      MUL
*/

    ret = 0;
    if (strcmp("hg18",bld) == 0) 
    {
        sprintf(q,"SELECT chrom,chromStart,chromEnd from hg18.snp129 where name='%s'",rs);
        strcpy(locDB,"hg18");
    }
    else if (strcmp("hg18",bld) == 0) 
    {
        sprintf(q,"SELECT chrom,chromStart,chromEnd from hg19.snp132 where name='%s'",rs);
        strcpy(locDB,"hg19");
    }
    else if (strcmp("hg38",bld) == 0) 
    {
        sprintf(q,"SELECT chrom,chromStart,chromEnd from hg38.snp132 where name='%s'",rs);
        strcpy(locDB,"hg19");
    }
    conn = mysql_init(conn);
    if (conn == (MYSQL *)0)
    {
        sprintf(m, "ERROR: CAN NOT sql init");
        jdebug(m);
        return 0;
    }

    mysql_options(conn,MYSQL_OPT_COMPRESS,0);
    mysql_options(conn,MYSQL_READ_DEFAULT_GROUP,"odbc");
    if (!mysql_real_connect(conn,"localhost","root","ncigenome",locDB,0,NULL,0))
    {
        sprintf(m, "Failed to connect to database: Error: %s\n", mysql_error(conn));
        jdebug(m);
        mysql_close(conn);
        return 0;
    }

    if (mysql_real_query(conn, q, strlen(q)) != 0)
    {
        sprintf(m, "ERROR: Failed to query to database: q=\"%s\" mysqlerror: %s\n", q,mysql_error(conn)); jdebug(m);
        return 0;
    }
    result = mysql_store_result(conn);
    if (result == (void *)0)
    {
        sprintf(m, "ERROR: no results \n"); jdebug(m); printf( "ERROR: no results <br>\n");
        return 0;
    }
    count = 1;
    while ((row = mysql_fetch_row(result)) != (void *)0)
    {
       sprintf(pos,"%s:%d-%d",row[0],atoi(row[1])-100,atoi(row[2]) +100);
    }
    mysql_free_result(result);
    mysql_close(conn);
    return ret;
}

int check_4_bld_again(char bamfn[])
{
    int j,lastslash;
    int count = 1;
    int ret;
    char sfn[MAXBUFF];
    char m[MAXBUFF];
    char q[MAXBUFF];
    char locDB[MAXBUFF];
    MYSQL_ROW row;
    MYSQL_RES *result = (MYSQL_RES *)0;
    MYSQL *conn = (MYSQL *)0;

    for (lastslash=j=0; bamfn[j]; j++)
        if (bamfn[j] == '/') lastslash = j;
   strcpy(sfn,&bamfn[lastslash+1]); 
/*
chrom   varchar(31)     NO      MUL
chromStart      int(10) unsigned        NO
chromEnd        int(10) unsigned        NO
name    varchar(15)     NO      MUL
*/

    ret = 0;
    sprintf(q,"SELECT genome_build from SeqCoverage where BAMFileName='%s'",sfn);
    strcpy(locDB,"hg19");
    conn = mysql_init(conn);
    if (conn == (MYSQL *)0)
    {
        sprintf(m, "ERROR: CAN NOT sql init");
        jdebug(m);
        return 0;
    }

    mysql_options(conn,MYSQL_OPT_COMPRESS,0);
    mysql_options(conn,MYSQL_READ_DEFAULT_GROUP,"odbc");
    if (!mysql_real_connect(conn,"localhost","root","ncigenome",locDB,0,NULL,0))
    {
        sprintf(m, "Failed to connect to database: Error: %s\n", mysql_error(conn));
        jdebug(m);
        mysql_close(conn);
        return 0;
    }

    if (mysql_real_query(conn, q, strlen(q)) != 0)
    {
        sprintf(m, "ERROR: Failed to query to database: q=\"%s\" mysqlerror: %s\n", q,mysql_error(conn)); jdebug(m);
        return 0;
    }
    result = mysql_store_result(conn);
    if (result == (void *)0)
    {
        sprintf(m, "ERROR: no results \n"); jdebug(m); printf( "ERROR: no results <br>\n");
        return 0;
    }
    count = 1;
    while ((row = mysql_fetch_row(result)) != (void *)0)
    {
       if (strstr(row[0],"hg38")) ret = 38;
       else if (strstr(row[0],"hg19")) ret = 19;
       else ret = 18;
       break;
    }
    mysql_free_result(result);
    mysql_close(conn);
    return ret;
}

void getrs(char *rs, char *pos18,char *pos19) // "rs" is for dbsnp
{
    getrs_withbuild(rs, "hg18",pos18);
    getrs_withbuild(rs, "hg19",pos19);
}


int getfldfromfileznum(int idx)
{
    int j,k,lastslash;
    char s[512];
    char fn[512];

// sprintf(m,"in getfldfromfileznum 1 %d ",idx); jdebug(m); 
    fn[0] = (char)0;
    lastslash = k = 0;
    strcpy(s,filez[idx].fullpath);
    j = strlen(s) -1 ;
    if (j<=0) 
    {
printf("error: in getfldfromfileznum() idx = %d %p\n",idx,filez[idx].fullpath); fflush(stdout);
if (filez[idx].fullpath) printf("error: in getfldfromfileznum() idx = %d s=%s\n",idx,filez[idx].fullpath);
        return 0;
    }

    for (j=0;s[j];j++)
        if (s[j] == '/') lastslash = j;
   strcpy(fn,&s[lastslash+1]); 
   for (k=0;k<numwigs;k++)
   {
// if (k == 735) printf("k735:[%s,%s] \n",wigs[k].BAMFileName,fn);
        if (wigs[k].BAMFileName)
        {
            if (strstr(wigs[k].BAMFileName,fn))
            {
               if (strstr(wigs[k].genome_build,"hg18")) return 18;
               if (strstr(wigs[k].genome_build,"hg19")) return 19;
               if (strstr(wigs[k].genome_build,"hg38")) return 38;
               return 19;
            }
        }
   }
//printf("warn: could not find %s in wigs , checked %d of %d j=%d fn=[%s]\n",s,k,numwigs,j,fn);
    return 0;
}

int do_subsamps(char *frce,char *pos, char *pos18, char *pos19)
{
    int i,bld;
    char s[512];
    char m[512];


// sprintf(m,"in do_subsamps() %p %p %p %p",frce,pos,pos18,pos19);   jdebug(m);
// sprintf(m,"in do_subsamps() %s %s %s %s",frce,pos,pos18,pos19);   jdebug(m);
    if (frce == (void *)0)
    {
sprintf(m,"ERROR: in do_subsamps(), null passed in for frce");   jdebug(m);
        return -1;
    }
    if (pos == (void *)0)
    {
sprintf(m,"ERROR: in do_subsamps(), null passed in for pos");   jdebug(m);
        return -1;
    }
    if (pos18 == (void *)0)
    {
sprintf(m,"ERROR: in do_subsamps(), null passed in for pos18");   jdebug(m);
        return -1;
    }

    for (i=0 ; i < num_filez; i++) 
    {
       if (strstr(filez[i].fullpath,frce))
       {
           if (pos18[0])
           {
               bld = getfldfromfileznum(i);
               if (bld == 0) bld = check_4_bld_again(filez[i].fullpath);
               if (bld == 19) strcpy(s,pos19);
               else if (bld == 18) strcpy(s,pos18);
           }
           else 
           {
               strcpy(s,pos);
           }
       }
    }
    return -1;
}
#endif


void getsn4pj(char *s, char *sn)
{
     if (strstr(s,"STAD")) { strcpy(sn,"STAD"); return; }
     if (strstr(s,"Stomach")) { strcpy(sn,"STAD"); return; }
     if (strstr(s,"THCA")) { strcpy(sn,"THCA"); return; }
     if (strstr(s,"Thyroid")) { strcpy(sn,"THCA"); return; }
     if (strstr(s,"LAML")) { strcpy(sn,"LAML"); return; }
     if (strstr(s,"Luekemia")) { strcpy(sn,"LAML"); return; }
     if (strstr(s,"LGG")) { strcpy(sn,"LGG"); return; }
     if (strstr(s,"Lower Grade")) { strcpy(sn,"LGG"); return; }
     if (strstr(s,"BRCA")) { strcpy(sn,"BRCA"); return; }
     if (strstr(s,"Breast")) { strcpy(sn,"BRCA"); return; }
     if (strstr(s,"LAML")) { strcpy(sn,"LAML"); return; }
     if (strstr(s,"Leukemia")) { strcpy(sn,"LAML"); return; }
     if (strstr(s,"Colom")) { strcpy(sn,"COAD"); return; }
     if (strstr(s,"COAD")) { strcpy(sn,"COAD"); return; }
     if (strstr(s,"GBM")) { strcpy(sn,"GBM"); return; }
     if (strstr(s,"Glioblastoma")) { strcpy(sn,"GBM"); return; }
     if (strstr(s,"clear cell")) { strcpy(sn,"KIRC"); return; }
     if (strstr(s,"KIRC")) { strcpy(sn,"KIRC"); return; }
     if (strstr(s,"KIRP")) { strcpy(sn,"KIRP"); return; }
     if (strstr(s,"papillary")) { strcpy(sn,"KIRP"); return; }
     if (strstr(s,"LIHC")) { strcpy(sn,"KIRP"); return; }
     if (strstr(s,"Liver")) { strcpy(sn,"KIRP"); return; }
     if (strstr(s,"Lung squamous")) { strcpy(sn,"LUSC"); return; }
     if (strstr(s,"LUSC")) { strcpy(sn,"LUSC"); return; }
     if (strstr(s,"Lung adenocarinoma")) { strcpy(sn,"LUAD"); return; }
     if (strstr(s,"LUAD")) { strcpy(sn,"LUAD"); return; }
     if (strstr(s,"Ovarian")) { strcpy(sn,"OV"); return; }
     if (strstr(s,"OVARIAN")) { strcpy(sn,"OV"); return; }
     if (strstr(s,"OV")) { strcpy(sn,"OV"); return; }
/*
    {"Acute Myeloid Leukemia", "LAML" } , 
    {"Bladder Urothelial Carcinoma", "BLCA" } , 
    {"Brain Lower Grade Glioma", "LGG" } , 
    {"Breast invasive carcinoma", "BRCA" } , 
    {"Cervical Squamous Cell Carcinoma", "CESC" } , 
    {"Colon adenocarcinoma", "COAD" } , 
    {"Glioblastoma multiforme", "GBM" } , 
    {"Head and Neck squamous cell carcinoma", "HNSC" } , 
    {"Kidney renal clear cell carcinoma", "KIRC" } , 
    {"Kidney renal papillary cell carcinoma", "KIRP" } , 
    {"Liver hepatocellular carcinoma", "LIHC" } , 
    {"Lung adenocarcinoma", "LUAD" } , 
    {"Lung squamous cell carcinoma", "LUSC" } , 
    {"Ovarian serous cystadenocarcinoma", "OV" } , 
    {"Pancreatic adenocarcinoma", "PAAD" } , 
    {"Prostate adenocarcinoma", "PRAD" } , 
    {"Rectum adenocarcinoma", "READ" } , 
    {"Stomach adenocarcinoma", "STAD" } , 
    {"Thyroid carcinoma", "THCA" } , 
*/
    strcpy(sn,"UNKNONWN");
}

#if USE_MYSQL
void do_pick_q()
{
    // char m[512];
    char pos18[512];
    char pos19[512];
    int i;
    int numtokens;
    char *t = (char *)0;
    char sn[1024];
    char t1[1024];
    char t2[1024];
    char t3[1024];
    char t4[1024];
    char t5[1024];
    char t6[1024];
    char disease[1024];
    char sample[1024];
    char rs[1024];
    char snp[1024];
    char gene[1024];
    char pos[1024];


    pos18[0] = pos19[0] = (char)0;
    disease[0] = sample[0] = rs[0] = snp[0] = gene[0] = pos[0] = (char)0;

    printf("Content-type: text/html\n\n");
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>CBIIT ALVIEW</title>\n"); 
    print_some_css();
    printf("</head>\n"); 
    printf("<body style=\"margin:10px;padding:0px;\">\n"); 
    printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?menu=1\">For ALVIEW Main Page, click here.</a></font>&nbsp;\n"); 
#if PICK
    printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?pickmenu=1\">For ALVIEW PICK Page, click here.</a></font><br>\n"); 
#endif

                  // incoming like this : "OV TCGA-B6-A0I6-01A rs256000 SNP TP53 chr1:1-5"
    numtokens = sscanf(pickq,"%s %s %s %s %s %s",t1,t2,t3,t4,t5,t6);

// printf("pickq = %s , tokens=%d <br>\n",pickq,numtokens);

    for (i=0;i<numtokens;i++)
    {
        switch (i)
        {
            case 0 : t = &t1[0]; break; 
            case 1 : t = &t2[0]; break; 
            case 2 : t = &t3[0]; break; 
            case 3 : t = &t4[0]; break; 
            case 4 : t = &t5[0]; break; 
            case 5 : t = &t6[0]; break; 
            default : break;
        }
        if (t)
        {
            if (isdisease(t) >= 0) {strcpy(disease,t);  }
            else if (strncmp("rs",t,2) == 0) {strcpy(rs,t);  }
            else if (strncmp("chr",t,3) == 0){strcpy(pos,t);   }
            else if (strncmp("SNP",t,3) == 0){strcpy(snp,t);   }
            else if (ishugo(t)) {strcpy(gene,t);   }
            else {strcpy(sample,t); }
        }
    }

    if (disease[0]) printf("disease = %s \n",disease); 
    if (sample[0]) printf("sample = %s \n",sample); 
    if (rs[0]) printf("rs = %s \n",rs); 
    if (gene[0]) printf("gene = %s \n",gene); 
    if (pos[0]) printf("pos = %s \n",pos); 
    if (snp[0]) printf("snp = %s \n",snp); 
    printf("<br>\n");

    if (snp[0]) 
    {
       getrs_4affysnp(snp,rs);
       printf("NOTE: set rsID to %s for affy %s, <br>\n",rs,snp);
    }
    if (gene[0]) 
    {
        strcpy(pos,gene);
        strcpy(pos18,pos);
        strcpy(pos19,pos);
    }
    else if (rs[0]) 
    {
        getrs(rs,pos18,pos19);
printf("RS.dbsnp: %s hg18:%s hg19:%s <br>\n",rs,pos18,pos19);
    }
    if (pos18[0] == (char)0) strcpy(pos18,pos);
    if (pos19[0] == (char)0) strcpy(pos19,pos);

    load_wigs() ;

// dump_wigs();
    int good = 0;

printf("Results:<br>\n"); 
    if (sample[0])
    {
         do_subsamps(sample,pos,pos18,pos19); //may change pos
         fflush(stdout);
    }
    else if (disease[0])
    {
        for (i=0;i<numwigs;i++)
        {
            getsn4pj(wigs[i].ProjectName,sn);
            good = 0;
            if (strcmp(sn,disease) == 0)
            {
                 if (strstr(wigs[i].genome_build,"hg18")) strcpy(pos,pos18);
                 else                                     strcpy(pos,pos19);
                 do_subsamps(wigs[i].ProjectSampleName,pos,pos18,pos19); //may change pos
            }
        }
        fflush(stdout);
    }
    printf("</body>");
}


#endif

#if WEB_SERVER
void do_pick_form()
{
    printf("Content-type: text/html\n\n");
    did_content_type = 1;
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>CBIIT ALVIEW - Pick BAM File</title>\n"); 
    print_some_css();
    printf("</head>\n"); 
    printf("<body style=\"margin:10px;padding:0px;\">\n"); 
    printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?menu=1\">For ALVIEW Main Page, click here.</a></font><br>\n"); 

    pick_form();
    printf("</body>");
}
#endif


void end_javascript(char chrarg[],long int pos1, long int pos2)
{
printf("<script type=\"text/javascript\">\n"); 

printf("var jalv_chr = \"%s\";\n",chrarg); 
printf("var gx1 = %ld; // genmomic start position \n",pos1); 
printf("var gx2 = %ld; // genomic end position \n",pos2); 

printf("function imgsel_end_func(img,selection)\n"); 
printf("{\n"); 
printf("var startw = 0;\n"); 
printf("var starth = 0;\n"); 
//printf("alert(\"x1,y1,x2,y2 = \"+selection.x1 + \" \" +selection.y1 + \" \" +selection.x2 + \" \" +selection.y2);\n");
printf("    var url = document.URL; \n"); 
printf("    var myimg = document.getElementById('bamimg');\n"); 
printf("    if (myimg)\n"); 
printf("    {\n"); 
printf("        startw = myimg.width;\n"); 
printf("    }\n"); 
printf("      var nx1 = selection.x1;             \n");
printf("      var nx2 = selection.x2;             \n");
printf("      var x1 = gx1 + (nx1 * ((gx2 - gx1) / startw));             \n");
printf("      var x2 = gx1 + (nx2 * ((gx2 - gx1) / startw));             \n");
printf("      x1 = parseInt(x1);              \n");
printf("      x2 = parseInt(x2);              \n");
printf("      if (x1 == x2) { x2 = x1 + 1};;             \n");
printf("    var pos = url.indexOf(\"position=\");  \n"); 
printf("    if (pos == -1) { return false;} \n"); 
printf("    var part1 =   url.slice(0,pos); \n"); 
printf("    var rest =url.substr(pos); \n"); 
printf("    var pos = rest.indexOf(\"&\");  \n"); 
printf("    var part2 = rest.slice(pos); \n"); 
printf("    var position = \"position=\" + jalv_chr + \":\" + x1 + \"-\" + x2; \n");
printf("    u = part1 + position+part2;   \n"); 
printf("    self.location=u;\n"); 
printf("    return false;\n"); 
printf("}\n"); 
printf("function imgsel_sel_func (img,selection)\n"); 
printf("{\n"); 
printf("}\n"); 
printf("function rpf_load_func(img,selection)\n"); 
printf("{\n"); 
printf("    var myimg = document.getElementById('bamimg');\n"); 
printf("    if (myimg)\n"); 
printf("    {\n"); 
printf("        startw = myimg.width;\n"); 
printf("        starth = myimg.height;\n"); 
printf("        // debug alert(\"start : \" + startw + \" \" +starth);\n"); 
printf("    }\n"); 
printf("    else { alert(\"cant find img\"); }\n"); 
printf("}\n"); 
printf("function imgsel_start_func (img,selection)\n"); 
printf("{\n"); 
printf("}\n"); 
printf("$(document).ready(function () {\n"); 
printf("    rpf_load_func();\n"); 
printf("    $('img#bamimg').imgAreaSelect({\n"); 
printf("        handles: true,\n"); 
printf("        minHeight: starth,\n"); 
printf("        onSelectEnd: imgsel_end_func ,\n"); 
printf("        onSelectStart :  imgsel_start_func ,\n"); 
printf("        onSelectChange : imgsel_sel_func ,\n"); 
printf("    });\n"); 
printf("});\n"); 
printf("</script>\n"); 
}


void badhost()
{
    printf("Content-type: text/html\n\n");
    did_content_type = 1;
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>ERROR</title>\n"); 
    printf("</head>\n"); 
    printf("Only available inside NIH\n"); 
    printf("<body>\n"); 
    printf("</body>\n"); 
}


void initialize_ttf(void)
{
    return; 
#if 0
// this is OLD - it was used with gd library. gd is obsoleted to make more portable 

    if (FULL_PATH_TO_TTF[0]) return;     // already set, possibly by config file

#ifdef WIN32
    strcpy(FULL_PATH_TO_TTF,ced);
    strcat(FULL_PATH_TO_TTF,"\\GENOMEDATA\\luximb.ttf");
#else
    strcpy(FULL_PATH_TO_TTF,"/usr/X11R6/lib/X11/fonts/TTF/luximb.ttf");
#endif
    return;
#endif
}

#if WEB_SERVER
int web_server_main(int argc,char *argv[])
{
    int stati;
    int s,e;
    int status;
    int image_type = 0;
    char tcga[MAXBUFF];
    char short_image_filename[MAXBUFF];
    char long_image_filename[MAXBUFF];
    char m[MAXBUFF];
    char tmps[MAXBUFF];
    char mmspot_s[MAXBUFF];
    char exp_s[MAXBUFF];
    char nso_s[MAXBUFF];
    char spliceonly_s[MAXBUFF];
    char filenumstring[MAXBUFF];
    char tds_s[MAXBUFF];
    char tmp_chr[512];
    char url[512];
    char hostname[512];
    hostname[0] = (char)0; // you can hardwire things based on hostname
//     char altonly_s[MAXBUFF];


     strcpy(blds,"hg19"); // default 
     status =  gethostname(hostname,500);
     if (strcmp("lpgws511.nci.nih.gov",hostname))
     {
// checking for what machine thisis running on . "lpgws511" is a test server, so I can disable any secuirty   
// I hack this to see private BAM files  ***** you can do your own hacks on your server
     }

     initialize_ttf();

     tds_val = 0;

     tcga[0] = (char)0;
     spliceonly_s[0] = (char)0;
//     altonly_s[0] = (char)0;

    mmspot_s[0] = exp_s[0] = nso_s[0] = tds_s[0] = (char)0;
sprintf(m," ");  jdebug(m); 
sprintf(m,"-------------- in web_server_main alview -------------");  jdebug(m); 
// sprintf(m,"SIZE_REFFLAT_REC = %d ",SIZE_REFFLAT_REC);  jdebug(m);
    short_image_filename[0] = (char)0;
    srand(time(NULL) + getpid()); 
    alview_load_config();

sprintf(m,"before atypical_go_front_end");  jdebug(m); 
    atypical_go_front_end();                //    printf("Content-type: text/html\n\n");
sprintf(m,"after atypical_go_front_end, before parse");  jdebug(m); 
    un_escape((unsigned char *)buff); 
    parse_cgi_contents(buff);
sprintf(m,"after parse_cgi_contents");  jdebug(m); 

    strcpy(position,"chr17:7520037-7531588");
    get_cgi_param((char *)"position",position,0);
sprintf(m,"web_server_main position=[%s] ",position);  jdebug(m); 
    get_cgi_param((char *)"op",op,0);
    get_cgi_param((char *)"tcga",tcga,0);

  
    menu = -1;
    get_cgi_param((char *)"menu",menustring,0);
    if (menustring[0]) menu = atoi(menustring);

    pickmenu = -1;
    get_cgi_param((char *)"pickmenu",pickmenustring,0);
    if (pickmenustring[0]) pickmenu = atoi(pickmenustring);

    pickq[0] = (char)0; // pick query
    get_cgi_param((char *)"pickq",pickq,0);

    bamname[0] = (char)0; // find by bamname
    get_cgi_param((char *)"bamname",bamname,0);

    filez_id = -1;
    get_cgi_param((char *)"filenum",filenumstring,0);
    if (filenumstring[0]) filez_id = atoi(filenumstring);
sprintf(m,"after get_cgi_param() filenum  ,filez_id = %d [%s]" ,filez_id,filenumstring); jdebug(m); 

    filecgiarg[0] = (char)0;
    get_cgi_param((char *)"file",filecgiarg,0);
    if (filecgiarg[0])
    {
        tds_val = filez_id = get_full_path_for_short_file_name(filecgiarg);
    }

    jutlfn[0] = (char)0;
    get_cgi_param((char *)"jutlfn",jutlfn,0);
    if ( jutlfn[0] )
    {
        tds_val = filez_id = get_tds_for_lfn(jutlfn);
    }

    forcefilecgiarg[0] = (char)0;
    get_cgi_param((char *)"forcefile",forcefilecgiarg,0);

    get_cgi_param((char *)"mmspot",mmspot_s,0);
    if (strcmp(mmspot_s,"on")  == 0 )
        mmspot_flag = 1; 
    else
        mmspot_flag = 0; 

    get_cgi_param((char *)"exp",exp_s,0);
    if (strcmp(exp_s,"on")  == 0 )
        extrapolate_flag = 1; 

    get_cgi_param((char *)"nso",nso_s,0);
    if (strcmp(nso_s,"on")  == 0 )
        nso_flag = 1; 

    get_cgi_param((char *)"spliceonly",spliceonly_s,0);
    if (strcmp(spliceonly_s,"on") == 0 )
        spliceonly_flag = 1; 

#if 0
    get_cgi_param((char *)"altonly",altonly_s,0);
    if (strcmp(altonly_s,"on") == 0 ) altonly_flag = 1; 
#endif 

    tmps[0] = (char)0;
    get_cgi_param((char *)"basecolors",tmps,0);
    if (strcmp(tmps,"on") == 0 )
        basecolors_flag = 1; 

    tmps[0] = (char)0;
    get_cgi_param((char *)"uniq",tmps,0);
    if (strcmp(tmps,"on") == 0 )
        uniq_flag = 1; 

    tmps[0] = (char)0;
    get_cgi_param((char *)"qual",tmps,0);
    if (strcmp(tmps,"on") == 0 )
        qual_flag = 1; 

    jobstring[0] = (char)0;
    get_cgi_param((char *)"job",jobstring,0);

    tmps[0] = (char)0;
    get_cgi_param((char *)"iw",tmps,0);
    iw = atoi(tmps);
    if (iw < 50) iw = 1000; 
    if (iw > 5000) iw = 5000; 
    tmps[0] = (char)0;
    get_cgi_param((char *)"ih",tmps,0);
jdebug("after getting params "); 

    ih = atoi(tmps);
    if (ih < 50) ih = 400;  
    if (ih > 10000) ih = 10000;

    ih = ih;
    iw = iw;

    load_filez();


    if (bamname[0])
    {      // HACK !!!! , this fixes calls from microsoft spreadheet "Excel" url to Mozilla/Firefox browser
       char *tmptr;
       tmptr = strstr(bamname,"bamFirefoxHTML");  //  \Shell\Open\CommandFirefoxHTML\Shell\Open\Command
       if (tmptr) *(tmptr+3) = (char)0 ;
       filez_id = find_fileid_for_bamname(bamname); 
    }
     
#if 0
    if (tcga[0] != (char)0) filez_id = findtcga(tcga);
    else if (forcefilecgiarg[0] != (char)0)
#endif

    if (forcefilecgiarg[0] != (char)0)
    {
         fix_doubleslash_to_singleslash(forcefilecgiarg); 
         filez_id = findforcefile(forcefilecgiarg);  // returns id, not index
         if (filez_id >= 0) forcefilecgiarg[0] = (char)0;
    }

sprintf(m,"after load_filez() , filez = %p , num_filez=%d",(void *)filez,num_filez); jdebug(m);
// sprintf(m,"filez_id = %d ",filez_id);  jdebug(m); 
    fn_bam[0] = (char)0;
    if (forcefilecgiarg[0]) strcpy(fn_bam,filecgiarg); 
    else 
    {
        if (filez_id >= 0)
        {
            int idx;
            idx = get_filez_index_for_id(filez_id);
            if (idx < 0) 
            {
                sprintf(m,"invalid bam file index = %d, id = %d",idx,filez_id); 
                badness(m); // prints message and dies
            }
            if (filez[idx].fullpath) 
            {
                strcpy(fn_bam,filez[idx].fullpath);
sprintf(m,"set fn_bam to [%s] ",fn_bam); jdebug(m);
            }
        }
    }

// sprintf(m,"before find_build %s, filez_id = %d",fn_bam,filez_id);  jdebug(m); 
    stati = find_build(fn_bam); // 18=hg18, 19=hg19, else -1 is error , sets global string variable "blds"
// sprintf(m,"after build = %s, stati=%d",blds,stati);  jdebug(m); 
    if (status == -1) badbuild();
    if  (pickq[0] != (char)0) 
    {
#if USE_MYSQL
         do_pick_q();
#endif
         fflush(stdout);
         exit(0);
    }
    if  (pickmenu == 1) 
    {
         do_pick_form();
         fflush(stdout);
         exit(0);
     }
    if ( ((menu == 1) || (filez_id < 0)) && (filecgiarg[0] == (char)0) )
    {
sprintf(m,"in web_server_main, pumping out headers ");  jdebug(m); 
         printf("Content-type: text/html\n\n");
         did_content_type = 2;
         printf("<head>\n"); 
         printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
         printf("<title>CBIIT ALVIEW</title>\n"); 
         print_some_css();
         printf("</head>\n"); 
         printf("<body  bgcolor=\"#%s\" style=\"margin:10px;padding:0px;\">\n",FAVEBACKGROUNDCOLOR); 
         printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview\">For ALVIEW Main Menu, click here.</a></font>&nbsp;&nbsp;\n"); 
#if PICK
         printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?pickmenu=1\">For ALVIEW Pick Page, click here.</a></font><br>\n"); 
#endif

sprintf(m,"in web_server_main, before do_filez_form() ");  jdebug(m); 
         do_filez_form();
sprintf(m,"in web_server_main, after do_filez_form() ");  jdebug(m); 

         printf("</body>");
         fflush(stdout);
         exit(0);
     }
 
sprintf(m,"before postion, filez = %p ",(void *)filez); jdebug(m);
    if (position[0])
    {
        int ki;

        sprintf(long_image_filename,"images/tmpalview.%d.%d.png",rand(),getpid());
        global_khr[0] = (char)0;

        if (strncmp(position,"chr",3) == 0) 
        {
            parse_position(position,chr,&alv_start,&alv_end); // eg.: "position=chrX:37301314-37347604"
        }
        else
        {
// sprintf(m,"before  do_by_gene_name_from_refflat");  jdebug(m); 
            do_by_gene_name_from_refflat(position,chr,&alv_start,&alv_end); // eg.: "position=chrX:37301314-37347604"
// sprintf(m,"after do_by_gene_name_from_refflat %d %d",alv_start,alv_end);  jdebug(m); 
        }
        strcpy(global_khr,chr); 

        if (alv_end <= 0) 
        {
            alv_end = alv_start + 2;
        }
        else if (alv_start > alv_end)
        {
            unsigned int tmp;
            tmp = alv_end;
            alv_end = alv_start;
            alv_start = tmp;
        }
        s = alv_start;
        e = alv_end;
        ge = e;
        gs = s;
        gdiff = e - s;

        status = get_chr_lo_hi(chr,&s,&e,&ki);  // dont care about ki

        if ((strcmp(jobstring,"fasta") == 0) || (strcmp(jobstring,"fastq") == 0)) 
        {
            if (strcmp(jobstring,"fasta") == 0)
                 global_fastafastq_option =  0;
            else
                 global_fastafastq_option =  1;
#define TEXT 1
#if TEXT
            // printf("Content-type: text/plain; charset=us-ascii\n\n"); 
            printf("Content-type: text/plain; \n\n"); 
            did_content_type = 3;
#else
            printf("Content-type: text/html\n\n");
            did_content_type = 4;
            printf("<head>\n"); 
            printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
            printf("<title>CBIIT ALVIEW</title>\n"); 
            print_some_css();
            printf("</head>\n"); 
            printf("<body style=\"margin:10px;padding:0px;\">\n"); 

            printf("<pre>"); 
#endif

            (void)do_fasta(alv_start,alv_end,image_type,short_image_filename,chr,alv_start,alv_end) ;

#if TEXT
#else
            printf("</pre>"); 
            printf("</body>");
#endif
            return 0;
        }
        else if ( (strcmp(jobstring,"sam") == 0) || (strcmp(jobstring,"sam0") == 0) ||
                  (strcmp(jobstring,"sam1") == 0) || (strcmp(jobstring,"sam2") == 0)  )
        {
            printf("Content-type: text/html\n\n");
            did_content_type = 5;
            printf("<head>\n"); 
            printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
            printf("<title>CBIIT ALVIEW</title>\n"); 
            print_some_css();
            printf("</head>\n"); 
            printf("<body style=\"margin:10px;padding:0px;\">\n"); 

            samoption = 0; 
            if ( (strcmp(jobstring,"sam") == 0) || (strcmp(jobstring,"sam0") == 0)) samoption = 0; 
            else if (strcmp(jobstring,"sam1") == 0) samoption = 1;
            else if (strcmp(jobstring,"sam2") == 0) samoption = 2;
            printf("<pre>"); 
/* BIG WORK HERE ... */
jdebug("before do_sam"); 
            (void)do_sam(alv_start,alv_end,image_type,short_image_filename,chr,alv_start,alv_end) ;
            printf("</pre>"); 
            printf("</body>");
            return 0;
        }
#if GOTOH
        else if (strcmp(jobstring,"align") == 0) 
        {
            printf("Content-type: text/html\n\n");
            did_content_type = 6;
            printf("<head>\n"); 
            printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
            printf("<title>CBIIT ALVIEW</title>\n"); 
            print_some_css();
            printf("</head>\n"); 
            printf("<body style=\"margin:10px;padding:0px;\">\n"); 

            printf("<pre>"); 
(void)do_gotoh(alv_start,alv_end,image_type,short_image_filename,chr,alv_start,alv_end) ;
            printf("</pre>"); 
            printf("</body>");
            return 0;
        }
#endif
        else
        {

////////////
// THE BIG WORK GETS DONE IN imgen()
////////////
            status = imgen(short_image_filename,chr,alv_start,alv_end,image_type) ;
sprintf(m,"after imgen() in webserver web_server_main status=%d",status); jdebug(m); 

            if ((status) || (global_bamerr > 0) )
            {
                sprintf(m,"Error- webserver probably can not generate image. imgen failed %s, status=%d , global_bamerr=%d ",short_image_filename,status,global_bamerr); jdebug(m); 
                printf("Content-type: text/html\n\n");
                did_content_type = 7;
                printf("<head>\n"); 
                printf("</head>\n"); 
                printf("<body>\n"); 
                if (global_bamerr > 0) 
                    printf(m,"<font color=\"red\">ERROR: failed to open bam file<br>\n");

                printf("ERROR after imgen (), possible reasons: bad genomic address. try chrX:5-999 for example OR source bam does not exist.<br>\n");
                printf("status code = %d , global_bamerr code = %d \n",status,global_bamerr);
                printf("Press back button on browser.<br>\n");
                printf("</body>\n"); 
                return 0;
            }
        }
    } // position[0]
sprintf(m,"after position, filez = %p",(void *)filez); jdebug(m);

sprintf(m,"in webserver web_server_main, before contents , status=%d contents=%d",status,did_content_type); jdebug(m); 

    printf("Content-type: text/html\n\n");
    did_content_type = 8;
    printf("<!DOCTYPE html>\n");
    printf("<head>\n"); 
    printf("<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">\n"); 
    printf("<title>CBIIT ALVIEW</title>\n"); 
    print_some_css();


printf("<link rel=\"stylesheet\" type=\"text/css\" href=\"../imgareaselect-default.css\" />\n"); 
printf("<script type=\"text/javascript\" src=\"../js/alview.jquery-1.8.3.js\"></script>\n"); 
printf("<script type=\"text/javascript\" src=\"../js/alview.jquery-ui-1.9.2.custom.js\"></script>\n"); 
printf("<script type=\"text/javascript\" src=\"../js/alview.jquery.imgareaselect.js\"></script>\n"); 
    printf("</head>\n"); 
    printf("<body style=\"margin:10px;padding:0px;  background-color:#%s\">\n",FAVEBACKGROUNDCOLOR); 

// sprintf(m,"filez_id = %d ",filez_id);  jdebug(m); 
// sprintf(m,"filez = %p ",filez);  jdebug(m); 
// sprintf(m,"fullpath = %p ",filez[filez_id].fullpath);  jdebug(m); 
// sprintf(m,"fullpath = %s ",filez[filez_id].fullpath);  jdebug(m); 

    do_describe_based_on_filename(filez[filez_id].fullpath);

    printf("<br>"); 
    printf("<img id=\"bamimg\" WIDTH=%d HEIGHT=%d src=\"../%s\">\n",iw,ih,short_image_filename);

    print_some_javascript_for_popups();
sprintf(m,"after print_some_javascript_for_popups");  jdebug(m); 

    do_form(chr,alv_start,alv_end,s,e);
// sprintf(m,"after do_form");  jdebug(m); 

printf("Using Dataset:%d Position : %s %d %d <br>\n",tds_val,chr,alv_start,alv_end);

    strcpy(tmp_chr,chr);
    if (strcmp(chr,"chr23") == 0) strcpy(tmp_chr,"chrX"); 
    else if (strcmp(chr,"chr24") == 0) strcpy(tmp_chr,"chrY"); 
    else if (strcmp(chr,"chr25") == 0) strcpy(tmp_chr,"chrM"); 

printf("<a href=\"../cgi-bin/fetchdna?position=%s:%d-%d\" target=\"_blank\">Get DNA</a> | \n",tmp_chr,alv_start,alv_end);
printf("Range is %s %d %d (%d base pairs)\n",tmp_chr,alv_start,alv_end,alv_end-alv_start);
printf("<br>\n");

printf("<a href=\"../cgi-bin/hgTracks?db=%s&position=%s:%d-%d\" target=\"_blank\">CGWB Link</a> | \n",blds,tmp_chr,alv_start,alv_end);
printf("<a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=%s&position=%s:%d-%d\" target=\"_blank\">UCSC Link</a> | \n",blds,tmp_chr,alv_start,alv_end);

sprintf(url,"../cgi-bin/bambino?local_file=%s;chr=%s;start=%d;end=%d",fn_bam,tmp_chr,alv_start,alv_end);
printf("<a href=\"%s\" target=\"_blank\">bambino </a> | \n",url);

printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=sam0\" target=\"_blank\">sam</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=sam1\" target=\"_blank\">sam1</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=sam2\" target=\"_blank\">sam2</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=fasta\" target=\"_blank\">fasta</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=fastq\" target=\"_blank\">FQ</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&job=align\" target=\"_blank\">align</a> | \n",tmp_chr,alv_start,alv_end,filez_id);
printf("<a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat\" target=\"_blank\">blat</a> | \n");

 printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&geneannot=on\" >P</a> | \n",tmp_chr,alv_start,alv_end,filez_id-1);
 printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&geneannot=on\" >N</a> | \n",tmp_chr,alv_start,alv_end,filez_id+1);
 printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&geneannot=on\" target=\"_blank\">Pl</a> | \n",tmp_chr,alv_start,alv_end,filez_id-1);
 printf("<a href=\"../cgi-bin/alview?position=%s:%d-%d&filenum=%d&geneannot=on\" target=\"_blank\">Nl</a> | \n",tmp_chr,alv_start,alv_end,filez_id+1);

printf("\n"); 
        printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?menu=1\">For ALVIEW Main Page, click here.</a></font>&nbsp;\n"); 
#if PICK
        printf("<font size=\"3\" color=\"red\"><a href=\"../cgi-bin/alview?pickmenu=1\">For ALVIEW PICK Page, click here.</a></font><br>\n"); 
#endif
printf("<br>\n");
fflush(stdout);  // why not ?

sprintf(m,"before do_gotgenes");  jdebug(m); 
    do_gotgenes();
sprintf(m,"after  do_gotgenes");  jdebug(m); 
    do_gotgenes_trawler();

    end_javascript(chr,alv_start,alv_end);
sprintf(m,"after  end_javascript");  jdebug(m); 

printf("<hr> \n");
printf("<center> \n");
printf("<font size=\"1\"> \n");
printf("CBIIT ALVIEW by Richard Finney : National Cancer Institute - National Institutes of Health - Bethesda, Maryland.\n"); 
printf("</font> \n");
printf("</center>\n");
printf("<hr> \n");
printf("</body>\n"); 
printf("</html>\n"); 
fflush(stdout); 

sprintf(m,"before zaps ");  jdebug(m); 
    zap_tmpalview_files();                                        // *************************************************
sprintf(m,"before freedom_for_memory ");  jdebug(m); 
    freedom_for_memory();
    
sprintf(m,"end web program");  jdebug(m); 
    return 0;
}
#endif


#if CMD_LINE

int fast_snp_call( char fnbam[], char chr[],int khrstart, int khrend)
{
    int s,e;
    int len;
    int chrom_index;
    int status;
    char chrchr[1024];  
    char fn[1024];
    char m[2048];


// debug fprintf(stderr,"in fast_snp_call  0 khrstart=%d khrend=%d chr=[%s] - bam=[%s]\n",khrstart,khrend,chr,fnbam);  fflush(stderr);  

    strcpy(chrchr,chr);
    fix_num_to_chrstyle(chr,chrchr); 

    strcpy(fn,fnbam);

    s = e = len = status = 0;
    totalcnt = 0;
    len = khrstart - khrend;
    if (len > MAX_CHROMSIZE)
    {
        sprintf(m,"ERROR: in fast_snp_call() size is too big "); jdebug(m);
        return -1;
    }

    pxwidth = (double)iw /(double)(len);
    gs = khrstart;
    ge = khrend;
    gdiff = ge - gs;
    ppp_firsttime = 1; 

    status = get_chr_spot(chr,&s,&e,&chrom_index); // not really using this right - fix 
    if (status < 0)
    {
        sprintf(m,"ERROR: Can not get location from get_chr_lo_hi for %s %d %d %s ",chr,gs,ge,blds); 
        jdebug(m);
        return -2;
    }
    snp_call_flag = 1; 
    strcpy(snp_call_chr,chr); 
    snp_call_spot = gs; 
    gs = snp_call_spot - 5000 ;
    ge = snp_call_spot + 200 ;
    snp_call_dnaat = snp_call_spot - gs;
    setup_reference_dna(chrchr,gs,ge,0); // 0 in arg4 means don't paint 
    setup_dnacnts_and_dnamms(ge-gs+1);
    dobam(fn_bam,chrom_index,gs,ge,chr);
    return 0;
}


int setup_and_do_fast_snp_caller(char fnbam[])
{
    double d;
    int cov,mx;
    int s,e;
    int status;

// debug fprintf(stderr,"in setup_and_do_fast_snp_caller(%s) \n",fnbam); fflush(stderr); 
    if ( GENOMEDATADIR[0]  == (char)0 )
    {
        fprintf(stderr,"ERROR, must set up GENOMEDATADIR in configure file \"alview.conf\"\n");
        exit(0);
    }
    status = parse_position(position,chr,&s,&e);


    status = fast_snp_call(fnbam,chr,s,e);
    if ((status) || (global_bamerr > 0) )
    {
        fprintf(stderr,"Error- in setup_and_do_fast_snp_caller() after fast_snp_call() status=%d ",status);
        fprintf(stderr,"status code = %d , global_bamerr code = %d \n",status,global_bamerr);
        return -1;
    }

 cov = ( snp_call_A_cnt  + snp_call_C_cnt  + snp_call_G_cnt  + 
                     snp_call_T_cnt  + snp_call_Ins_cnt  + snp_call_Del_cnt);
 d = 0.0;
 mx = 0;
 if (snp_call_A_cnt  > mx) mx = snp_call_A_cnt;
 if (snp_call_C_cnt  > mx) mx = snp_call_C_cnt;
 if (snp_call_G_cnt  > mx) mx = snp_call_G_cnt;
 if (snp_call_T_cnt  > mx) mx = snp_call_T_cnt;
 if (snp_call_Ins_cnt  > mx) mx = snp_call_Ins_cnt;
 if (snp_call_Del_cnt  > mx) mx = snp_call_Del_cnt;
 if (cov) d = (double)mx / (double) cov;

printf("A: %d ", snp_call_A_cnt ); 
printf("C: %d ", snp_call_C_cnt ); 
printf("G: %d ", snp_call_G_cnt ); 
printf("T: %d ", snp_call_T_cnt ); 
printf("Ins: %d ", snp_call_Ins_cnt ); 
printf("Del: %d ", snp_call_Del_cnt ); 
printf("cov: %d ", cov ); 
printf("max: %d ", mx ); 
printf("readcount: %d ", globalreadscount); 
printf("percmac: %10.7f ", d ); 
printf("ref: %c ", snp_call_referece ); 
printf("\n");
    return 0;
}


int command_line_main(int argc,char *argv[])
{
    int s,e;
    int status;
    char outimg[MAXBUFF];
    char m[MAXBUFF];
    char hostname[512];


    hostname[0] = (char)0;
    if (argc < 3) 
    { 
        fprintf(stderr,"ERROR: usage inbam outpng position build [imageheight] [ imagewidth ]\n"); 
        fprintf(stderr,"or usage2 (experimental \"snp caller\"):  CALL inbam position build \n"); 
        exit(0); 
    }
    if (argc == 4) 
    {
        alview_load_config();
        strcpy(fn_bam,argv[1]);
        strcpy(position,argv[2]);
        strcpy(blds,argv[3]);
        setup_and_do_fast_snp_caller(fn_bam);
        return 0;
    }

    strcpy(fn_bam,argv[1]);
    strcpy(outimg,argv[2]);
    strcpy(position,argv[3]);
    strcpy(blds,argv[4]);
    ih = 450; //default 
    iw = 1000; //deafult 
    if (argc >= 6) ih = atoi(argv[5]); 
    if (argc >= 7) iw = atoi(argv[6]); 


#if 0 // don't use GD graphcis anymore
// defaults 
    strcpy(FULL_PATH_TO_TTF,"/usr/X11R6/lib/X11/fonts/TTF/luximb.ttf");
#endif

    tds_val = 0;
    alview_load_config();

    if ( GENOMEDATADIR[0]  == (char)0 )
    {
        fprintf(stderr,"ERROR, must set up GENOMEDATADIR in configure file \"alview.conf\"\n");
        exit(0);
    }
    status = parse_position(position,chr,&s,&e);

fprintf(stderr,"command_line_main() chr=[%s] start=%d end=%d position=[%s] \n",chr,s,e,position); 

    status = imgen(outimg,chr,s,e,0);
    if ((status) || (global_bamerr > 0) )
    {
        fprintf(stderr,"Error- command_line_main() probably can not generate image. imgen failed %s, status=%d , global_bamerr=%d",outimg,status,global_bamerr);
        fprintf(stderr,"ERROR after imgen (), possible reasons: bad genomic address. try chrX:5-999 for example OR source bam does not exist.<br>\n");
        fprintf(stderr,"status code = %d , global_bamerr code = %d \n",status,global_bamerr);
        return 0;
    }
sprintf(m,"command_line_main done CMD_LINE alview");  jdebug(m); 
    return 0;
}
#endif

#if WX_GUI
// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/sizer.h"
#include "wx/gbsizer.h"
#include "wx/statline.h"
#include "wx/notebook.h"
#include "wx/spinctrl.h"
#include "wx/wrapsizer.h"
#include "wx/generic/stattextg.h"

#include "layout.h"

#ifndef wxHAS_IMAGES_IN_RESOURCES
    //rpf #include "../sample.xpm"
    #include "sample.xpm"
#endif


// ----------------------------------------------------------------------------
// MyApp
// ----------------------------------------------------------------------------

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
  if ( !wxApp::OnInit() )
      return false;

  // Create the main frame window
  MyFrame *frame = new MyFrame;

  frame->Show(true);

  return true;
}

// ----------------------------------------------------------------------------
// MyFrame
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
  EVT_MENU(LAYOUT_ABOUT, MyFrame::OnAbout)
  EVT_MENU(LAYOUT_QUIT, MyFrame::OnQuit)

  EVT_MENU(LAYOUT_TEST_PROPORTIONS, MyFrame::TestProportions)
  EVT_MENU(LAYOUT_TEST_SIZER, MyFrame::TestFlexSizers)
  EVT_MENU(LAYOUT_TEST_NB_SIZER, MyFrame::TestNotebookSizers)
  EVT_MENU(LAYOUT_TEST_GB_SIZER, MyFrame::TestGridBagSizer)
  EVT_MENU(LAYOUT_TEST_SET_MINIMAL, MyFrame::TestSetMinimal)
  EVT_MENU(LAYOUT_TEST_NESTED, MyFrame::TestNested)
  EVT_MENU(LAYOUT_TEST_WRAP, MyFrame::TestWrap)
END_EVENT_TABLE()

// Define my frame constructor
MyFrame::MyFrame()
       : wxFrame(NULL, wxID_ANY, wxT("wxWidgets Layout Demo"))
{
    SetIcon(wxICON(sample));

    // Make a menubar
    wxMenu *file_menu = new wxMenu;

    file_menu->Append(LAYOUT_TEST_PROPORTIONS, wxT("&Proportions demo...\tF1"));
    file_menu->Append(LAYOUT_TEST_SIZER, wxT("Test wx&FlexSizer...\tF2"));
    file_menu->Append(LAYOUT_TEST_NB_SIZER, wxT("Test &notebook sizers...\tF3"));
    file_menu->Append(LAYOUT_TEST_GB_SIZER, wxT("Test &gridbag sizer...\tF4"));
    file_menu->Append(LAYOUT_TEST_SET_MINIMAL, wxT("Test Set&ItemMinSize...\tF5"));
    file_menu->Append(LAYOUT_TEST_NESTED, wxT("Test nested sizer in a wxPanel...\tF6"));
    file_menu->Append(LAYOUT_TEST_WRAP, wxT("Test wrap sizers...\tF7"));

    file_menu->AppendSeparator();
    file_menu->Append(LAYOUT_QUIT, wxT("E&xit"), wxT("Quit program"));

    wxMenu *help_menu = new wxMenu;
    help_menu->Append(LAYOUT_ABOUT, wxT("&About"), wxT("About layout demo..."));

    wxMenuBar *menu_bar = new wxMenuBar;

    menu_bar->Append(file_menu, wxT("&File"));
    menu_bar->Append(help_menu, wxT("&Help"));

    // Associate the menu bar with the frame
    SetMenuBar(menu_bar);

#if wxUSE_STATUSBAR
    CreateStatusBar(2);
    SetStatusText(wxT("wxWidgets layout demo"));
#endif // wxUSE_STATUSBAR

    wxPanel* p = new wxPanel(this, wxID_ANY);

    // we want to get a dialog that is stretchable because it
    // has a text ctrl in the middle. at the bottom, we have
    // two buttons which.

    wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

    // 1) top: create wxStaticText with minimum size equal to its default size
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_RIGHT).") ),
        wxSizerFlags().Align(wxALIGN_RIGHT).Border(wxALL & ~wxBOTTOM, 5));
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_LEFT).") ),
        wxSizerFlags().Align(wxALIGN_LEFT).Border(wxALL & ~wxBOTTOM, 5));
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_CENTRE_HORIZONTAL).") ),
        wxSizerFlags().Align(wxALIGN_CENTRE_HORIZONTAL).Border(wxALL & ~wxBOTTOM, 5));

    // 2) top: create wxTextCtrl with minimum size (100x60)
    topsizer->Add(
        new wxTextCtrl( p, wxID_ANY, wxT("My text (wxEXPAND)."), wxDefaultPosition, wxSize(100,60), wxTE_MULTILINE),
        wxSizerFlags(1).Expand().Border(wxALL, 5));

    // 2.5) Gratuitous test of wxStaticBoxSizers
    wxBoxSizer *statsizer = new wxStaticBoxSizer(
        new wxStaticBox(p, wxID_ANY, wxT("A wxStaticBoxSizer")), wxVERTICAL );
    statsizer->Add(
        new wxStaticText(p, wxID_ANY, wxT("And some TEXT inside it")),
        wxSizerFlags().Border(wxALL, 30));
    topsizer->Add(
        statsizer,
        wxSizerFlags(1).Expand().Border(wxALL, 10));

    // 2.7) And a test of wxGridSizer
    wxGridSizer *gridsizer = new wxGridSizer(2, 5, 5);
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("Grid sizer demo")),
                wxSizerFlags(1).Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Another label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("More text")),
                wxSizerFlags(1).Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Final label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("And yet more text")),
                wxSizerFlags().Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    topsizer->Add(
        gridsizer,
        wxSizerFlags().Proportion(1).Expand().Border(wxALL, 10));


#if wxUSE_STATLINE
    // 3) middle: create wxStaticLine with minimum size (3x3)
    topsizer->Add(
        new wxStaticLine( p, wxID_ANY, wxDefaultPosition, wxSize(3,3), wxHORIZONTAL),
        wxSizerFlags().Expand());
#endif // wxUSE_STATLINE


    // 4) bottom: create two centred wxButtons
    wxBoxSizer *button_box = new wxBoxSizer( wxHORIZONTAL );
    button_box->Add(
        new wxButton( p, wxID_ANY, wxT("Two buttons in a box") ),
        wxSizerFlags().Border(wxALL, 7));
    button_box->Add(
        new wxButton( p, wxID_ANY, wxT("(wxCENTER)") ),
        wxSizerFlags().Border(wxALL, 7));

    topsizer->Add(button_box, wxSizerFlags().Center());

    p->SetSizer( topsizer );

    // don't allow frame to get smaller than what the sizers tell it and also set
    // the initial size as calculated by the sizers
    topsizer->SetSizeHints( this );
}

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
    Close(true);
}

void MyFrame::TestProportions(wxCommandEvent& WXUNUSED(event))
{
    (new MyProportionsFrame(this))->Show();
}

void MyFrame::TestFlexSizers(wxCommandEvent& WXUNUSED(event) )
{
    (new MyFlexSizerFrame(this))->Show();
}

void MyFrame::TestNotebookSizers(wxCommandEvent& WXUNUSED(event) )
{
    MySizerDialog dialog( this, wxT("Notebook Sizer Test Dialog") );

    dialog.ShowModal();
}

void MyFrame::TestSetMinimal(wxCommandEvent& WXUNUSED(event) )
{
    (new MySimpleSizerFrame(this))->Show();
}

void MyFrame::TestNested(wxCommandEvent& WXUNUSED(event) )
{
    (new MyNestedSizerFrame(this))->Show();
}

void MyFrame::TestWrap(wxCommandEvent& WXUNUSED(event) )
{
    (new MyWrapSizerFrame(this))->Show();
}


void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event) )
{
    (void)wxMessageBox(wxT("wxWidgets GUI library layout demo\n"),
            wxT("About Layout Demo"), wxOK|wxICON_INFORMATION);
}

void MyFrame::TestGridBagSizer(wxCommandEvent& WXUNUSED(event) )
{
    (new MyGridBagSizerFrame(this))->Show();
}

// ----------------------------------------------------------------------------
// MyProportionsFrame
// ----------------------------------------------------------------------------

MyProportionsFrame::MyProportionsFrame(wxFrame *parent)
                  : wxFrame(parent, wxID_ANY, wxT("Box Sizer Proportions Demo"))
{
    size_t n;

    // create the controls
    wxPanel *panel = new wxPanel(this, wxID_ANY);
    for ( n = 0; n < WXSIZEOF(m_spins); n++ )
    {
        m_spins[n] = new wxSpinCtrl(panel);
        m_spins[n]->SetValue(n);
    }

    // lay them out
    m_sizer = new wxStaticBoxSizer(wxHORIZONTAL, panel,
                wxT("Try changing elements proportions and resizing the window"));
    for ( n = 0; n < WXSIZEOF(m_spins); n++ )
        m_sizer->Add(m_spins[n], wxSizerFlags().Border());

    // put everything together
    panel->SetSizer(m_sizer);
    wxSizer *sizerTop = new wxBoxSizer(wxVERTICAL);
    sizerTop->Add(panel, wxSizerFlags(1).Expand().Border());
    UpdateProportions();
    SetSizerAndFit(sizerTop);

    // and connect the events
    Connect(wxEVT_TEXT,
                wxCommandEventHandler(MyProportionsFrame::OnProportionUpdated));
    Connect(wxEVT_SPINCTRL,
            wxSpinEventHandler(MyProportionsFrame::OnProportionChanged));
}

void MyProportionsFrame::UpdateProportions()
{
    for ( size_t n = 0; n < WXSIZEOF(m_spins); n++ )
    {
        m_sizer->GetItem(n)->SetProportion(m_spins[n]->GetValue());
    }

    m_sizer->Layout();
}

void MyProportionsFrame::OnProportionUpdated(wxCommandEvent& WXUNUSED(event))
{
    UpdateProportions();
}

void MyProportionsFrame::OnProportionChanged(wxSpinEvent& WXUNUSED(event))
{
    UpdateProportions();
}

// ----------------------------------------------------------------------------
//  MyFlexSizerFrame
// ----------------------------------------------------------------------------

void MyFlexSizerFrame::InitFlexSizer(wxFlexGridSizer *sizer, wxWindow* parent)
{
    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 3; j++ )
        {
            wxWindow * const cell = new wxGenericStaticText
                                        (
                                            parent,
                                            wxID_ANY,
                                            wxString::Format("(%d, %d)",
                                                             i + 1, j + 1)
                                        );
            if ( (i + j) % 2 )
                cell->SetBackgroundColour( *wxLIGHT_GREY );
            sizer->Add(cell, 0, wxEXPAND | wxALIGN_CENTER_VERTICAL | wxALL, 3);
        }
    }
}

MyFlexSizerFrame::MyFlexSizerFrame(wxFrame* parent)
            : wxFrame(parent, wxID_ANY, "Flex Sizer Test Frame")
{
    wxFlexGridSizer *sizerFlex;
    wxPanel* p = new wxPanel(this, wxID_ANY);

    // consttuct the first column
    wxSizer *sizerCol1 = new wxBoxSizer(wxVERTICAL);
    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Ungrowable:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle column:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle row:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableRow(1);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("All growable columns:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(0, 1);
    sizerFlex->AddGrowableCol(1, 2);
    sizerFlex->AddGrowableCol(2, 3);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    // the second one
    wxSizer *sizerCol2 = new wxBoxSizer(wxVERTICAL);
    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle row and column:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with horz flex direction")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with grow mode == \"none\"")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerFlex->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_NONE);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with grow mode == \"all\"")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerFlex->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_ALL);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    // add both columns to grid sizer
    wxGridSizer *sizerTop = new wxGridSizer(2, 0, 20);
    sizerTop->Add(sizerCol1, 1, wxEXPAND);
    sizerTop->Add(sizerCol2, 1, wxEXPAND);

    p->SetSizer(sizerTop);
    sizerTop->SetSizeHints(this);
}

// ----------------------------------------------------------------------------
// MySizerDialog
// ----------------------------------------------------------------------------

MySizerDialog::MySizerDialog(wxWindow *parent, const wxString &title)
             : wxDialog(parent, wxID_ANY, wxString(title))
{
    // Begin with first hierarchy: a notebook at the top and
    // and OK button at the bottom.

    wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

    wxNotebook *notebook = new wxNotebook( this, wxID_ANY );
    topsizer->Add( notebook, 1, wxGROW );

    wxButton *button = new wxButton( this, wxID_OK, wxT("OK") );
    topsizer->Add( button, 0, wxALIGN_RIGHT | wxALL, 10 );

    // First page: one big text ctrl
    wxTextCtrl *multi = new wxTextCtrl( notebook, wxID_ANY, wxT("TextCtrl."), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE );
    notebook->AddPage( multi, wxT("Page One") );

    // Second page: a text ctrl and a button
    wxPanel *panel = new wxPanel( notebook, wxID_ANY );
    notebook->AddPage( panel, wxT("Page Two") );

    wxSizer *panelsizer = new wxBoxSizer( wxVERTICAL );

    wxTextCtrl *text = new wxTextCtrl( panel, wxID_ANY, wxT("TextLine 1."), wxDefaultPosition, wxSize(250,wxDefaultCoord) );
    panelsizer->Add( text, 0, wxGROW|wxALL, 30 );
    text = new wxTextCtrl( panel, wxID_ANY, wxT("TextLine 2."), wxDefaultPosition, wxSize(250,wxDefaultCoord) );
    panelsizer->Add( text, 0, wxGROW|wxALL, 30 );
    wxButton *button2 = new wxButton( panel, wxID_ANY, wxT("Hallo") );
    panelsizer->Add( button2, 0, wxALIGN_RIGHT | wxLEFT|wxRIGHT|wxBOTTOM, 30 );

    panel->SetSizer( panelsizer );

    // Tell dialog to use sizer
    SetSizerAndFit( topsizer );
}

// ----------------------------------------------------------------------------
// MyGridBagSizerFrame
// ----------------------------------------------------------------------------

// some simple macros to help make the sample code below more clear
#define TEXTCTRL(text)   new wxTextCtrl(p, wxID_ANY, wxT(text))
#define MLTEXTCTRL(text) new wxTextCtrl(p, wxID_ANY, wxT(text), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE)
#define POS(r, c)        wxGBPosition(r,c)
#define SPAN(r, c)       wxGBSpan(r,c)

const wxChar gbsDescription[] =wxT("\
The wxGridBagSizer is similar to the wxFlexGridSizer except the items are explicitly positioned\n\
in a virtual cell of the layout grid, and column or row spanning is allowed.  For example, this\n\
static text is positioned at (0,0) and it spans 7 columns.");


// Some IDs
enum {
    GBS_HIDE_BTN = 1212,
    GBS_SHOW_BTN,
    GBS_MOVE_BTN1,
    GBS_MOVE_BTN2,

    GBS_MAX
};


BEGIN_EVENT_TABLE(MyGridBagSizerFrame, wxFrame)
    EVT_BUTTON( GBS_HIDE_BTN,  MyGridBagSizerFrame::OnHideBtn)
    EVT_BUTTON( GBS_SHOW_BTN,  MyGridBagSizerFrame::OnShowBtn)
    EVT_BUTTON( GBS_MOVE_BTN1, MyGridBagSizerFrame::OnMoveBtn)
    EVT_BUTTON( GBS_MOVE_BTN2, MyGridBagSizerFrame::OnMoveBtn)
END_EVENT_TABLE()


MyGridBagSizerFrame::MyGridBagSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "wxGridBagSizer Test Frame")
{
    wxPanel* p = new wxPanel(this, wxID_ANY);
    m_panel = p;
    m_gbs = new wxGridBagSizer();


    m_gbs->Add( new wxStaticText(p, wxID_ANY, gbsDescription),
                POS(0,0), SPAN(1, 7),
                wxALIGN_CENTER | wxALL, 5);

    m_gbs->Add( TEXTCTRL("pos(1,0)"),   POS(1,0) );
    m_gbs->Add( TEXTCTRL("pos(1,1)"),   POS(1,1) );
    m_gbs->Add( TEXTCTRL("pos(2,0)"),   POS(2,0) );
    m_gbs->Add( TEXTCTRL("pos(2,1)"),   POS(2,1) );
    m_gbs->Add( MLTEXTCTRL("pos(3,2), span(1,2)\nthis row and col are growable"),
              POS(3,2), SPAN(1,2), wxEXPAND );
    m_gbs->Add( MLTEXTCTRL("pos(4,3)\nspan(3,1)"),
              POS(4,3), SPAN(3,1), wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(5,4)"),   POS(5,4), wxDefaultSpan, wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(6,5)"),   POS(6,5), wxDefaultSpan, wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(7,6)"),   POS(7,6) );

    //m_gbs->Add( TEXTCTRL("bad position"), POS(4,3) );  // Test for assert
    //m_gbs->Add( TEXTCTRL("bad position"), POS(5,3) );  // Test for assert


    m_moveBtn1 = new wxButton(p, GBS_MOVE_BTN1, wxT("Move this to (3,6)"));
    m_moveBtn2 = new wxButton(p, GBS_MOVE_BTN2, wxT("Move this to (3,6)"));
    m_gbs->Add( m_moveBtn1, POS(10,2) );
    m_gbs->Add( m_moveBtn2, POS(10,3) );

    m_hideBtn = new wxButton(p, GBS_HIDE_BTN, wxT("Hide this item -->"));
    m_gbs->Add(m_hideBtn, POS(12, 3));

    m_hideTxt = new wxTextCtrl(p, wxID_ANY, wxT("pos(12,4), size(150, wxDefaultCoord)"),
                                wxDefaultPosition, wxSize(150,wxDefaultCoord));
    m_gbs->Add( m_hideTxt, POS(12,4) );

    m_showBtn = new wxButton(p, GBS_SHOW_BTN, wxT("<-- Show it again"));
    m_gbs->Add(m_showBtn, POS(12, 5));
    m_showBtn->Disable();

    m_gbs->Add(10,10, POS(14,0));

    m_gbs->AddGrowableRow(3);
    m_gbs->AddGrowableCol(2);

    p->SetSizerAndFit(m_gbs);
    SetClientSize(p->GetSize());
}


void MyGridBagSizerFrame::OnHideBtn(wxCommandEvent&)
{
    m_gbs->Hide(m_hideTxt);
    m_hideBtn->Disable();
    m_showBtn->Enable();
    m_gbs->Layout();
}

void MyGridBagSizerFrame::OnShowBtn(wxCommandEvent&)
{
    m_gbs->Show(m_hideTxt);
    m_hideBtn->Enable();
    m_showBtn->Disable();
    m_gbs->Layout();
}


void MyGridBagSizerFrame::OnMoveBtn(wxCommandEvent& event)
{
    wxButton* btn = (wxButton*)event.GetEventObject();
    wxGBPosition curPos = m_gbs->GetItemPosition(btn);

    // if it's already at the "other" spot then move it back
    if (curPos == wxGBPosition(3,6))
    {
        m_gbs->SetItemPosition(btn, m_lastPos);
        btn->SetLabel(wxT("Move this to (3,6)"));
    }
    else
    {
        if ( m_gbs->CheckForIntersection(wxGBPosition(3,6), wxGBSpan(1,1)) )
            wxMessageBox(
wxT("wxGridBagSizer will not allow items to be in the same cell as\n\
another item, so this operation will fail.  You will also get an assert\n\
when compiled in debug mode."), wxT("Warning"), wxOK | wxICON_INFORMATION);

        if ( m_gbs->SetItemPosition(btn, wxGBPosition(3,6)) )
        {
            m_lastPos = curPos;
            btn->SetLabel(wxT("Move it back"));
        }
    }
    m_gbs->Layout();
}

// ----------------------------------------------------------------------------
// MySimpleSizerFrame
// ----------------------------------------------------------------------------

// Some IDs
enum {
    ID_SET_SMALL = 1300,
    ID_SET_BIG
};

BEGIN_EVENT_TABLE(MySimpleSizerFrame, wxFrame)
    EVT_MENU( ID_SET_SMALL, MySimpleSizerFrame::OnSetSmallSize)
    EVT_MENU( ID_SET_BIG, MySimpleSizerFrame::OnSetBigSize)
END_EVENT_TABLE()

MySimpleSizerFrame::MySimpleSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Simple Sizer Test Frame")
{
    wxMenu *menu = new wxMenu;

    menu->Append(ID_SET_SMALL, wxT("Make text control small\tF4"));
    menu->Append(ID_SET_BIG, wxT("Make text control big\tF5"));

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, wxT("&File"));

    SetMenuBar( menu_bar );

    wxBoxSizer *main_sizer = new wxBoxSizer( wxHORIZONTAL );

    m_target = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 80, wxDefaultCoord ) );
    main_sizer->Add( m_target, 1, wxALL, 5 );

    main_sizer->Add( new wxStaticText( this, wxID_ANY, wxT("Set alternating sizes using F4 and F5") ), 0, wxALL, 5 );

    SetSizer( main_sizer);

    Layout();
    GetSizer()->Fit( this );
}

void MySimpleSizerFrame::OnSetSmallSize( wxCommandEvent& WXUNUSED(event))
{
    GetSizer()->SetItemMinSize( m_target, 40, -1 );
    Layout();
    GetSizer()->Fit( this );
}

void MySimpleSizerFrame::OnSetBigSize( wxCommandEvent& WXUNUSED(event))
{
    GetSizer()->SetItemMinSize( m_target, 140, -1 );
    Layout();
    GetSizer()->Fit( this );
}


// ----------------------------------------------------------------------------
// MyNestedSizerFrame
// ----------------------------------------------------------------------------


MyNestedSizerFrame::MyNestedSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Nested Sizer Test Frame")
{
    wxMenu *menu = new wxMenu;

    menu->Append(wxID_ABOUT, wxT("Do nothing"));

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, wxT("&File"));

    SetMenuBar( menu_bar );

    wxBoxSizer *main_sizer = new wxBoxSizer( wxVERTICAL );

    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );

    wxPanel *panel = new wxPanel( this, -1, wxDefaultPosition, wxDefaultSize,
                                  wxTAB_TRAVERSAL | wxSUNKEN_BORDER );
    main_sizer->Add( panel, 0, wxALIGN_CENTER );
    wxBoxSizer *panel_sizer = new wxBoxSizer( wxVERTICAL );
    panel->SetSizer( panel_sizer );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );

    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );

    m_target = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 80, wxDefaultCoord ) );
    main_sizer->Add( m_target, 1, wxALL|wxGROW, 5 );

    SetSizerAndFit( main_sizer);
}


// ----------------------------------------------------------------------------
// MyWrapSizerFrame
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(MyWrapSizerFrame, wxFrame)
    EVT_MENU(wxID_ADD, MyWrapSizerFrame::OnAddCheckbox)
    EVT_MENU(wxID_REMOVE, MyWrapSizerFrame::OnRemoveCheckbox)
END_EVENT_TABLE()

MyWrapSizerFrame::MyWrapSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Wrap Sizer Test Frame",
              wxDefaultPosition, wxSize(200,-1))
{
    wxMenu *menu = new wxMenu;

    menu->Append(wxID_ADD, "&Add a checkbox\tCtrl-+");
    menu->Append(wxID_REMOVE, "&Remove a checkbox\tCtrl--");

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, "&Wrap sizer");

    SetMenuBar( menu_bar );

    wxBoxSizer *root = new wxBoxSizer( wxVERTICAL );

    wxStaticBoxSizer *topSizer = new wxStaticBoxSizer( wxVERTICAL, this, "Wrapping check-boxes" );
    m_checkboxParent = topSizer->GetStaticBox();
    m_wrapSizer = new wxWrapSizer(wxHORIZONTAL);

    // A number of checkboxes inside a wrap sizer
    for( int i = 0; i < 6; i++ )
        DoAddCheckbox();

    topSizer->Add( m_wrapSizer, wxSizerFlags(1).Expand());
    root->Add( topSizer, wxSizerFlags().Expand().Border());

    // A shaped item inside a box sizer
    wxSizer *bottomSizer = new wxStaticBoxSizer( wxVERTICAL, this, "With wxSHAPED item" );
    wxSizer *horzBoxSizer = new wxBoxSizer(wxHORIZONTAL);
    bottomSizer->Add( horzBoxSizer, 100, wxEXPAND );
    horzBoxSizer->Add( new wxListBox(this,wxID_ANY,wxPoint(0,0),wxSize(70,70)), 0, wxEXPAND|wxSHAPED );
    horzBoxSizer->Add( 10,10 );
    horzBoxSizer->Add( new wxCheckBox(this,wxID_ANY,"A much longer option..."), 100, 0, 5 );

    root->Add( bottomSizer, 1, wxEXPAND | wxALL, 5 );

    // Set sizer for window
    SetSizerAndFit( root );
}

void MyWrapSizerFrame::DoAddCheckbox()
{
    m_wrapSizer->Add(new wxCheckBox
                         (
                            m_checkboxParent,
                            wxID_ANY,
                            wxString::Format
                            (
                                "Option %d",
                                (int)m_wrapSizer->GetItemCount() + 1
                            )
                         ),
                     wxSizerFlags().Centre().Border());
}

void MyWrapSizerFrame::OnAddCheckbox(wxCommandEvent& WXUNUSED(event))
{
    DoAddCheckbox();
    Layout();
}

void MyWrapSizerFrame::OnRemoveCheckbox(wxCommandEvent& WXUNUSED(event))
{
    if ( m_wrapSizer->IsEmpty() )
    {
        wxLogStatus(this, "No more checkboxes to remove.");
        return;
    }

    delete m_wrapSizer->GetItem(m_wrapSizer->GetItemCount() - 1)->GetWindow();
    Layout();
}


IMPLEMENT_APP(MyApp)
     // NO wx_gui_main(argc,argv);
#else

#if CMD_LINE 
int main(int argc,char *argv[])
{
    command_line_main(argc,argv);
    return 0;
}
#endif

#if WEB_SERVER
int main(int argc,char *argv[])
{
    web_server_main(argc,argv);
    return 0;
}
#endif

/*
=HYPERLINK("https://cgwb-test.nci.nih.gov:8443/cgi-bin/alview?position="&D1&":"&E1-100&"-"&E1+100&"&iw=1000&ih=400&file="&A1,"n")
*/

#endif

#if WX_GUI
// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/sizer.h"
#include "wx/gbsizer.h"
#include "wx/statline.h"
#include "wx/notebook.h"
#include "wx/spinctrl.h"
#include "wx/wrapsizer.h"
#include "wx/generic/stattextg.h"

#include "layout.h"

#ifndef wxHAS_IMAGES_IN_RESOURCES
    //rpf #include "../sample.xpm"
    #include "sample.xpm"
#endif


// ----------------------------------------------------------------------------
// MyApp
// ----------------------------------------------------------------------------

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
  if ( !wxApp::OnInit() )
      return false;

  // Create the MAIN frame window
  MyFrame *frame = new MyFrame;

  frame->Show(true);

  return true;
}

// ----------------------------------------------------------------------------
// MyFrame
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
  EVT_MENU(LAYOUT_ABOUT, MyFrame::OnAbout)
  EVT_MENU(LAYOUT_QUIT, MyFrame::OnQuit)

  EVT_MENU(LAYOUT_TEST_PROPORTIONS, MyFrame::TestProportions)
  EVT_MENU(LAYOUT_TEST_SIZER, MyFrame::TestFlexSizers)
  EVT_MENU(LAYOUT_TEST_NB_SIZER, MyFrame::TestNotebookSizers)
  EVT_MENU(LAYOUT_TEST_GB_SIZER, MyFrame::TestGridBagSizer)
  EVT_MENU(LAYOUT_TEST_SET_MINIMAL, MyFrame::TestSetMinimal)
  EVT_MENU(LAYOUT_TEST_NESTED, MyFrame::TestNested)
  EVT_MENU(LAYOUT_TEST_WRAP, MyFrame::TestWrap)
END_EVENT_TABLE()

// Define my frame constructor
MyFrame::MyFrame()
       : wxFrame(NULL, wxID_ANY, wxT("wxWidgets Layout Demo"))
{
    SetIcon(wxICON(sample));

    // Make a menubar
    wxMenu *file_menu = new wxMenu;

    file_menu->Append(LAYOUT_TEST_PROPORTIONS, wxT("&Proportions demo...\tF1"));
    file_menu->Append(LAYOUT_TEST_SIZER, wxT("Test wx&FlexSizer...\tF2"));
    file_menu->Append(LAYOUT_TEST_NB_SIZER, wxT("Test &notebook sizers...\tF3"));
    file_menu->Append(LAYOUT_TEST_GB_SIZER, wxT("Test &gridbag sizer...\tF4"));
    file_menu->Append(LAYOUT_TEST_SET_MINIMAL, wxT("Test Set&ItemMinSize...\tF5"));
    file_menu->Append(LAYOUT_TEST_NESTED, wxT("Test nested sizer in a wxPanel...\tF6"));
    file_menu->Append(LAYOUT_TEST_WRAP, wxT("Test wrap sizers...\tF7"));

    file_menu->AppendSeparator();
    file_menu->Append(LAYOUT_QUIT, wxT("E&xit"), wxT("Quit program"));

    wxMenu *help_menu = new wxMenu;
    help_menu->Append(LAYOUT_ABOUT, wxT("&About"), wxT("About layout demo..."));

    wxMenuBar *menu_bar = new wxMenuBar;

    menu_bar->Append(file_menu, wxT("&File"));
    menu_bar->Append(help_menu, wxT("&Help"));

    // Associate the menu bar with the frame
    SetMenuBar(menu_bar);

#if wxUSE_STATUSBAR
    CreateStatusBar(2);
    SetStatusText(wxT("wxWidgets layout demo"));
#endif // wxUSE_STATUSBAR

    wxPanel* p = new wxPanel(this, wxID_ANY);

    // we want to get a dialog that is stretchable because it
    // has a text ctrl in the middle. at the bottom, we have
    // two buttons which.

    wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

    // 1) top: create wxStaticText with minimum size equal to its default size
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_RIGHT).") ),
        wxSizerFlags().Align(wxALIGN_RIGHT).Border(wxALL & ~wxBOTTOM, 5));
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_LEFT).") ),
        wxSizerFlags().Align(wxALIGN_LEFT).Border(wxALL & ~wxBOTTOM, 5));
    topsizer->Add(
        new wxStaticText( p, wxID_ANY, wxT("An explanation (wxALIGN_CENTRE_HORIZONTAL).") ),
        wxSizerFlags().Align(wxALIGN_CENTRE_HORIZONTAL).Border(wxALL & ~wxBOTTOM, 5));

    // 2) top: create wxTextCtrl with minimum size (100x60)
    topsizer->Add(
        new wxTextCtrl( p, wxID_ANY, wxT("My text (wxEXPAND)."), wxDefaultPosition, wxSize(100,60), wxTE_MULTILINE),
        wxSizerFlags(1).Expand().Border(wxALL, 5));

    // 2.5) Gratuitous test of wxStaticBoxSizers
    wxBoxSizer *statsizer = new wxStaticBoxSizer(
        new wxStaticBox(p, wxID_ANY, wxT("A wxStaticBoxSizer")), wxVERTICAL );
    statsizer->Add(
        new wxStaticText(p, wxID_ANY, wxT("And some TEXT inside it")),
        wxSizerFlags().Border(wxALL, 30));
    topsizer->Add(
        statsizer,
        wxSizerFlags(1).Expand().Border(wxALL, 10));

    // 2.7) And a test of wxGridSizer
    wxGridSizer *gridsizer = new wxGridSizer(2, 5, 5);
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("Grid sizer demo")),
                wxSizerFlags(1).Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Another label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("More text")),
                wxSizerFlags(1).Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxStaticText(p, wxID_ANY, wxT("Final label")),
                wxSizerFlags().Align(wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL));
    gridsizer->Add(new wxTextCtrl(p, wxID_ANY, wxT("And yet more text")),
                wxSizerFlags().Align(wxGROW | wxALIGN_CENTER_VERTICAL));
    topsizer->Add(
        gridsizer,
        wxSizerFlags().Proportion(1).Expand().Border(wxALL, 10));


#if wxUSE_STATLINE
    // 3) middle: create wxStaticLine with minimum size (3x3)
    topsizer->Add(
        new wxStaticLine( p, wxID_ANY, wxDefaultPosition, wxSize(3,3), wxHORIZONTAL),
        wxSizerFlags().Expand());
#endif // wxUSE_STATLINE


    // 4) bottom: create two centred wxButtons
    wxBoxSizer *button_box = new wxBoxSizer( wxHORIZONTAL );
    button_box->Add(
        new wxButton( p, wxID_ANY, wxT("Two buttons in a box") ),
        wxSizerFlags().Border(wxALL, 7));
    button_box->Add(
        new wxButton( p, wxID_ANY, wxT("(wxCENTER)") ),
        wxSizerFlags().Border(wxALL, 7));

    topsizer->Add(button_box, wxSizerFlags().Center());

    p->SetSizer( topsizer );

    // don't allow frame to get smaller than what the sizers tell it and also set
    // the initial size as calculated by the sizers
    topsizer->SetSizeHints( this );
}

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
    Close(true);
}

void MyFrame::TestProportions(wxCommandEvent& WXUNUSED(event))
{
    (new MyProportionsFrame(this))->Show();
}

void MyFrame::TestFlexSizers(wxCommandEvent& WXUNUSED(event) )
{
    (new MyFlexSizerFrame(this))->Show();
}

void MyFrame::TestNotebookSizers(wxCommandEvent& WXUNUSED(event) )
{
    MySizerDialog dialog( this, wxT("Notebook Sizer Test Dialog") );

    dialog.ShowModal();
}

void MyFrame::TestSetMinimal(wxCommandEvent& WXUNUSED(event) )
{
    (new MySimpleSizerFrame(this))->Show();
}

void MyFrame::TestNested(wxCommandEvent& WXUNUSED(event) )
{
    (new MyNestedSizerFrame(this))->Show();
}

void MyFrame::TestWrap(wxCommandEvent& WXUNUSED(event) )
{
    (new MyWrapSizerFrame(this))->Show();
}


void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event) )
{
    (void)wxMessageBox(wxT("wxWidgets GUI library layout demo\n"),
            wxT("About Layout Demo"), wxOK|wxICON_INFORMATION);
}

void MyFrame::TestGridBagSizer(wxCommandEvent& WXUNUSED(event) )
{
    (new MyGridBagSizerFrame(this))->Show();
}

// ----------------------------------------------------------------------------
// MyProportionsFrame
// ----------------------------------------------------------------------------

MyProportionsFrame::MyProportionsFrame(wxFrame *parent)
                  : wxFrame(parent, wxID_ANY, wxT("Box Sizer Proportions Demo"))
{
    size_t n;

    // create the controls
    wxPanel *panel = new wxPanel(this, wxID_ANY);
    for ( n = 0; n < WXSIZEOF(m_spins); n++ )
    {
        m_spins[n] = new wxSpinCtrl(panel);
        m_spins[n]->SetValue(n);
    }

    // lay them out
    m_sizer = new wxStaticBoxSizer(wxHORIZONTAL, panel,
                wxT("Try changing elements proportions and resizing the window"));
    for ( n = 0; n < WXSIZEOF(m_spins); n++ )
        m_sizer->Add(m_spins[n], wxSizerFlags().Border());

    // put everything together
    panel->SetSizer(m_sizer);
    wxSizer *sizerTop = new wxBoxSizer(wxVERTICAL);
    sizerTop->Add(panel, wxSizerFlags(1).Expand().Border());
    UpdateProportions();
    SetSizerAndFit(sizerTop);

    // and connect the events
    Connect(wxEVT_TEXT,
                wxCommandEventHandler(MyProportionsFrame::OnProportionUpdated));
    Connect(wxEVT_SPINCTRL,
            wxSpinEventHandler(MyProportionsFrame::OnProportionChanged));
}

void MyProportionsFrame::UpdateProportions()
{
    for ( size_t n = 0; n < WXSIZEOF(m_spins); n++ )
    {
        m_sizer->GetItem(n)->SetProportion(m_spins[n]->GetValue());
    }

    m_sizer->Layout();
}

void MyProportionsFrame::OnProportionUpdated(wxCommandEvent& WXUNUSED(event))
{
    UpdateProportions();
}

void MyProportionsFrame::OnProportionChanged(wxSpinEvent& WXUNUSED(event))
{
    UpdateProportions();
}

// ----------------------------------------------------------------------------
//  MyFlexSizerFrame
// ----------------------------------------------------------------------------

void MyFlexSizerFrame::InitFlexSizer(wxFlexGridSizer *sizer, wxWindow* parent)
{
    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 3; j++ )
        {
            wxWindow * const cell = new wxGenericStaticText
                                        (
                                            parent,
                                            wxID_ANY,
                                            wxString::Format("(%d, %d)",
                                                             i + 1, j + 1)
                                        );
            if ( (i + j) % 2 )
                cell->SetBackgroundColour( *wxLIGHT_GREY );
            sizer->Add(cell, 0, wxEXPAND | wxALIGN_CENTER_VERTICAL | wxALL, 3);
        }
    }
}

MyFlexSizerFrame::MyFlexSizerFrame(wxFrame* parent)
            : wxFrame(parent, wxID_ANY, "Flex Sizer Test Frame")
{
    wxFlexGridSizer *sizerFlex;
    wxPanel* p = new wxPanel(this, wxID_ANY);

    // consttuct the first column
    wxSizer *sizerCol1 = new wxBoxSizer(wxVERTICAL);
    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Ungrowable:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle column:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle row:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableRow(1);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol1->Add(new wxStaticText(p, wxID_ANY, wxT("All growable columns:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(0, 1);
    sizerFlex->AddGrowableCol(1, 2);
    sizerFlex->AddGrowableCol(2, 3);
    sizerCol1->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    // the second one
    wxSizer *sizerCol2 = new wxBoxSizer(wxVERTICAL);
    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Growable middle row and column:")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with horz flex direction")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with grow mode == \"none\"")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerFlex->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_NONE);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    sizerCol2->Add(new wxStaticText(p, wxID_ANY, wxT("Same with grow mode == \"all\"")), 0, wxCENTER | wxTOP, 20);
    sizerFlex = new wxFlexGridSizer(3, 3, wxSize(5, 5));
    InitFlexSizer(sizerFlex, p);
    sizerFlex->AddGrowableCol(1);
    sizerFlex->AddGrowableRow(1);
    sizerFlex->SetFlexibleDirection(wxHORIZONTAL);
    sizerFlex->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_ALL);
    sizerCol2->Add(sizerFlex, 1, wxALL | wxEXPAND, 10);

    // add both columns to grid sizer
    wxGridSizer *sizerTop = new wxGridSizer(2, 0, 20);
    sizerTop->Add(sizerCol1, 1, wxEXPAND);
    sizerTop->Add(sizerCol2, 1, wxEXPAND);

    p->SetSizer(sizerTop);
    sizerTop->SetSizeHints(this);
}

// ----------------------------------------------------------------------------
// MySizerDialog
// ----------------------------------------------------------------------------

MySizerDialog::MySizerDialog(wxWindow *parent, const wxString &title)
             : wxDialog(parent, wxID_ANY, wxString(title))
{
    // Begin with first hierarchy: a notebook at the top and
    // and OK button at the bottom.

    wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

    wxNotebook *notebook = new wxNotebook( this, wxID_ANY );
    topsizer->Add( notebook, 1, wxGROW );

    wxButton *button = new wxButton( this, wxID_OK, wxT("OK") );
    topsizer->Add( button, 0, wxALIGN_RIGHT | wxALL, 10 );

    // First page: one big text ctrl
    wxTextCtrl *multi = new wxTextCtrl( notebook, wxID_ANY, wxT("TextCtrl."), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE );
    notebook->AddPage( multi, wxT("Page One") );

    // Second page: a text ctrl and a button
    wxPanel *panel = new wxPanel( notebook, wxID_ANY );
    notebook->AddPage( panel, wxT("Page Two") );

    wxSizer *panelsizer = new wxBoxSizer( wxVERTICAL );

    wxTextCtrl *text = new wxTextCtrl( panel, wxID_ANY, wxT("TextLine 1."), wxDefaultPosition, wxSize(250,wxDefaultCoord) );
    panelsizer->Add( text, 0, wxGROW|wxALL, 30 );
    text = new wxTextCtrl( panel, wxID_ANY, wxT("TextLine 2."), wxDefaultPosition, wxSize(250,wxDefaultCoord) );
    panelsizer->Add( text, 0, wxGROW|wxALL, 30 );
    wxButton *button2 = new wxButton( panel, wxID_ANY, wxT("Hallo") );
    panelsizer->Add( button2, 0, wxALIGN_RIGHT | wxLEFT|wxRIGHT|wxBOTTOM, 30 );

    panel->SetSizer( panelsizer );

    // Tell dialog to use sizer
    SetSizerAndFit( topsizer );
}

// ----------------------------------------------------------------------------
// MyGridBagSizerFrame
// ----------------------------------------------------------------------------

// some simple macros to help make the sample code below more clear
#define TEXTCTRL(text)   new wxTextCtrl(p, wxID_ANY, wxT(text))
#define MLTEXTCTRL(text) new wxTextCtrl(p, wxID_ANY, wxT(text), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE)
#define POS(r, c)        wxGBPosition(r,c)
#define SPAN(r, c)       wxGBSpan(r,c)

const wxChar gbsDescription[] =wxT("\
The wxGridBagSizer is similar to the wxFlexGridSizer except the items are explicitly positioned\n\
in a virtual cell of the layout grid, and column or row spanning is allowed.  For example, this\n\
static text is positioned at (0,0) and it spans 7 columns.");


// Some IDs
enum {
    GBS_HIDE_BTN = 1212,
    GBS_SHOW_BTN,
    GBS_MOVE_BTN1,
    GBS_MOVE_BTN2,

    GBS_MAX
};


BEGIN_EVENT_TABLE(MyGridBagSizerFrame, wxFrame)
    EVT_BUTTON( GBS_HIDE_BTN,  MyGridBagSizerFrame::OnHideBtn)
    EVT_BUTTON( GBS_SHOW_BTN,  MyGridBagSizerFrame::OnShowBtn)
    EVT_BUTTON( GBS_MOVE_BTN1, MyGridBagSizerFrame::OnMoveBtn)
    EVT_BUTTON( GBS_MOVE_BTN2, MyGridBagSizerFrame::OnMoveBtn)
END_EVENT_TABLE()


MyGridBagSizerFrame::MyGridBagSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "wxGridBagSizer Test Frame")
{
    wxPanel* p = new wxPanel(this, wxID_ANY);
    m_panel = p;
    m_gbs = new wxGridBagSizer();


    m_gbs->Add( new wxStaticText(p, wxID_ANY, gbsDescription),
                POS(0,0), SPAN(1, 7),
                wxALIGN_CENTER | wxALL, 5);

    m_gbs->Add( TEXTCTRL("pos(1,0)"),   POS(1,0) );
    m_gbs->Add( TEXTCTRL("pos(1,1)"),   POS(1,1) );
    m_gbs->Add( TEXTCTRL("pos(2,0)"),   POS(2,0) );
    m_gbs->Add( TEXTCTRL("pos(2,1)"),   POS(2,1) );
    m_gbs->Add( MLTEXTCTRL("pos(3,2), span(1,2)\nthis row and col are growable"),
              POS(3,2), SPAN(1,2), wxEXPAND );
    m_gbs->Add( MLTEXTCTRL("pos(4,3)\nspan(3,1)"),
              POS(4,3), SPAN(3,1), wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(5,4)"),   POS(5,4), wxDefaultSpan, wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(6,5)"),   POS(6,5), wxDefaultSpan, wxEXPAND );
    m_gbs->Add( TEXTCTRL("pos(7,6)"),   POS(7,6) );

    //m_gbs->Add( TEXTCTRL("bad position"), POS(4,3) );  // Test for assert
    //m_gbs->Add( TEXTCTRL("bad position"), POS(5,3) );  // Test for assert


    m_moveBtn1 = new wxButton(p, GBS_MOVE_BTN1, wxT("Move this to (3,6)"));
    m_moveBtn2 = new wxButton(p, GBS_MOVE_BTN2, wxT("Move this to (3,6)"));
    m_gbs->Add( m_moveBtn1, POS(10,2) );
    m_gbs->Add( m_moveBtn2, POS(10,3) );

    m_hideBtn = new wxButton(p, GBS_HIDE_BTN, wxT("Hide this item -->"));
    m_gbs->Add(m_hideBtn, POS(12, 3));

    m_hideTxt = new wxTextCtrl(p, wxID_ANY, wxT("pos(12,4), size(150, wxDefaultCoord)"),
                                wxDefaultPosition, wxSize(150,wxDefaultCoord));
    m_gbs->Add( m_hideTxt, POS(12,4) );

    m_showBtn = new wxButton(p, GBS_SHOW_BTN, wxT("<-- Show it again"));
    m_gbs->Add(m_showBtn, POS(12, 5));
    m_showBtn->Disable();

    m_gbs->Add(10,10, POS(14,0));

    m_gbs->AddGrowableRow(3);
    m_gbs->AddGrowableCol(2);

    p->SetSizerAndFit(m_gbs);
    SetClientSize(p->GetSize());
}


void MyGridBagSizerFrame::OnHideBtn(wxCommandEvent&)
{
    m_gbs->Hide(m_hideTxt);
    m_hideBtn->Disable();
    m_showBtn->Enable();
    m_gbs->Layout();
}

void MyGridBagSizerFrame::OnShowBtn(wxCommandEvent&)
{
    m_gbs->Show(m_hideTxt);
    m_hideBtn->Enable();
    m_showBtn->Disable();
    m_gbs->Layout();
}


void MyGridBagSizerFrame::OnMoveBtn(wxCommandEvent& event)
{
    wxButton* btn = (wxButton*)event.GetEventObject();
    wxGBPosition curPos = m_gbs->GetItemPosition(btn);

    // if it's already at the "other" spot then move it back
    if (curPos == wxGBPosition(3,6))
    {
        m_gbs->SetItemPosition(btn, m_lastPos);
        btn->SetLabel(wxT("Move this to (3,6)"));
    }
    else
    {
        if ( m_gbs->CheckForIntersection(wxGBPosition(3,6), wxGBSpan(1,1)) )
            wxMessageBox(
wxT("wxGridBagSizer will not allow items to be in the same cell as\n\
another item, so this operation will fail.  You will also get an assert\n\
when compiled in debug mode."), wxT("Warning"), wxOK | wxICON_INFORMATION);

        if ( m_gbs->SetItemPosition(btn, wxGBPosition(3,6)) )
        {
            m_lastPos = curPos;
            btn->SetLabel(wxT("Move it back"));
        }
    }
    m_gbs->Layout();
}

// ----------------------------------------------------------------------------
// MySimpleSizerFrame
// ----------------------------------------------------------------------------

// Some IDs
enum {
    ID_SET_SMALL = 1300,
    ID_SET_BIG
};

BEGIN_EVENT_TABLE(MySimpleSizerFrame, wxFrame)
    EVT_MENU( ID_SET_SMALL, MySimpleSizerFrame::OnSetSmallSize)
    EVT_MENU( ID_SET_BIG, MySimpleSizerFrame::OnSetBigSize)
END_EVENT_TABLE()

MySimpleSizerFrame::MySimpleSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Simple Sizer Test Frame")
{
    wxMenu *menu = new wxMenu;

    menu->Append(ID_SET_SMALL, wxT("Make text control small\tF4"));
    menu->Append(ID_SET_BIG, wxT("Make text control big\tF5"));

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, wxT("&File"));

    SetMenuBar( menu_bar );

    wxBoxSizer *main_sizer = new wxBoxSizer( wxHORIZONTAL );

    m_target = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 80, wxDefaultCoord ) );
    main_sizer->Add( m_target, 1, wxALL, 5 );

    main_sizer->Add( new wxStaticText( this, wxID_ANY, wxT("Set alternating sizes using F4 and F5") ), 0, wxALL, 5 );

    SetSizer( main_sizer);

    Layout();
    GetSizer()->Fit( this );
}

void MySimpleSizerFrame::OnSetSmallSize( wxCommandEvent& WXUNUSED(event))
{
    GetSizer()->SetItemMinSize( m_target, 40, -1 );
    Layout();
    GetSizer()->Fit( this );
}

void MySimpleSizerFrame::OnSetBigSize( wxCommandEvent& WXUNUSED(event))
{
    GetSizer()->SetItemMinSize( m_target, 140, -1 );
    Layout();
    GetSizer()->Fit( this );
}


// ----------------------------------------------------------------------------
// MyNestedSizerFrame
// ----------------------------------------------------------------------------


MyNestedSizerFrame::MyNestedSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Nested Sizer Test Frame")
{
    wxMenu *menu = new wxMenu;

    menu->Append(wxID_ABOUT, wxT("Do nothing"));

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, wxT("&File"));

    SetMenuBar( menu_bar );

    wxBoxSizer *main_sizer = new wxBoxSizer( wxVERTICAL );

    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );
    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );

    wxPanel *panel = new wxPanel( this, -1, wxDefaultPosition, wxDefaultSize,
                                  wxTAB_TRAVERSAL | wxSUNKEN_BORDER );
    main_sizer->Add( panel, 0, wxALIGN_CENTER );
    wxBoxSizer *panel_sizer = new wxBoxSizer( wxVERTICAL );
    panel->SetSizer( panel_sizer );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );
    panel_sizer->Add( new wxStaticText( panel, -1, wxT("Hello inside") ) );

    main_sizer->Add( new wxStaticText( this, -1, wxT("Hello outside") ), 0, wxALIGN_CENTER );

    m_target = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 80, wxDefaultCoord ) );
    main_sizer->Add( m_target, 1, wxALL|wxGROW, 5 );

    SetSizerAndFit( main_sizer);
}


// ----------------------------------------------------------------------------
// MyWrapSizerFrame
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(MyWrapSizerFrame, wxFrame)
    EVT_MENU(wxID_ADD, MyWrapSizerFrame::OnAddCheckbox)
    EVT_MENU(wxID_REMOVE, MyWrapSizerFrame::OnRemoveCheckbox)
END_EVENT_TABLE()

MyWrapSizerFrame::MyWrapSizerFrame(wxFrame* parent)
    : wxFrame(parent, wxID_ANY, "Wrap Sizer Test Frame",
              wxDefaultPosition, wxSize(200,-1))
{
    wxMenu *menu = new wxMenu;

    menu->Append(wxID_ADD, "&Add a checkbox\tCtrl-+");
    menu->Append(wxID_REMOVE, "&Remove a checkbox\tCtrl--");

    wxMenuBar *menu_bar = new wxMenuBar;
    menu_bar->Append(menu, "&Wrap sizer");

    SetMenuBar( menu_bar );

    wxBoxSizer *root = new wxBoxSizer( wxVERTICAL );

    wxStaticBoxSizer *topSizer = new wxStaticBoxSizer( wxVERTICAL, this, "Wrapping check-boxes" );
    m_checkboxParent = topSizer->GetStaticBox();
    m_wrapSizer = new wxWrapSizer(wxHORIZONTAL);

    // A number of checkboxes inside a wrap sizer
    for( int i = 0; i < 6; i++ )
        DoAddCheckbox();

    topSizer->Add( m_wrapSizer, wxSizerFlags(1).Expand());
    root->Add( topSizer, wxSizerFlags().Expand().Border());

    // A shaped item inside a box sizer
    wxSizer *bottomSizer = new wxStaticBoxSizer( wxVERTICAL, this, "With wxSHAPED item" );
    wxSizer *horzBoxSizer = new wxBoxSizer(wxHORIZONTAL);
    bottomSizer->Add( horzBoxSizer, 100, wxEXPAND );
    horzBoxSizer->Add( new wxListBox(this,wxID_ANY,wxPoint(0,0),wxSize(70,70)), 0, wxEXPAND|wxSHAPED );
    horzBoxSizer->Add( 10,10 );
    horzBoxSizer->Add( new wxCheckBox(this,wxID_ANY,"A much longer option..."), 100, 0, 5 );

    root->Add( bottomSizer, 1, wxEXPAND | wxALL, 5 );

    // Set sizer for window
    SetSizerAndFit( root );
}

void MyWrapSizerFrame::DoAddCheckbox()
{
    m_wrapSizer->Add(new wxCheckBox
                         (
                            m_checkboxParent,
                            wxID_ANY,
                            wxString::Format
                            (
                                "Option %d",
                                (int)m_wrapSizer->GetItemCount() + 1
                            )
                         ),
                     wxSizerFlags().Centre().Border());
}

void MyWrapSizerFrame::OnAddCheckbox(wxCommandEvent& WXUNUSED(event))
{
    DoAddCheckbox();
    Layout();
}

void MyWrapSizerFrame::OnRemoveCheckbox(wxCommandEvent& WXUNUSED(event))
{
    if ( m_wrapSizer->IsEmpty() )
    {
        wxLogStatus(this, "No more checkboxes to remove.");
        return;
    }

    delete m_wrapSizer->GetItem(m_wrapSizer->GetItemCount() - 1)->GetWindow();
    Layout();
}

IMPLEMENT_APP(MyApp)
#else

#endif


void freedom_for_memory(void)
{
    free_filez();

    if (dnaspace) free(dnaspace);
    dnaspace = (char *)0;
    if (dnacnts) free(dnacnts);
    dnacnts = (unsigned short int *)0;
    if (dnamms) free(dnamms);
    dnamms = (unsigned short int *)0;
    if (ppp) free(ppp);
    ppp = (unsigned char *)0;
}


