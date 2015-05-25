
/*
vi +164 nsfp2mydb.c ;gcc -Wall -O3 -o nsfp2mydb nsfp2mydb.c
*/

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "errno.h"

/*
12 rs730880256
*/
// new way
//

struct snp_type   // ****** input and output as packed in this order !!!! >>> 
{
    unsigned int s;  // start 4 bytes 
    unsigned int e;  // end 4 bytes 
    unsigned char chrindex; // chr index
    char mask; // info on sift , bit 1 = is sift , bit2 = nsfp origin
    char rs[12]; //  null term 
};

#if 0
// old way 
#define MAX_CHR_NAME 16
struct snp_type 
{
    char chr[MAX_CHR_NAME];
    unsigned int s; // 4 bytes 
    unsigned int e; // 4 bytes 
    float sift;
    char rs[14] ; // "rs539566132" // filter out "rs200354637;rs2257402"
};
#endif


int numoutz = 0;
#define MAXREC 80000000
struct snp_type outz[MAXREC];

int snpcmp(const void *aarg, const void *barg)
{
     struct snp_type *a;
     struct snp_type *b;

     a = (struct snp_type *)aarg;
     b = (struct snp_type *)barg;

     if (a->chrindex < b->chrindex ) return -1;
     else if (a->chrindex > b->chrindex ) return 1;

     if (a->s <b->s ) return -1;
     else if (a->s >b->s ) return 1;
     return 0;
}

char prevs[16] = "";
unsigned int previ;
static unsigned char chr2index(char *s)
{
    unsigned char ret = 0xff;
    if (strcmp(s,prevs) == 0) return previ;

    
    if (strcmp(s,"1") == 0) { ret =  0; }  
    else if (strcmp(s,"2") == 0) { ret =  1;  } 
    else if (strcmp(s,"3") == 0) { ret =  2;  } 
    else if (strcmp(s,"4") == 0) { ret =  3;  } 
    else if (strcmp(s,"5") == 0) { ret =   4;  } 
    else if (strcmp(s,"6") == 0) { ret =   5;  } 
    else if (strcmp(s,"7") == 0) { ret =   6;  } 
    else if (strcmp(s,"8") == 0) { ret =   7;  } 
    else if (strcmp(s,"9") == 0) { ret =   8;  } 
    else if (strcmp(s,"10") == 0) { ret =   9;  } 
    else if (strcmp(s,"11") == 0) { ret =   10;  } 
    else if (strcmp(s,"12") == 0) { ret =   11;  } 
    else if (strcmp(s,"13") == 0) { ret =   12;  } 
    else if (strcmp(s,"14") == 0) { ret =   13;  } 
    else if (strcmp(s,"15") == 0) { ret =   14;  } 
    else if (strcmp(s,"16") == 0) { ret =   15;  } 
    else if (strcmp(s,"17") == 0) { ret =   16;  } 
    else if (strcmp(s,"18") == 0) { ret =   17;  } 
    else if (strcmp(s,"19") == 0) { ret =   18;  } 
    else if (strcmp(s,"20") == 0) { ret =   19;  } 
    else if (strcmp(s,"21") == 0) { ret =   20;  } 
    else if (strcmp(s,"22") == 0) { ret =   21;  } 
    else if (strcmp(s,"X") == 0) { ret =   22;  } 
    else if (strcmp(s,"Y") == 0) { ret =   23;  } 
    else if (strcmp(s,"M") == 0) { ret =   24;  } 
    else if (strcmp(s,"chr1") == 0) { ret =   0;  } 
    else if (strcmp(s,"chr2") == 0) { ret =   1;  } 
    else if (strcmp(s,"chr3") == 0) { ret =   2;  } 
    else if (strcmp(s,"chr4") == 0) { ret =   3;  } 
    else if (strcmp(s,"chr5") == 0) { ret =   4;  } 
    else if (strcmp(s,"chr6") == 0) { ret =   5;  } 
    else if (strcmp(s,"chr7") == 0) { ret =   6;  } 
    else if (strcmp(s,"chr8") == 0) { ret =   7;  } 
    else if (strcmp(s,"chr9") == 0) { ret =   8;  } 
    else if (strcmp(s,"chr10") == 0) { ret =   9;  } 
    else if (strcmp(s,"chr11") == 0) { ret =   10;  } 
    else if (strcmp(s,"chr12") == 0) { ret =   11;  } 
    else if (strcmp(s,"chr13") == 0) { ret =   12;  } 
    else if (strcmp(s,"chr14") == 0) { ret =   13;  } 
    else if (strcmp(s,"chr15") == 0) { ret =   14;  } 
    else if (strcmp(s,"chr16") == 0) { ret =   15;  } 
    else if (strcmp(s,"chr17") == 0) { ret =   16;  } 
    else if (strcmp(s,"chr18") == 0) { ret =   17;  } 
    else if (strcmp(s,"chr19") == 0) { ret =   18;  } 
    else if (strcmp(s,"chr20") == 0) { ret =   19;  } 
    else if (strcmp(s,"chr21") == 0) { ret =   20;  } 
    else if (strcmp(s,"chr22") == 0) { ret =   21;  } 
    else if (strcmp(s,"chrX") == 0) { ret =   22;  } 
    else if (strcmp(s,"chrY") == 0) { ret =   23;  } 
    else if (strcmp(s,"chrM") == 0) { ret =   24;  } 
    else 
    {
         fprintf(stderr,"ERROR - invalid chromosome [%s]\n",s);  
         prevs[0] = (char)0;
         return 0xff;
    }

    strcpy(prevs,s);
    previ = ret;
    return ret;
}


// ***** pass in zero based, not one based genomic locations
void addrec(unsigned int chrindex_arg, int s,  int e,char *rs, double sift, int force)
{
    unsigned int s_ui;
    unsigned int e_ui;

    if ((force == 0) && (sift > 0.05) ) return;
    if (numoutz >= MAXREC)
    {
        fprintf(stderr,"ERROR in addrec() overflow %d \n",numoutz); 
        exit(0);
    }
    s_ui = s;
    e_ui = e;
    if (!force) 
    {
       if (numoutz)           // not the first one
       {
           if ((s_ui == outz[numoutz-1].s) && ((unsigned char)(chrindex_arg && 0xff) == outz[numoutz-1].chrindex) )      // if same snp
           {
               if ((sift >= 0.0) && (sift <= 0.05) ) 
                   outz[numoutz-1].mask &= 0x01;
               return;
           }
       }
    }

    outz[numoutz].chrindex = chrindex_arg;
    outz[numoutz].s = s_ui;
    outz[numoutz].e = e_ui;
    strcpy(outz[numoutz].rs,rs);
    outz[numoutz].mask = 0;
    if ((sift >= 0.0) && (sift <= 0.05) ) 
    {
        outz[numoutz].mask |= 0x01; // set low sift bit 
    }
    if (!force)
    {
        outz[numoutz].mask |= 0x02; // set 2nd bit for force
    }
    numoutz++;
    return;
}


// test routine
void output_to_text()
{
    int stat = 0;
    int cnt;
    FILE *fp;
    struct snp_type orec;
/* example siftscore values :  .  0.0 0.001 SIFT_score - last "SIFT_score is header line only */

    fp = fopen("out.hg19","r");
    cnt = 0;
    while (1)
    {
        stat = fread(&orec,sizeof(struct snp_type),1,fp);
        if (stat < 1) break;
// printf("%u %u %u %u %s\n", orec.chrindex, orec.s, orec.e, orec.mask, orec.rs);
        cnt++;
    }
    fclose(fp);
    fprintf(stderr,"dumped %d recs \n",cnt);
}


void output_to_binary(char *outfn)
{
    int fixes = 0;
    int outcnt = 0;
    int i,j;
    FILE *fp;

    fp = fopen(outfn,"w");
    for (i=0;i<numoutz;) 
    {
        for (j=i+1;j<numoutz;) // look ahead 
        {
             if (outz[i].chrindex != outz[j].chrindex) break; // different chromsomes
             if (outz[i].s != outz[j].s) break;   // different start , break;
             fixes++;
             if (outz[j].mask & 0x01) // if bit set
                outz[i].mask |= 0x01; // set sift bit
             j++;
        }
        fwrite(&outz[i].s,sizeof(unsigned int ),1,fp); // 4 
        fwrite(&outz[i].e,sizeof(unsigned int ),1,fp); // 4 
        fwrite(&outz[i].chrindex,1,1,fp); // 1
        fwrite(&outz[i].mask,1,1,fp); // 1 
        fwrite(&outz[i].rs,12,1,fp); // 12 
        outcnt++;
printf("%u %u %u %u %s\n",outz[i].chrindex,outz[i].s,outz[i].e,outz[i].mask,outz[i].rs);
        i = j;
    }
    fclose(fp);

fprintf(stderr,"end output_to_binary(), numoutz=%d, fixes=%d outcnt=%d \n", numoutz,fixes,outcnt);  

}


char *filez[] = 
{
"dbNSFP3.0b2a_variant.chr1" , 
"dbNSFP3.0b2a_variant.chr2" , 
"dbNSFP3.0b2a_variant.chr3" , 
"dbNSFP3.0b2a_variant.chr4" , 
"dbNSFP3.0b2a_variant.chr5" , 
"dbNSFP3.0b2a_variant.chr6" , 
"dbNSFP3.0b2a_variant.chr7" , 
"dbNSFP3.0b2a_variant.chr8" , 
"dbNSFP3.0b2a_variant.chr9" , 
"dbNSFP3.0b2a_variant.chr10" , 
"dbNSFP3.0b2a_variant.chr11" , 
"dbNSFP3.0b2a_variant.chr12" , 
"dbNSFP3.0b2a_variant.chr13" , 
"dbNSFP3.0b2a_variant.chr14" , 
"dbNSFP3.0b2a_variant.chr15" , 
"dbNSFP3.0b2a_variant.chr16" , 
"dbNSFP3.0b2a_variant.chr17" , 
"dbNSFP3.0b2a_variant.chr18" , 
"dbNSFP3.0b2a_variant.chr19" , 
"dbNSFP3.0b2a_variant.chr20" , 
"dbNSFP3.0b2a_variant.chr21" , 
"dbNSFP3.0b2a_variant.chr22" , 
"dbNSFP3.0b2a_variant.chrM" , 
"dbNSFP3.0b2a_variant.chrX" , 
"dbNSFP3.0b2a_variant.chrY" , 
(char *)0  
};


char buff[62][6]; // max of 60 sifts,  strlen is max of 5

int drive_nsfp(int bld)
{
    unsigned char uc;
    FILE *fp;
    double d;
    double mmin;
    int sift_rec_cnt;
    int lineno = 0;
    int i;
    int j;
    char s[5120];
    char t[5120];
    char cmd[5120];
    char chr[512]; // 1	
    char pos[512];	// 2	(1-based)
    char rs_dbSNP142[512];// 	7	
    char hg19_chr[512]; //8	
    char hg19_pos[512]; // 9	(1-based)
    char hg18_chr[512]; // 10	
    char hg18_pos[512]; // 11	(1-based)
    char SIFT_score[512]; // 24	

    
#define DO16 0
#if DO16 
    for (i=15; i<16 /* filez[i] */ ;i++) 
#else
    for (i=0; filez[i] ;i++) 
#endif
    {
        sprintf(cmd,"cat %s | cut -f1,2,7,8,9,10,11,24",filez[i]);
        fp = popen(cmd,"r");
        lineno = 0;
        while (fgets(s,5000,fp))
        {
            if (lineno++ == 0) continue; // skip first line
            sscanf(s, "%s %s %s %s %s %s %s %s", 
                   chr, pos, rs_dbSNP142, hg19_chr, hg19_pos, hg18_chr, hg18_pos, SIFT_score);
            for (j=0;SIFT_score[j]; j++) 
                if (SIFT_score[j]  == ';') SIFT_score[j]= ' '; 
//            memset(buff,0,sizeof(buff));
            sift_rec_cnt = sscanf(SIFT_score,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", 
&buff[0][0],&buff[1][0],&buff[2][0],&buff[3][0],&buff[4][0],&buff[5][0],&buff[6][0],&buff[7][0],&buff[8][0],&buff[9][0],
&buff[10][0],&buff[11][0],&buff[12][0],&buff[13][0],&buff[14][0],&buff[15][0],&buff[16][0],&buff[17][0],&buff[18][0],&buff[19][0],
&buff[20][0],&buff[21][0],&buff[22][0],&buff[23][0],&buff[24][0],&buff[25][0],&buff[26][0],&buff[27][0],&buff[28][0],&buff[29][0],
&buff[30][0],&buff[31][0],&buff[32][0],&buff[33][0],&buff[34][0],&buff[35][0],&buff[36][0],&buff[37][0],&buff[38][0],&buff[39][0],
&buff[40][0],&buff[41][0],&buff[42][0],&buff[43][0],&buff[44][0],&buff[45][0],&buff[46][0],&buff[47][0],&buff[48][0],&buff[49][0],
&buff[50][0],&buff[51][0],&buff[52][0],&buff[53][0],&buff[54][0],&buff[55][0],&buff[56][0],&buff[57][0],&buff[58][0],&buff[59][0]
                  );
              mmin = 1.0;
              d = -1.0;
              for (j=0;j<sift_rec_cnt;j++)
              {
                  if (buff[j][0] == '.') continue; // skip this
                  d = -1.0;
                  d = atof(&buff[j][0]);
                  if (d < mmin) mmin = d;
              }
              for (j=0;rs_dbSNP142[j]; j++) 
                 if (rs_dbSNP142[j]  == ';') rs_dbSNP142[j]= (char)0; 
              if ((rs_dbSNP142[0] != '.') || ((mmin>=0.0) && (mmin < 0.05)) )
              {
                  if (bld == 19) 
                  {
                      if (hg19_chr[0] == '.') continue;
                      sprintf(t,"chr%s",hg19_chr);
                      uc = chr2index(hg19_chr); 
                      if (uc != (unsigned char)0xff) 
                          addrec(uc, atoi(hg19_pos)-1, atoi(hg19_pos),rs_dbSNP142, d,0);
// subract one because one based
                  }
                  else
                  {
                      if (hg18_chr[0] == '.') continue;
                      sprintf(t,"chr%s",hg18_chr);
                      uc = chr2index(hg18_chr); 
                      if (uc != (unsigned char)0xff) 
                          addrec(uc, atoi(hg18_pos)-1, atoi(hg18_pos),rs_dbSNP142, d,0);
// subract one because one based
                  }
              }
        }
        if (fp) pclose(fp);
        fp = (FILE *)0;
    }
fprintf(stderr,"numoutz = %d , sizeof struct snp_type=%ld\n",numoutz,sizeof(struct snp_type )); 
    return 0;   
}


 // dbsnp file is from ucsc
int drive_dbsnp(char *fn)
{
    unsigned char uc;
    FILE *fp;
    double d;
    int lineno = 0;
    char s[5120];
    char cmd[5120];
    char chr[512]; // 2	
    char pos[512];	// 3	(1-based)
    char pose[512];	// 4	(1-based)
    char rs[512];// 	5	


#if DO16
    sprintf(cmd,"grep chr16 %s | sed 's/chr//' | cut -f2,3,4,5",fn); // chr s e rs
#else
    sprintf(cmd,"cat %s | cut -f2,3,4,5",fn); // chr s e rs
#endif
// xxx 
    fp = popen(cmd,"r");
    lineno = 0;
    d = -1.0;
    while (fgets(s,5000,fp))
    {
        sscanf(s,"%s %s %s %s", chr, pos, pose , rs);
        if (atoi(pos)+1 == atoi(pose)) // if single base snp
        {
            uc = chr2index(chr);
            if (uc != (unsigned char)0xff) 
                addrec(uc, atoi(pos), atoi(pose),rs, d,1);
        }
        lineno++;
    }
    if (fp) pclose(fp);
    fp = (FILE *)0;
fprintf(stderr,"lineno = %d numoutz = %d ,  sizeof struct snp_type=%ld\n",
              lineno,numoutz,sizeof(struct snp_type )); 
    return 0;   
}


int main()
{
// printf("size of unsigned int = %ld \n",sizeof(unsigned int)); exit(0);
fprintf(stderr,"main 1 , drive_nsfp\n"); fflush(stderr); 
    numoutz = 0;
    drive_nsfp(19);
fprintf(stderr,"main 2 drive dbsnp\n"); fflush(stderr); 
    drive_dbsnp("snp142Common.txt");
fprintf(stderr,"main 3 qsort \n"); fflush(stderr); 
    qsort(outz,numoutz,sizeof(struct snp_type),snpcmp); 
fprintf(stderr,"main 4 uotuput \n"); fflush(stderr); 
    output_to_binary("hg19.snps.dat");

fprintf(stderr,"main 5\n"); fflush(stderr); 

    memset(outz,0,sizeof(outz));
    numoutz = 0;
fprintf(stderr,"main 6, drive ns 18 \n"); fflush(stderr); 
    drive_nsfp(18);
fprintf(stderr,"main 7, drive dbsnp 18 \n"); fflush(stderr); 
    drive_dbsnp("hg18snp130.txt");
fprintf(stderr,"main 8, drive before qsort 18 \n"); fflush(stderr); 

    qsort(outz,numoutz,sizeof(struct snp_type),snpcmp); 
fprintf(stderr,"main 9, output hg18 \n"); fflush(stderr); 
    output_to_binary("hg18.snps.dat"); fprintf(stderr,"done \n"); fflush(stderr); 

    return 0;   
}


/*
     1	#chr
     2	pos(1-based)
     3	ref
     4	alt
     5	aaref
     6	aaalt
     7	rs_dbSNP142
     8	hg19_chr
     9	hg19_pos(1-based)
    10	hg18_chr
    11	hg18_pos(1-based)
    12	genename
    13	cds_strand
    14	refcodon
    15	codonpos
    16	codon_degeneracy
    17	Ancestral_allele
    18	AltaiNeandertal
    19	Denisova
    20	Ensembl_geneid
    21	Ensembl_transcriptid
    22	Ensembl_proteinid
    23	aapos
    24	SIFT_score
    25	SIFT_converted_rankscore
    26	SIFT_pred
    27	Uniprot_acc_Polyphen2
    28	Uniprot_id_Polyphen2
    29	Uniprot_aapos_Polyphen2
    30	Polyphen2_HDIV_score
    31	Polyphen2_HDIV_rankscore
    32	Polyphen2_HDIV_pred
    33	Polyphen2_HVAR_score
    34	Polyphen2_HVAR_rankscore
    35	Polyphen2_HVAR_pred
    36	LRT_score
    37	LRT_converted_rankscore
    38	LRT_pred
    39	LRT_Omega
    40	MutationTaster_score
    41	MutationTaster_converted_rankscore
    42	MutationTaster_pred
    43	MutationTaster_model
    44	MutationTaster_AAE
    45	Uniprot_id_MutationAssessor
    46	Uniprot_variant_MutationAssessor
    47	MutationAssessor_score
    48	MutationAssessor_rankscore
    49	MutationAssessor_pred
    50	FATHMM_score
    51	FATHMM_converted_rankscore
    52	FATHMM_pred
    53	PROVEAN_score
    54	PROVEAN_converted_rankscore
    55	PROVEAN_pred
    56	Transcript_id_VEST3
    57	Transcript_var_VEST3
    58	VEST3_score
    59	VEST3_rankscore
    60	CADD_raw
    61	CADD_raw_rankscore
    62	CADD_phred
    63	MetaSVM_score
    64	MetaSVM_rankscore
    65	MetaSVM_pred
    66	MetaLR_score
    67	MetaLR_rankscore
    68	MetaLR_pred
    69	Reliability_index
    70	GERP++_NR
    71	GERP++_RS
    72	GERP++_RS_rankscore
    73	phyloP7way_vertebrate
    74	phyloP7way_vertebrate_rankscore
    75	phastCons7way_vertebrate
    76	phastCons7way_vertebrate_rankscore
    77	SiPhy_29way_pi
    78	SiPhy_29way_logOdds
    79	SiPhy_29way_logOdds_rankscore
    80	1000Gp3_AC
    81	1000Gp3_AF
    82	1000Gp3_AFR_AC
    83	1000Gp3_AFR_AF
    84	1000Gp3_EUR_AC
    85	1000Gp3_EUR_AF
    86	1000Gp3_AMR_AC
    87	1000Gp3_AMR_AF
    88	1000Gp3_EAS_AC
    89	1000Gp3_EAS_AF
    90	1000Gp3_SAS_AC
    91	1000Gp3_SAS_AF
    92	TWINSUK_AC
    93	TWINSUK_AF
    94	ALSPAC_AC
    95	ALSPAC_AF
    96	ESP6500_AA_AC
    97	ESP6500_AA_AF
    98	ESP6500_EA_AC
    99	ESP6500_EA_AF
   100	ExAC_AC
   101	ExAC_AF
   102	ExAC_Adj_AC
   103	ExAC_Adj_AF
   104	ExAC_AFR_AC
   105	ExAC_AFR_AF
   106	ExAC_AMR_AC
   107	ExAC_AMR_AF
   108	ExAC_EAS_AC
   109	ExAC_EAS_AF
   110	ExAC_FIN_AC
   111	ExAC_FIN_AF
   112	ExAC_NFE_AC
   113	ExAC_NFE_AF
   114	ExAC_SAS_AC
   115	ExAC_SAS_AF
   116	clinvar_rs
   117	clinvar_clnsig
   118	clinvar_trait
   119	Interpro_domain
*/

