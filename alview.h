 
#ifdef ALVIEW_H  
#else
#define ALVIEW_H  
struct image_type
{
    unsigned char *data;   // RGB space (3 * height * widht)
    int height;
    int width;
};
#endif

unsigned ImageSaveAsPNG(struct image_type *im, char *fn);
void ImageFilledRectangle(struct image_type *im,int x, int y, int x1, int y1 ,int color);
void ImageString(struct image_type *im, int x, int y, unsigned char *s, int color);
void ImageString5x5(struct image_type *im, int x, int y, unsigned char *s, int color);

int do_by_gene_name_from_refflat(const char gene[],char chr[],int *start,int *end);

void set_ucsc_url(char *url);  // uses blds chr, start , end to make linkt to ucsc browser 
void set_blat_url(char *url);
void set_alviewhelp_url(char *url); 

void open_url(char *url); // invoke URL using your system's default(?) browser

#ifdef UNIX
void linux_call_url(char *url); // calls gtk_link_button_new() stuff
#endif
#ifdef MAC
void cocoa_call_url(char *url);  // [[NSWorkspace sharedWorkspace] openFile:mad:"http://www.mylink.com"];
#endif
#ifdef WIN32
void windows_call_url(char *url); 
#endif

unsigned int ALVIEWWRAP_lodepng_encode_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);
unsigned int lodepng_encode24_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);
void gtk_set_atcg_info_label(char s[]);

void init_colors(void);


#define dna_a_R "9f" 
#define dna_a_G "00" 
#define dna_a_B "5f" 
#define dna_c_R "ff" 
#define dna_c_G "5f" 
#define dna_c_B "00" 
#define dna_g_R "bf" 
#define dna_g_G "ff" 
#define dna_g_B "00" 
#define dna_t_R "00" 
#define dna_t_G "3f" 
#define dna_t_B "bf" 
#define dna_I_R "00" 
#define dna_I_G "ff" 
#define dna_I_B "ff" 
#define dna_D_R "ff" 
#define dna_D_G "b0" 
#define dna_D_B "00" 

#define Hdna_a_R 0x9f 
#define Hdna_a_G 0x00 
#define Hdna_a_B 0x5f 
#define Hdna_c_R 0xff 
#define Hdna_c_G 0x5f 
#define Hdna_c_B 0x00 
#define Hdna_g_R 0xbf 
#define Hdna_g_G 0xff 
#define Hdna_g_B 0x00 
#define Hdna_t_R 0x00 
#define Hdna_t_G 0x3f 
#define Hdna_t_B 0xbf 
#define Hdna_I_R 0x00 
#define Hdna_I_G 0xff 
#define Hdna_I_B 0xff 
#define Hdna_D_R 0xff 
#define Hdna_D_G 0xb0 
#define Hdna_D_B 0x00 
