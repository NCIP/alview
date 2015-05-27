 
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

void set_ucsc_url(char *url); // uses blds chr, start , end to make linkt to ucsc browser 
void set_blat_url(char *url) ;
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

