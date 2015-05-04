 
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

void ImageFilledRectangle(struct image_type *im,int x, int y, int x1, int y1 ,int color);
void ImageString(struct image_type *im, int x, int y, unsigned char *s, int color);
int do_by_gene_name_from_refflat(const char gene[],char chr[],int *start,int *end);

unsigned int ALVIEWWRAP_lodepng_encode_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);
unsigned int lodepng_encode24_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);

// horror  !!!  - can't get it to eat unsigned this prototype: lodepng_encode24(unsigned char** out, size_t* outsize, const unsigned char* image, unsigned w, unsigned h)

