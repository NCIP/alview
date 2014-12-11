
struct image_type
{
    unsigned char *data; // RGB space (3 * height * widht)
    int height;
    int width;
};

void ImageFilledRectangle(struct image_type *im,int x, int y, int x1, int y1 ,int color); 
void ImageString(struct image_type *im, int x, int y, unsigned char *s, int color);

