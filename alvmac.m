
#import <Cocoa/Cocoa.h>
#include "stdio.h"

#define DEFAULT_WIDTH 1000
#define DEFAULT_HEIGHT 400

#include "alview.h"

void get_params_and_draw(int arg);
void calc_and_draw_new_img(void);
void drawclip(void);
void alview_load_config(void);
int end_select_x = 0;
int start_select_x = 0;
unsigned char *imgen_mem(char fn[], char chr[],int khrstart, int khrend, int h, int w,int *ret_status) ;
int xoron = 0;
void freedom_for_memory(void);
void alview_init(void);
int parse_position(const char argposition[],char chr[],int *start,int *end);

extern char fn_bam[];    // size=FILENAME_MAX -the bam file name
extern struct image_type *im;
extern char GENOMEDATADIR[512];

unsigned char *img;
unsigned char *img2;
int diw;     // "darea" pixbuff width 
int dih;     // "darea" pixbuff height 
int darea_on = 0; // flag, means RGB image is loaded in "im"
int khrstart =0 ;
int khrend = 0 ;
char khr[1024];
char pos[1024]; // genomic position or gene name


// --- IMAGE view ENTRY
@interface MyImageViewClass : NSImageView 
{
@public
NSPoint lastDragLocation;
NSPoint startDragLocation;
}
@end
@implementation MyImageViewClass
- (void)mouseUp:(NSEvent *)event {
printf("mouse up"); fflush(stdout); 
NSPoint newDragLocation = [event locationInWindow];
end_select_x = newDragLocation.x;
// printf("drag location %f %f\n",newDragLocation.x, newDragLocation.y);
calc_and_draw_new_img();
//      if ([[self target] respondsToSelector:[self action]]) { [NSApp sendAction:[self action] to:[self target] from:self]; }
}

- (void)mouseDown:(NSEvent *)event 
{
    //  printf("mouse down"); fflush(stdout); 
    startDragLocation = [event locationInWindow];
    // printf("start location %f %f \n",startDragLocation.x, startDragLocation.y);
    lastDragLocation = [event locationInWindow];
    start_select_x = startDragLocation.x;
//if ([[self target] respondsToSelector:[self action]]) { [NSApp sendAction:[self action] to:[self target] from:self]; }
}

- (void)mouseDragged:(NSEvent *)event 
{
NSPoint newDragLocation = [event locationInWindow];
// printf("drag location %f %f \n",newDragLocation.x, newDragLocation.y);
   end_select_x = newDragLocation.x;
   drawclip();

#if 0
NSPoint thisOrigin = [self frame].origin;
thisOrigin.x += (-lastDragLocation.x + newDragLocation.x);
thisOrigin.y += (-lastDragLocation.y + newDragLocation.y);
[self setFrameOrigin:thisOrigin];
lastDragLocation = newDragLocation;
#endif
}
@end


// --- TEXT ENTRY
@interface MyTextClass:NSTextField
{
@public
int myid;
}
@end

MyImageViewClass *theImageView;
NSImage *image1;
NSImage *image2; // the being selected image
NSRect superframe;

@implementation MyTextClass
-(void)keyUp:(NSEvent *)event   // can't override keydown NSTextField - so override keyUp
{
    if ([event keyCode] == 36) // enter key pressed 
    {
        NSString *userinput = [self stringValue];
        char *sptr = (char *)[userinput UTF8String];
printf("got [%s]\n",sptr);
        get_params_and_draw(1);
// printf("keydown , event keycode is [%d]  \n",[event keyCode]);
    }
    [super keyUp:event];
}
- (BOOL)control:(NSControl *)control textView:(NSTextView *)fieldEditor doCommandBySelector:(SEL)commandSelector
{           
    BOOL retval = NO;
    if (commandSelector == @selector(insertNewline:)) {
        retval = YES; // causes Apple to NOT fire the default enter action
        // Do your special handling of the "enter" key here
    }
//    [super control:event];
    // NSLog(@"Selector = %@", NSStringFromSelector( commandSelector ) );
    return retval;  
}
@end

struct my_text_type
{
    int id; 
    int x; int y; int w; int h; // label location
    int x2; int y2; int w2; int h2; // entry box 
    char *label;
    char *default_value;
    MyTextClass *object;
};

struct my_text_type texts[] =
{       // xxx
 { 100  ,  10 , 65 , 70, 20,      70 , 65 , 170, 20 , "Position:" ,    "chr1:11778-15130"} ,
 { 101  , 320 , 65 , 50, 20,     362 , 65 , 50,  20 , "Width:"  ,      "900" } ,
 { 112  , 422 , 65 , 50, 20,     470 , 65 , 50,  20 , "Height:"  ,     "400" } ,
 { -1 ,   -1 , -1,  -1, -1,      -1 , -1 , -1 , -1 , "ERRoverflow"  , "(enter position or gene)"} 
};
NSImage *rbg2NSImage(unsigned char *data, int width, int height)
{
    size_t bufferLength = width * height * 3;
    CGDataProviderRef provider = CGDataProviderCreateWithData(NULL, data, bufferLength, NULL);
    size_t bitsPerComponent = 8;
    size_t bitsPerPixel = 24;
    size_t bytesPerRow = 3 * width;
    CGColorSpaceRef colorSpaceRef = CGColorSpaceCreateDeviceRGB();
    CGBitmapInfo bitmapInfo = kCGBitmapByteOrderDefault ; //  | kCGImageAlphaPremultipliedLast;
    CGColorRenderingIntent renderingIntent = kCGRenderingIntentDefault;
    CGImageRef iref = CGImageCreate(width, 
                                height, 
                                bitsPerComponent, 
                                bitsPerPixel, 
                                bytesPerRow, 
                                colorSpaceRef, 
                                bitmapInfo, 
                                provider,   // data provider
                                NULL,       // decode
                                YES,        // should interpolate
                                renderingIntent);
    return [[NSImage alloc] initWithCGImage:iref size:NSMakeSize(width, height)];
}

void sanity_khr()
{
    if (khrstart < 1) khrstart = 1;
    if (khrend < 2) khrend = 2;
    if (khrend <= khrstart) khrend = khrstart + 2;
}

void draw_mac(unsigned char *data, int width, int height,int arg)
{
   if (!data)
   {
printf("draw_mac data is null - arg = %d\n",arg); 
return;
   }
   [theImageView setImage:image1];
   [image1 release];

   image1 = rbg2NSImage(data, width, height);
   [theImageView setImage:image1];

// resize ... 
    NSRect imageRect = NSMakeRect(4, superframe.size.height-150-[image1 size].height,
         [image1 size].width, [image1 size].height);
   theImageView.frame = imageRect;
#if 0
printf("xframe %d %d %d %d \n", 
   (int) theImageView.frame.origin.x , (int) theImageView.frame.origin.y , 
   (int) theImageView.frame.size.width , (int) theImageView.frame.size.height);  fflush(stdout); 
printf("zframe %d %d %d %d \n", 
   (int) theImageView.frame.origin.x ,   (int) theImageView.frame.origin.y , 
   (int) theImageView.frame.size.width , (int) theImageView.frame.size.height);  fflush(stdout); 

   NSRect imageRect = NSMakeRect(width,height,[image1 size].width, [image1 size].height);
   theImageView.frame = NSMakeRect(4, superframe.size.height - 320, 
            theImageView.image.size.width, theImageView.image.size.height);
#endif
// ??? Adjust the size of your scroll views view accordingly.
   [theImageView setNeedsDisplay:true] ;
}

void do_img_and_draw(char *fn, char khr[], int s, int e, int arg)
{
    int status;

#if 1
printf("in do_img_and_draw() file name = %p \n", fn);  fflush(stdout); 
#endif
    if (fn) strcpy(fn_bam,(char *)fn); 
// sanity 
    if (dih < 0) dih = DEFAULT_HEIGHT;
    if (diw < 0) diw = DEFAULT_WIDTH;
    if (khr[0] == (char)0)
    {
         strcpy(khr,"chr1"); 
         khrstart = 12982;
         khrend = 15849;
    }
    if (img) { free(img); }
#if 1
printf("in do_img_and_draw() before imgen_mem %s %d %d %d %d\n",khr,khrstart,khrend,dih,diw);  fflush(stdout); 
#endif

    img = (void *)imgen_mem(fn_bam,khr,khrstart,khrend,dih, diw,&status);
printf("in do_img_and_draw after imgen_mem() img =%p, status=%d \n",img,status); fflush(stdout); 
    if (img)
    {
        darea_on = 1;
        draw_mac(img,diw,dih,1);
#if 0
printf("lode after imgen_mem() img=%p\n",img); fflush(stdout); 
unsigned lodepng_encode24_file(const char* filename, const unsigned char* image, unsigned w, unsigned h);
        lodepng_encode24_file("tmp.png", img, diw, dih);
printf("after lodepng_encode24_file\n"); fflush(stdout); 
printf("after lodepng w=%d h=%d\n",diw,dih); fflush(stdout); 
#endif
        return;
    }
printf("ERROR: imgdata is null , no lode\n"); fflush(stdout); 
    return;
}

void get_params_and_draw(int arg)
{
    int status;
    char *sptr;

    NSString *userinput = [texts[0].object stringValue];
    sptr = (char *)[userinput UTF8String];
printf("xxx in get_params_and_draw() - sptr=[%s]\n",sptr); 
    strcpy(pos,sptr);
    NSString *userinput2 = [texts[1].object stringValue];
    sptr = (char *)[userinput2 UTF8String];
    diw = atoi(sptr);
    NSString *userinput3 = [texts[2].object stringValue];
    sptr = (char *)[userinput3 UTF8String];
    dih = atoi(sptr);

printf("xxx in get_params_and_draw() before parse_position - pos=[%s]\n",pos); 
    status = parse_position(pos,khr,&khrstart,&khrend); 
printf("xxx in get_params_and_draw() after parse_position - pos=[%s]\n",pos); 

// put pos back to screen
    NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
    [texts[0].object setStringValue:defval];
    do_img_and_draw(fn_bam, khr,khrstart,khrend,arg);
}

void calc_and_draw_new_img(void)
{
    char *sptr;
    int status;
    int size;
    int tmp;

    xoron = 0;
    if (start_select_x > end_select_x)
    {
        tmp = start_select_x; 
        start_select_x =  end_select_x; 
        end_select_x = tmp; 
    }

    tmp = khrstart;
    khrstart = (int)(khrstart + (double)((start_select_x * ((khrend - khrstart)) / (double)diw)));
    khrend =   (int)(tmp + (double)((end_select_x   * ((khrend - tmp)) / (double)diw)));             
// fprintf(stderr,"in release , END %d %d \n",khrstart,khrend); fflush(stderr);  

    if (khrstart < 1) khrstart = 1;
    if (khrend < 2) khrend = 2;
    if (khrend <= khrstart) khrend = khrstart + 2;
// set pos
    sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
    start_select_x = end_select_x = 1; // reset

// get pos width heigh ...
    NSString *userinput = [texts[0].object stringValue];
    sptr = (char *)[userinput UTF8String];
    strcpy(pos,sptr);
    NSString *userinput2 = [texts[1].object stringValue];
    sptr = (char *)[userinput2 UTF8String];
    diw = atoi(sptr);
    NSString *userinput3 = [texts[2].object stringValue];
    sptr = (char *)[userinput3 UTF8String];
    dih = atoi(sptr);
    sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
    status = parse_position(pos,khr,&khrstart,&khrend); 
// put pos back to screen
    NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
    [texts[0].object setStringValue:defval];
    do_img_and_draw(fn_bam, khr,khrstart,khrend,2);
}



int on_button_move(int movek)
{
    int page_nts = 0;
    int status = 0;

    if (fn_bam[0] == (char)0)
    {
        printf("No bam specified.  Plase pick one. \n");
        return -1;
    }

printf("in on_button_moveSTART movek=%d , khrstart = %d khrend = %d diw=%d dih=%d fn=%s\n",movek,khrstart,khrend,diw,dih,fn_bam); fflush(stdout);  
    page_nts = khrend - khrstart;
    if (movek == 1)   // page right special code 
    {
         khrstart = khrend;
         khrend += page_nts;
    }
    else if (movek == -1)  // page left special code 
    {
         khrend  = khrstart;
         khrstart -= page_nts;
    }
    else
    {
        khrstart += movek;
        khrend   += movek;
    }
    sanity_khr();
    sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);

// set pos field - wtheck? done in get_params_and_draw
    NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
    [texts[0].object setStringValue:defval];

// xxx
printf("onmove before get_params,, pos = [%s]\n",pos); fflush(stdout);  
    get_params_and_draw(1);
    return 0;
}


void setup_xorea(int x, int w)
{
    struct image_type im2;
    if (img) 
    {
        if (img2) free(img2); 
        img2 = malloc(dih*diw*3);
        if (!img2) return;
        memcpy(img2,img,dih*diw*3);
        im2.data = img2;
        im2.width = diw;
        im2.height = dih;
        void ImageXRectangle(struct image_type *im,int x, int y, int x1, int y1 );
        ImageXRectangle(&im2,x, 0, x+w, dih );
        xoron = 1;
    }
}

void drawclip(void)
{
    char *sptr;
    int status;
    int clipx;
    int clipw;


    if (end_select_x > start_select_x) clipx = start_select_x;
    else clipx = end_select_x;
    if ( (start_select_x - end_select_x) < 0 ) clipw = end_select_x - start_select_x;
    else                                       clipw = start_select_x - end_select_x;
    setup_xorea(clipx,clipw);

    NSString *userinput = [texts[0].object stringValue];
    sptr = (char *)[userinput UTF8String];
    strcpy(pos,sptr);
    NSString *userinput2 = [texts[1].object stringValue];
    sptr = (char *)[userinput2 UTF8String];
    diw = atoi(sptr);
    NSString *userinput3 = [texts[2].object stringValue];
    sptr = (char *)[userinput3 UTF8String];
    dih = atoi(sptr);
    sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
    status = parse_position(pos,khr,&khrstart,&khrend); 

// put pos back to screen
    NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
    [texts[0].object setStringValue:defval];

    draw_mac(img2,diw,dih,2);
}
    


struct my_button_type
{
    int id; int x; int y; int w; int h;
    char *label;
};
struct my_button_type buttons[] =
{ // use usual x,y at northwest, fix in code
 { 1  ,  250 , 70,  60, 30, "Submit" } ,
 { 2  ,   10 ,110,  60, 30, "base" } ,
 { 3  ,   70 ,110,  60, 30, "<Page" } ,
 { 4  ,  130 ,110,  60, 30, "Page>" } ,
 { 5  ,  190 ,110,  60, 30, "lefthalf" } ,
 { 6  ,  250 ,110,  60, 30, "righthalf" } ,
 { 7  ,  310 ,110,  60, 30, "ZoomIN" } ,
 { 8  ,  370 ,110,  60, 30, "ZoomOut" } ,
 { 9  ,   10 ,140,  60, 30, "<10 " } ,
 { 10 ,   70 ,140,  60, 30, "10>" } ,
 { 11 ,  130 ,140,  60, 30, "<100" } ,
 { 12 ,  190 ,140,  60, 30, "100>" } ,
 { 13 ,  250 ,140,  60, 30, "<1000" } ,
 { 14 ,  310 ,140,  60, 30, "1000>" } ,
 { 15 ,  370 ,140,  60, 30, "<10000" } ,
 { 16 ,  430 ,140,  60, 30, ">10000" } ,
 { 17 ,  490 ,140,  60, 30, "<100000" } ,
 { 18 ,  550 ,140,  60, 30, ">100000" } ,
 { -1 ,   -1 , -1,  -1, -1, "ERRoverflow" } 
};

/*
 Position: Width: Height:    Splice Only On:   BaseColor:   Uniq:   Quality:
  7571199-7572901 , mismatch A  C  G  T  Ins  Del  mainmenu ins=0 del=0 tot=71 blds=hg19
Using Dataset:0 Position : chr17 7571199 7572901
Get DNA | Range is chr17 7571199 7572901 (1702 base pairs)
CGWB Link | UCSC Link | bambino | sam | sam1 | sam2 | fasta | FQ | align | blat | P | N | Pl | Nl | For ALVIEW Main Page, click here. 
NCBI: TP53 | all |
TRAWLER: TP53 | 
*/

// --- BUTTON 
@interface MyButtonClass:NSButton
{
@public
int myid;
}
@end

@implementation MyButtonClass
int myid;
// override mouseDown
- (void)mouseDown:(NSEvent *)theEvent;
{
    int size;
    int mid;
    int status; 

printf("start mousedown pos=[%s] dih = %d, diw = %d khr=[%s] s=%d e=%d\n",pos,dih,diw,khr,khrstart,khrend);

    if (self->myid == 1) // submit button 
    {
// xxx
        get_params_and_draw(1);
printf("in mousedown after get_params_and_drawpos=[%s] dih = %d, diw = %d khr=[%s] s=%d e=%d\n",pos,dih,diw,khr,khrstart,khrend);
    }
    else if (self->myid == 2) // { 2  ,   10 ,110,  60, 30, "base" } ,
    {
        mid = (khrend+khrstart)/2;
        khrstart = mid - (diw/2);
        khrend = khrstart + diw;
        sanity_khr();
        sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
        get_params_and_draw(1);
    }
    else if (self->myid == 3) // { 3  ,   70 ,110,  60, 30, "<Page" } ,
    {
        status = on_button_move(-1);
    }
    else if (self->myid == 4) // { 4  ,  130 ,110,  60, 30, "Page>" } ,
    {
        status = on_button_move(1);
    }
    else if (self->myid == 5) // { 5  ,  190 ,110,  60, 30, "lefthalf" } ,
    {
        mid = (khrend+khrstart)/2;
        khrend = mid;
        sanity_khr();
        sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
// put pos back 
        NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
        [texts[0].object setStringValue:defval];
        get_params_and_draw(1);
    }
    else if (self->myid == 6) // { 6  ,  250 ,110,  60, 30, "righthalf" } ,
    {
        mid = (khrend+khrstart)/2;
        khrstart = mid;
        sanity_khr();
        sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
// put pos back 
        NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
        [texts[0].object setStringValue:defval];
        get_params_and_draw(1);
    }
    else if (self->myid == 7) // { 7  ,  310 ,110,  60, 30, "ZoomIN" } ,
    {
        status = on_button_move(1);
        int third;

        size = khrend-khrstart;
        third=size/3;
        khrend = khrend - third;
        khrstart = khrstart + third;
        sanity_khr();
        sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
// put pos back 
        NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
        [texts[0].object setStringValue:defval];
        get_params_and_draw(1);
    }
    else if (self->myid == 8) // { 8  ,  370 ,110,  60, 30, "ZoomOut" } ,
    {
        int half;

        size = khrend-khrstart;
printf("Zoomout %s %d %d, size=%d\n",khr,khrstart,khrend,size); 
        half= size/2;
        khrend = khrend + half;
        khrstart = khrstart - half;
        sanity_khr();
        sprintf(pos,"%s:%d-%d",khr,khrstart,khrend);
// put pos back 
        NSString* defval = [NSString stringWithFormat:@"%s"  ,pos];
        [texts[0].object setStringValue:defval];

printf("Zoomout end %s %d %d pos=[%s]\n",khr,khrstart,khrend,pos); 
        get_params_and_draw(1);
    }
    else if (self->myid == 9) // { 9  ,   10 ,140,  60, 30, "<10 " } ,
    {
        status = on_button_move(-10);
    }
    else if (self->myid == 10) // { 10 ,   70 ,140,  60, 30, "10>" } ,
    {
        status = on_button_move(10);
    }
    else if (self->myid == 11) // { 11 ,  130 ,140,  60, 30, "<100" } ,
    {
        status = on_button_move(-100);
    }
    else if (self->myid == 12) // { 12 ,  190 ,140,  60, 30, "100>" } ,
    {
        status = on_button_move(100);
    }
    else if (self->myid == 13) // { 13 ,  250 ,140,  60, 30, "<1000" } ,
    {
        status = on_button_move(-1000);
    }
    else if (self->myid == 14) // { 14 ,  310 ,140,  60, 30, "1000>" } ,
    {
        status = on_button_move(1000);
    }
    else if (self->myid == 15) // { 15 ,  370 ,140,  60, 30, "<10000" } ,
    {
        status = on_button_move(-10000);
    }
    else if (self->myid == 16) // { 16 ,  430 ,140,  60, 30, ">10000" } ,
    {
        status = on_button_move(10000);
    }
    else if (self->myid == 17) // { 17 ,  490 ,140,  60, 30, "<100000" } ,
    {
        status = on_button_move(-100000);
    }
    else if (self->myid == 18) // { 18 ,  550 ,140,  60, 30, ">100000" } ,
    {
        status = on_button_move(100000);
    }
printf("end mousedown pos=[%s] dih = %d, diw = %d khr=[%s] s=%d e=%d\n",pos,dih,diw,khr,khrstart,khrend);
    [super mouseDown:theEvent];
}
@end


void open_file_dialog() // swiped from Rick's world
{
    int i; // Loop counter.
printf("in open_file_dialog\n"); 
    NSOpenPanel* openDlg = [NSOpenPanel openPanel]; // Create the File Open Dialog class.
    [openDlg setCanChooseFiles:YES]; // Enable the selection of files in the dialog.
    [openDlg setCanChooseDirectories:YES]; // Enable the selection of directories in the dialog.

// Display the dialog.  If the OK button was pressed, process the files.
    if ( [openDlg runModalForDirectory:nil file:nil] == NSOKButton )
    {
        // Get an array containing the full filenames of all files and directories selected.
        NSArray* files = [openDlg filenames];
        // Loop through all the files and process them.
        for( i = 0; i < [files count]; i++ )
        {
            NSString* fileName = [files objectAtIndex:i];
char *s = (char *)[fileName UTF8String];
printf("File name is %s\n",s);
if (i == 0) 
{
strcpy(fn_bam,s); 
printf("openfiledialog before get_params_and_draw [%s] ",fn_bam);
get_params_and_draw(1);
}
            // Do something with the filename.
        }
    }
}

MyButtonClass* create_mybutton(char *label, int id,float x, float y, float w, float h) 
{
    MyButtonClass *butt = [[MyButtonClass alloc] initWithFrame:CGRectMake(x,y,w,h)];
    butt->myid = id; 
    NSString* stringlabel = [NSString stringWithFormat:@"%s" , label];
    [butt setTitle:stringlabel];
    [butt setTarget:butt];
//    [butt setAction:@selector(someAction:)]; // left button click
    return butt;
}

MyTextClass* create_mytextview(int idarg,float x, float y, float w, float h, char *default_value) 
{
    MyTextClass  *tf = [[MyTextClass  alloc] initWithFrame:CGRectMake(x,y,w,h)];
    tf->myid = idarg; 
    // [tf.cell.setUsesSingleLineMode:true]; [tf.setWraps = false];
    NSString* defval = [NSString stringWithFormat:@"%s"  ,default_value];
    [tf setStringValue:defval];
    [tf setEditable: true];
    return tf;
}


NSTextField* create_label(char *label,float x, float y, float w, float h) 
{
    NSString* stringlabel = [NSString stringWithFormat:@"%s" , label];
    NSTextField  *lf = [[NSTextField  alloc] initWithFrame:CGRectMake(x,y,w,h)];
// printf("create_label %d %d %d %d\n",(int)x,(int)y,(int)w,(int)h); 
    [lf setStringValue: stringlabel];
    [lf setBezeled: false];
    [lf setDrawsBackground:false];
    [lf setSelectable:false];
    return lf;
}

@interface MyDelegate: NSObject
-(void)filestuff:(id)sender;
@end

@implementation MyDelegate: NSObject
-(void)filestuff:(id)sender;
{
    open_file_dialog();
}
@end

#include <sys/stat.h>

void check_genomedatadir(void)
{
    int bad = 0;
    struct stat sb;

    if (GENOMEDATADIR[0] == (char)0) 
        bad = 1;
    else
    {
        if (stat(GENOMEDATADIR, &sb) == 0 && S_ISDIR(sb.st_mode))
            bad = 0;
        else
            bad = 2;
    }
    if (bad)
    {
NSAlert *alert = [[[NSAlert alloc] init] autorelease];
[alert setMessageText:@"ERROR: no genome data directory.\nYou must set up the GENOMEDATADIR in text file alview.conf."];
[alert runModal];
    }
}


int main ()
{
    int i;
    [NSAutoreleasePool new];      // supports Cocoa's reference-counted memory management system.
    [NSApplication sharedApplication]; // manages an app's main event loop and resources use by all objects
          // You can use the global variabele NSApp to retieve the NSApplication instance

    [NSApp setActivationPolicy:NSApplicationActivationPolicyRegular]; 
    alview_load_config();
    check_genomedatadir();

    id menubar = [[NSMenu new] autorelease];
    id appMenuItem = [[NSMenuItem new] autorelease];
    [menubar addItem:appMenuItem];
    [NSApp setMainMenu:menubar];
    id appMenu = [[NSMenu new] autorelease];
    id appName = @"Alview Mac"; //id appName = [[NSProcessInfo processInfo] processName];

// file menu item
    id fileTitle = [@"File " stringByAppendingString:appName];
    id fileMenuItem = [[[NSMenuItem alloc] 
             initWithTitle:fileTitle action:@selector(filestuff:) keyEquivalent:@"f"] autorelease];
    id MyObject = [MyDelegate alloc];
    [fileMenuItem setTarget: MyObject]; 
    [appMenu addItem:fileMenuItem];
// quit menu item 
    id quitTitle = [@"Quit " stringByAppendingString:appName];
    id quitMenuItem = [[[NSMenuItem alloc] 
                 initWithTitle:quitTitle action:@selector(terminate:) keyEquivalent:@"q"] autorelease];
    [appMenu addItem:quitMenuItem];

    [appMenuItem setSubmenu:appMenu];

    id window = [[[NSWindow alloc] initWithContentRect:NSMakeRect(0, 0, 1000, 600)
        styleMask:NSTitledWindowMask | NSClosableWindowMask | NSMiniaturizableWindowMask | NSResizableWindowMask 
        backing:NSBackingStoreBuffered defer:NO] autorelease];
    [window cascadeTopLeftFromPoint:NSMakePoint(20,20)]; // rpf - careful dont want cascade  -right?
    [window setTitle:@"Alview for Mac"];
    [window makeKeyAndOrderFront:nil]; // makes it the "key" window

    NSView* superview=[[NSView alloc] initWithFrame:NSMakeRect(0,0,1000,800)];
    NSScrollView *scrollView = [[NSScrollView alloc] initWithFrame:[[window contentView] frame]];
    [scrollView setHasHorizontalScroller:YES];
    [scrollView setHasVerticalScroller:YES];
    [scrollView setBorderType:NSNoBorder];
    [scrollView setAutoresizingMask:NSViewWidthSizable|NSViewHeightSizable];
    [scrollView setDocumentView:superview];
    superframe = [ superview frame];
    NSRect subframe = [ superview frame]; // just to copy it !!!t
    image1  = [NSImage imageNamed:@"logo.png"];
    NSRect imageRect = NSMakeRect(0.0,500.0,[image1 size].width, [image1 size].height);
    theImageView = [[MyImageViewClass alloc] initWithFrame:imageRect];
[theImageView setTarget:theImageView];
    [theImageView setEnabled:true];
    [theImageView setBounds:imageRect];
    [theImageView setImage:image1];
    subframe.origin.x = 4;
    subframe.origin.y = superframe.size.height - 320;
    subframe.size.height = [image1 size].height;
    subframe.size.width =  [image1 size].width;
    [theImageView  setFrame:subframe];
    [superview addSubview:theImageView];

    NSRect subframe2 = [ superview frame];
    for (i=0 ; buttons[i].id != -1 ; i++)
    {
        MyButtonClass *butt2 = create_mybutton(
buttons[i].label, buttons[i].id,buttons[i].x,superframe.size.height - buttons[i].y,buttons[i].w,buttons[i].h);
        [superview addSubview:butt2];
    }

// label top
    NSTextField *tfhead = create_label(
               "Alview - Use Menu to Select BAM File to View - enter genome position (or gene name)",
               10,(int)(superframe.size.height-45),800,30);
    [superview addSubview:tfhead];

    for (i=0;texts[i].id != -1 ; i++)
    { // xxx 
        NSTextField *tf = create_label(texts[i].label,
                     texts[i].x,superframe.size.height - texts[i].y,texts[i].w,texts[i].h);
        [superview addSubview:tf];
        MyTextClass *txt1 = create_mytextview(
                     texts[i].id,texts[i].x2,
                     superframe.size.height - texts[i].y2,
                     texts[i].w2,texts[i].h2,
                     texts[i].default_value);
        texts[i].object = txt1;
        [superview addSubview:txt1];
    }
    [theImageView release];
    [window setContentView:scrollView];
    [scrollView release];
    [window makeKeyAndOrderFront:nil];

// *** force scroll to top left (this is because of the flipped coordinate system and nsscroll brain freeze)
    NSPoint pt = NSMakePoint(0.0, [[scrollView documentView] bounds].size.height);
    [[scrollView documentView] scrollPoint:pt];

    [NSApp activateIgnoringOtherApps:YES];
    [NSApp run];
    return 0;
}

// -- explain why this below is here ... ?
#include "alvmisc.c"
#include "alviewcore.cpp"

