
March 10, 2015     
ALVIEW NOW USES NATIVE GUIS ... GTK for Linux, Win32 for Windows and Cocoa for Mac.  It's fast, trim and zero licensing hassles.   Alview alsoruns as a web server or as a command line tool in batch. 

<p>The WEB DEMO is here : <a href="https://cgwb.nci.nih.gov/cgi-bin/alview">https://cgwb.nci.nih.gov/cgi-bin/alview</a></p>

Check BINARIES directory (above) for useable binaries.  Linux may require recompile from source; not a big deal though.

PCs and Unix compters are typically locked down these days in company, college and research environments.  This Alview distribution allows you to bypass installation restrictions on your work computer.  The downside is you have to a ~little~ bit of hassle work.  Not a big deal for 
wizards; but maybe a some careful work by novices.  You do have to download some big reference genome information.

Alview Runtime binaries rely only on standard runtime libraries.   No system updates necessary.
Just download the binary, some support files, edit alview.conf,  check permissions .... and GO VIEW YOUR BAMS.
Is that SNP real?   Check and make sure before you schedule a costly verification.  Make sure your snp call is not really noise.

If you wish to compile, please follow instructions in the ALVIEW_USERS_MANUAL.txt file. 

Aside from samtools, compilation and linking uses fairly standard libraries (zlib and run time libraries) .
Short, not complicated "go" files are provided instead of makefiles.  It's that simple: no cmake, 
no vjproj, no qmake, no Interface Builder, no nmake.exe etc. (You can use them if you want, of course).
Everything is pure vanilla low level C , except for the front ends : C++ for windows and Objective C for Mac Cocoa.

Want to port to a new OS?  Just re-implement the GUI front end.  Want to enhance it?  Go ahead.  You can edit and modify as you wish; it's public domain software.

_____

The Details.

ALVIEW is a PUBLIC DOMAIN/OPEN SOURCE bam file viewer.  It supports Windows, Linux and Mac.  It can run as 1) a native GUI program, 2) as a command line tool or 3) as a webserver (Apache/CGI).  Source is free and open.  The code is public domain (mostly).  Executables for Windows(64bit), MAC(64bit) and Linux(64bit) are provided.</p>

<p>The WEB DEMO is here : <a href="https://cgwb.nci.nih.gov/cgi-bin/alview">https://cgwb.nci.nih.gov/cgi-bin/alview</a></p>

For now, you must also download the large genome data support files targzipped in the file named 
alview_hg18hg19_genomedata.tar.bz2 from my dropbox account at <a href="https://www.dropbox.com/sh/ymszmksj6v83rmt/OCLsYTjgSu">https://www.dropbox.com/sh/ymszmksj6v83rmt/OCLsYTjgSu</a> .  This is large (&gt;1GB) file and takes a while to bzip2 decompress and "untar".  

Grab your executable BINARIES directory hereabouts on Github.  

<p>The TWO command lines to prepare this file (on LINUX, Mac command line or cygwin [ or cmd.exe with unix utilities ] ) ... 

<p>bzip2 -d alview_hg18hg19_genomedata.tar.bz2</p>

<p>tar xvf alview_hg18hg19_genomedata.tar </p>

Put this stuff in a convenient location.  In a directory called GENOMEDATA in the place you put the executable is a good place. 
Set up the "alview.conf" file.
Check permissions.
Run alview.

<p>Read the instructions for installing for your particular operating system (Linux, Mac or Windows).
The instructions are available in ALVIEW_USERS_MANUAL.txt file in the list of "root" level alview files.</p>

You are welcome to customize and re-brand this as you see fit. NCI and NIH retain the many trademarks associated with the National Institutes of Health and the National Cancer Institute.</p>

