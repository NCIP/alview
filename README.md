
December 11, 2014 MAJOR UPDATE IN PROGRESS, WORKABLE BETAS READY BY MONDAY.

ALVIEW NOW USES NATIVE GUIS ... GTK for Linux, Win32 for Windows and Cocoa for Mac !!!
Usable Windows and Linux (sort of) binaries are up now.
Need to test some Mac stuff before uploading; though source is there.
 

 Distribution of binaries (EXEs) is much now much simpler. No install program.  These are binaries that only rely on standard libraries. 
Making the program is easy.  Compilation uses standard libraries and tools (link, gcc, clang, cl.exe, rc.exe; depending on your OS)
Compiling to a NATIVE binary is just a few keystrokes.  No complicated qmake, cmake, configure, ant, vcproj files, etc.
Everything is pure vanilla.

You can edit and modify as you wish; it's public domain software.
No complicated Cmakes, makefiles, vcprojs, qmakes, etc.  QT and DLL hell is going away.
Just download the binary, some support files, check permissions .... and GO VIEW YOUR BAMS.</p>

_____

 ALVIEW is a PUBLIC DOMAIN/OPEN SOURCE bam file viewer.  It supports Windows, Linux and Mac.  It can run as 1) a native GUI program, 2) as a command line tool or 3) as a webserver (Apache/CGI).  Source is free/open.  The code is public domain (mostly).  Executables for Windows(64bit), MAC(64bit) and Linux(64bit) are provided.</p>

<p>The WEB DEMO is here : <a href="https://cgwb.nci.nih.gov/cgi-bin/alview">https://cgwb.nci.nih.gov/cgi-bin/alview</a></p>

FOr now must also download the large genome data support files targzipped in the file named 
alview_hg18hg19_genomedata.tar.bz2 from my dropbox account at <a href="https://www.dropbox.com/sh/ymszmksj6v83rmt/OCLsYTjgSu">https://www.dropbox.com/sh/ymszmksj6v83rmt/OCLsYTjgSu</a> .  This is large (&gt;1GB) file and takes a while to 
decompress and "untar".  

Grab your exectuable BINARIES directory hereabouts on Github.  
Set up the "alview.conf" file.
 

<p>The TWO command lines to prepare this file (on LINUX, Mac command line or cygwin [ or cmd.exe with unix utilites ] ) ... 

<p>bzip2 -d alview_hg18hg19_genomedata.tar.bz2</p>

<p>tar xvf alview_hg18hg19_genomedata.tar </p>

<p>Read the instructions for installing for your particular operating system (Linux, Mac or Windows).
Instructions are vailable in ALVIEW_USERS_MANAUAL.txt file in the list of "root" level alview files.</p>

 You are welcome to customize and re-brand this as you see fit.    NCI and NIH retain the many trademarks associated with the  National Institutes of Health and the National Cancer Institute.</p>
