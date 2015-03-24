
March 24, 2015  

<p>The Alview WEB DEMO is here : <a href="https://cgwb.nci.nih.gov/cgi-bin/alview">https://cgwb.nci.nih.gov/cgi-bin/alview</a></p>
<p>Installers for Mac/Windows/Linux are here : <a href="http://45.56.125.191/">http://45.56.125.191/</a></p>
<p>User's Manual is here : <a href="https://raw.githubusercontent.com/NCIP/alview/master/ALVIEW_USERS_MANUAL.txt">
 https://raw.githubusercontent.com/NCIP/alview/master/ALVIEW_USERS_MANUAL.txt . </p>
 
Current hg18 and hg19 are supported.  More genomes are coming.

If your system is locked down and you cannot install the package for your Operating System, then do this: .... create a directory to store alview program and data.  Then grab a binary from BINARIES directory ( from this NCIP Alview Github Front Page)  <a href="https://github.com/NCIP/alview/tree/master/BINARIES">  https://github.com/NCIP/alview/tree/master/BINARIES </a> , download and untar the Genome support data ("AlviewGenomeData.tar") in your new directory (this is very large) <a href="http://45.56.125.191/AlviewGenomeData.tar> http://45.56.125.191/AlviewGenomeData.tar</a>, Then download and edit the sample alview.conf file (sample at <a href="http://45.56.125.191/alview.conf"> http://45.56.125.191/alview.conf </a> to point to GENOMEDATADIR directory.
Then check permissions and you can run Alview.  Set all files to readable and set the binary to executable.   Many recent Linuxes are compatible with older binaries; but some Linux users might need to re-compile Alview.  Alview Linux relies on GTK2 which is often installed on systems.

Alview is a very fast bam file viewer.  GUI versions now support GTK for Linux, Win32 for Windows and Cocoa for Mac.  It's designed and implemented for speed and there are zero licensing hassles, no restrictions and no registration.   Alview also runs as a web server or as a command line tool to generate short read alignment images. This can of course be scripted in a batch file.  Power users will find this very useful for reviewing thousands of snp calls.
 
Alview  binaries rely on reasonably standard runtime libraries.  System updates are seldom necessary.

If you wish to compile, please follow instructions in the ALVIEW_USERS_MANUAL.txt file. 
Aside from samtools, compilation and linking uses very common libraries (zlib and OS run time libraries and gtk2 if on Linux) .
Short, not complicated "go" files are provided instead of makefiles.  It's that simple: no cmake, 
no vjproj, no qmake, no Interface Builder, no nmake.exe etc. (You can use them if you want, of course).
Everything is pure vanilla low level C , except for the front ends : C++ for windows and Objective C for Mac Cocoa.

Want to port to a new OS?  Just re-implement the GUI front end.  Want to enhance it?  Go ahead.  You can edit and modify as you wish; it's public domain software.

<p>Detailed instructions are available in the ALVIEW_USERS_MANUAL.txt file in the list of "root" level alview files.</p>

You are welcome to customize and re-brand this as you see fit. NCI and NIH retain the many trademarks associated with the National Institutes of Health and the National Cancer Institute.</p>

