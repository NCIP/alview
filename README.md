
March 17, 2015  

<p>The WEB DEMO is here : <a href="https://cgwb.nci.nih.gov/cgi-bin/alview">https://cgwb.nci.nih.gov/cgi-bin/alview</a></p>
<p>Installers for Mac/Windows/Linux are here : <a href="http://45.56.125.191/">http://45.56.125.191/</a></p>

If your system is locked down and you cannot install the package for your Operating System, then do this: .... create a directory to store alview.  grab binary from BINARIES directons ( from this NCIP Alview Front Page)  <a href="https://github.com/NCIP/alview/tree/master/BINARIES">  https://github.com/NCIP/alview/tree/master/BINARIES </a> , download and untar the Geneome suport data (this is very large) <a href="http://45.56.125.191/AlviewGenomeData.tar> http://45.56.125.191/AlviewGenomeData.tar</a>, edit the sample alview.conf file (sample at <a href="http://45.56.125.191/alview.conf"> http://45.56.125.191/alview.conf </a> to point to GENOMEDATADIR directory.
Then check permissions and you can run Alview.   Linux users will occasionally need to re-compile.

Alview is a very fast bam file viewer.  GUI versions now support GTK for Linux, Win32 for Windows and Cocoa for Mac.  It's fast, trim and there are zero licensing hassles.   Alview also runs as a web server or as a command line tool in batch. 
 
Alview Runtime binaries only typically only rely on standard runtime libraries.  System updates are seldom necessary.

If you wish to compile, please follow instructions in the ALVIEW_USERS_MANUAL.txt file. 
Aside from samtools, compilation and linking uses fairly standard libraries (zlib and run time libraries) .
Short, not complicated "go" files are provided instead of makefiles.  It's that simple: no cmake, 
no vjproj, no qmake, no Interface Builder, no nmake.exe etc. (You can use them if you want, of course).
Everything is pure vanilla low level C , except for the front ends : C++ for windows and Objective C for Mac Cocoa.

Want to port to a new OS?  Just re-implement the GUI front end.  Want to enhance it?  Go ahead.  You can edit and modify as you wish; it's public domain software.

<p>Detailed instructsion are available in the ALVIEW_USERS_MANUAL.txt file in the list of "root" level alview files.</p>

You are welcome to customize and re-brand this as you see fit. NCI and NIH retain the many trademarks associated with the National Institutes of Health and the National Cancer Institute.</p>

