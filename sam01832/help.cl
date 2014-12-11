                         C/C++ COMPILER OPTIONS


                              -OPTIMIZATION-

/O1 minimize space                      /O2 maximize speed
/Ob<n> inline expansion (default n=0)   /Od disable optimizations (default)
/Og enable global optimization          /Oi[-] enable intrinsic functions
/Os favor code space                    /Ot favor code speed
/Ox maximum optimizations               
/favor:<blend|AMD64|INTEL64|ATOM> select processor to optimize for, one of:
    blend - a combination of optimizations for several different x64 processors
    AMD64 - 64-bit AMD processors                                 
    INTEL64 - Intel(R)64 architecture processors                  
    ATOM - Intel(R) Atom(TM) processors                           

                             -CODE GENERATION-

/GF enable read-only string pooling     /Gm[-] enable minimal rebuild
/Gy[-] separate functions for linker    /GS[-] enable security checks
/GR[-] enable C++ RTTI                  /GX[-] enable C++ EH (same as /EHsc)
/EHs enable C++ EH (no SEH exceptions)  /EHa enable C++ EH (w/ SEH exceptions)
/EHc extern "C" defaults to nothrow     
/fp:<except[-]|fast|precise|strict> choose floating-point model:
    except[-] - consider floating-point exceptions when generating code
    fast - "fast" floating-point model; results are less predictable
    precise - "precise" floating-point model; results are predictable
    strict - "strict" floating-point model (implies /fp:except)
/Qfast_transcendentals generate inline FP intrinsics even with /fp:except
/Qpar[-] enable parallel code generation
/Qpar-report:1 auto-parallelizer diagnostic; indicate parallelized loops
/Qpar-report:2 auto-parallelizer diagnostic; indicate loops not parallelized
/Qvec-report:1 auto-vectorizer diagnostic; indicate vectorized loops
/Qvec-report:2 auto-vectorizer diagnostic; indicate loops not vectorized
/GL[-] enable link-time code generation 
/volatile:<iso|ms> choose volatile model:
    iso - Acquire/release semantics not guaranteed on volatile accesses
    ms  - Acquire/release semantics guaranteed on volatile accesses
/GA optimize for Windows Application    /Ge force stack checking for all funcs
/Gs[num] control stack checking calls   /Gh enable _penter function call
/GH enable _pexit function call         /GT generate fiber-safe TLS accesses
/RTC1 Enable fast checks (/RTCsu)       /RTCc Convert to smaller type checks
/RTCs Stack Frame runtime checking      /RTCu Uninitialized local usage checks
/clr[:option] compile for common language runtime, where option is:
    pure - produce IL-only output file (no native executable code)
    safe - produce IL-only verifiable output file
    oldSyntax - accept the Managed Extensions syntax from Visual C++ 2002/2003
    initialAppDomain - enable initial AppDomain behavior of Visual C++ 2002
    noAssembly - do not produce an assembly
    nostdlib - ignore the default \clr directory
/homeparams Force parameters passed in registers to be written to the stack
/GZ Enable stack checks (/RTCs)         
/arch:AVX enable use of Intel(R) Advanced Vector Extensions instructions

                              -OUTPUT FILES-

/Fa[file] name assembly listing file    /FA[scu] configure assembly listing
/Fd[file] name .PDB file                /Fe<file> name executable file
/Fm[file] name map file                 /Fo<file> name object file
/Fp<file> name precompiled header file  /Fr[file] name source browser file
/FR[file] name extended .SBR file       /Fi[file] name preprocessed file
/doc[file] process XML documentation comments and optionally name the .xdc file

                              -PREPROCESSOR-

/AI<dir> add to assembly search path    /FU<file> forced using assembly/module 
/C don't strip comments                 /D<name>{=|#}<text> define macro
/E preprocess to stdout                 /EP preprocess to stdout, no #line
/P preprocess to file                   /Fx merge injected code to file
/FI<file> name forced include file      /U<name> remove predefined macro
/u remove all predefined macros         /I<dir> add to include search path
/X ignore "standard places"             

                                -LANGUAGE-

/Zi enable debugging information        /Z7 enable old-style debug info
/Zp[n] pack structs on n-byte boundary  /Za disable extensions
/Ze enable extensions (default)         /Zl omit default library name in .OBJ
/Zg generate function prototypes        /Zs syntax check only
/vd{0|1|2} disable/enable vtordisp      /vm<x> type of pointers to members
/Zc:arg1[,arg2] C++ language conformance, where arguments can be:
    forScope[-] - enforce Standard C++ for scoping rules
    wchar_t[-] - wchar_t is the native type, not a typedef
    auto[-] - enforce the new Standard C++ meaning for auto
    trigraphs[-] - enable trigraphs (off by default)
/ZW enable WinRT language extensions    
/openmp enable OpenMP 2.0 language extensions

                              -MISCELLANEOUS-

@<file> options response file           /?, /help print this help message
/bigobj generate extended object format /c compile only, no link
/errorReport:option Report internal compiler errors to Microsoft
    none - do not send report                
    prompt - prompt to immediately send report
    queue - at next admin logon, prompt to send report (default)
    send - send report automatically         
/FC use full pathnames in diagnostics   /H<num> max external name length
/J default char type is unsigned        
/MP[n] use up to 'n' processes for compilation
/nologo suppress copyright message      
/sdl enable additional security features and warnings
/showIncludes show include file names   /Tc<source file> compile file as .c
/Tp<source file> compile file as .cpp   /TC compile all files as .c
/TP compile all files as .cpp           /V<string> set version string
/w disable all warnings                 /wd<n> disable warning n
/we<n> treat warning n as an error      /wo<n> issue warning n once
/w<l><n> set warning level 1-4 for n    /W<n> set warning level (default n=1)
/Wall enable all warnings               /WL enable one line diagnostics
/WX treat warnings as errors            /Yc[file] create .PCH file
/Yd put debug info in every .OBJ        /Yl[sym] inject .PCH ref for debug lib
/Yu[file] use .PCH file                 /Y- disable all PCH options
/Zm<n> max memory alloc (% of default)  /Wp64 enable 64 bit porting warnings

                                -LINKING-

/LD Create .DLL                         /LDd Create .DLL debug library
/LN Create a .netmodule                 /F<num> set stack size
/link [linker options and libraries]    /MD link with MSVCRT.LIB
/MT link with LIBCMT.LIB                /MDd link with MSVCRTD.LIB debug lib
/MTd link with LIBCMTD.LIB debug lib    


