rem:  this batch file builds the cfitsio library 
rem:  using the Borland C++ v4.5 compiler
rem:
bcc32 -c buffers.c
bcc32 -c cfileio.c
bcc32 -c checksum.c
bcc32 -c fitscore.c
bcc32 -c getcol.c
bcc32 -c getcolb.c
bcc32 -c getcoli.c
bcc32 -c getcolj.c
bcc32 -c getcolui.c
bcc32 -c getcoljj.c
bcc32 -c getcolk.c
bcc32 -c getcole.c
bcc32 -c getcold.c
bcc32 -c getcoll.c
bcc32 -c getcols.c
bcc32 -c editcol.c
bcc32 -c edithdu.c
bcc32 -c getkey.c
bcc32 -c modkey.c
bcc32 -c putcol.c
bcc32 -c putcolb.c
bcc32 -c putcoli.c
bcc32 -c putcolj.c
bcc32 -c putcolui.c
bcc32 -c putcoluj.c
bcc32 -c putcolk.c
bcc32 -c putcole.c
bcc32 -c putcold.c
bcc32 -c putcols.c
bcc32 -c putcoll.c
bcc32 -c putcolu.c
bcc32 -c putkey.c
bcc32 -c scalnull.c
bcc32 -c swapproc.c
bcc32 -c wcsutil.c
del cfitsio.lib
tlib cfitsio +buffers +cfileio +checksum +fitscore
tlib cfitsio +getcol  +getcolb +getcoli +getcolj +getcolk +getcole +getcold
tlib cfitsio +getcoll +getcols +editcol +edithdu +getkey  +modkey  +putkey 
tlib cfitsio +putcol  +putcolb +putcoli +putcolj +putcolk +putcole +putcold
tlib cfitsio +putcoll +putcols +putcolu +scalnull +swapproc +wcsutil
tlib cfitsio +getcolui +getcoluj +putcolui +putcoluj
bcc32 -f testprog.c cfitsio.lib
bcc32 -f cookbook.c cfitsio.lib

