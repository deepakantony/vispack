$ ! Command file to build the CFITSIO library on a VMS systems (VAX or Alpha)
$ set verify
$ cc buffers.c
$ cc cfileio.c
$ cc checksum.c
$ cc fitscore.c
$ cc getcol.c
$ cc getcolb.c
$ cc getcoli.c
$ cc getcolj.c
$ cc getcolui.c
$ cc getcoluj.c
$ cc getcolk.c
$ cc getcole.c
$ cc getcold.c
$ cc getcoll.c
$ cc getcols.c
$ cc editcol.c
$ cc edithdu.c
$ cc getkey.c
$ cc modkey.c
$ cc putcol.c
$ cc putcolb.c
$ cc putcoli.c
$ cc putcolj.c
$ cc putcolui.c
$ cc putcoluj.c
$ cc putcolk.c
$ cc putcole.c
$ cc putcold.c
$ cc putcols.c
$ cc putcoll.c
$ cc putcolu.c
$ cc putkey.c
$ cc scalnull.c
$ cc swapproc.c
$ cc wcsutil.c
$ lib/create cfitsio buffers,cfileio,checksum,fitscore,getcol,getcolb
$ lib/insert cfitsio getcoli,getcolj,getcolk,getcole,getcold,getcoll
$ lib/insert cfitsio getcols,editcol,edithdu,getkey,modkey,putcol
$ lib/insert cfitsio putcolb,putcoli,putcolj,putcolk,putcole,putcold
$ lib/insert cfitsio getcolui,getcoluj,putcolui,putcoluj
$ lib/insert cfitsio putcols,putcoll,putcolu,putkey,scalnull,swapproc,wcsutil
$ ! 
$ if (F$GETSYI("ARCH_NAME") .eqs. "VAX") then goto VAX
$   set noverify
$   exit
$ !
$ VAX:
$ ! add macro routines not needed on Alpha and only used on VAX
$   macro vmsieeer.mar
$   macro vmsieeed.mar
$   lib/insert cfitsio vmsieeer,vmsieeed
$   set noverify
$   exit
