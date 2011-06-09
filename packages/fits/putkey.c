/*  This file, putkey.c, contains routines that write keywords to          */
/*  a FITS header.                                                         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.  Users shall not, without prior written   */
/*  permission of the U.S. Government,  establish a claim to statutory     */
/*  copyright.  The Government and others acting on its behalf, shall have */
/*  a royalty-free, non-exclusive, irrevocable,  worldwide license for     */
/*  Government purposes to publish, distribute, translate, copy, exhibit,  */
/*  and perform such material.                                             */

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "fitsio2.h"
/*--------------------------------------------------------------------------*/
int ffcrim(fitsfile *fptr,      /* I - FITS file pointer           */
           int bitpix,          /* I - bits per pixel              */
           int naxis,           /* I - number of axes in the array */
           long *naxes,         /* I - size of each axis           */
           int *status)         /* IO - error status               */
/*
  create an IMAGE extension following the current HDU. If the
  current HDU is empty (contains no header keywords), then simply
  write the required image (or primary array) keywords to the current
  HDU. 
*/
{
    if (*status > 0)
        return(*status);

    /* create new extension if current header is not empty */
    if (fptr->headend != fptr->headstart[fptr->curhdu] )
        ffcrhd(fptr, status);

    /* write the required header keywords */
    ffphpr(fptr, TRUE, bitpix, naxis, naxes, 0, 1, TRUE, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffcrtb(fitsfile *fptr,  /* I - FITS file pointer                        */
           int tbltype,     /* I - type of table to create                  */
           long naxis2,     /* I - number of rows in the table              */
           int tfields,     /* I - number of columns in the table           */
           char **ttype,    /* I - name of each column                      */
           char **tform,    /* I - value of TFORMn keyword for each column  */
           char **tunit,    /* I - value of TUNITn keyword for each column  */
           char *extnm,   /* I - value of EXTNAME keyword, if any         */
           int *status)     /* IO - error status                            */
/*
  Create a table extension in a FITS file. 
*/
{
    long naxis1, ncols, *tbcol;

    if (*status > 0)
        return(*status);

    /* create new extension if current header is not empty */
    if (fptr->headend != fptr->headstart[fptr->curhdu] )
        ffcrhd(fptr, status);

    if (tbltype == BINARY_TBL)
    {
      /* write the required header keywords. This will write PCOUNT = 0 */
      /* so variable length array columns are not supported             */
      ffphbn(fptr, naxis2, tfields, ttype, tform, tunit, extnm, 0, status);
    }
    else if (tbltype == ASCII_TBL)
    {
      /* allocate mem for tbcol; malloc can have problems allocating small */
      /* arrays, so allocate at least 20 bytes */

      ncols = maxvalue(5, tfields);
      tbcol = (long *) calloc(ncols, sizeof(long));

      if (tbcol)
      {
        /* calculate width of a row and starting position of each column. */
        /* Each column will be separated by 1 blank space */
        ffgabc(tfields, tform, 1, &naxis1, tbcol, status);

        /* write the required header keywords */
        ffphtb(fptr, naxis1, naxis2, tfields, ttype, tbcol, tform, tunit,
               extnm, status);

        free(tbcol);
      }
    }
    else
      *status = NOT_TABLE;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpky( fitsfile *fptr,     /* I - FITS file pointer        */
           int  datatype,      /* I - datatype of the value    */
           char *keyname,      /* I - name of keyword to write */
           void *value,        /* I - keyword value            */
           char *comm,         /* I - keyword comment          */
           int  *status)       /* IO - error status            */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes a keyword value with the datatype specified by the 2nd argument.
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (datatype == TSTRING)
    {
        ffpkys(fptr, keyname, (char *) value, comm, status);
    }
    else if (datatype == TBYTE)
    {
        ffpkyj(fptr, keyname, (long) *(unsigned char *) value, comm, status);
    }
    else if (datatype == TUSHORT)
    {
        ffpkyj(fptr, keyname, (long) *(unsigned short *) value, comm, status);
    }
    else if (datatype == TSHORT)
    {
        ffpkyj(fptr, keyname, (long) *(short *) value, comm, status);
    }
    else if (datatype == TINT)
    {
        ffpkyj(fptr, keyname, (long) *(int *) value, comm, status);
    }
    else if (datatype == TLOGICAL)
    {
        ffpkyl(fptr, keyname, *(int *) value, comm, status);
    }
    else if (datatype == TULONG)
    {
        ffpkyj(fptr, keyname, (long) *(unsigned long *) value, comm, status);
    }
    else if (datatype == TLONG)
    {
        ffpkyj(fptr, keyname, *(long *) value, comm, status);
    }
    else if (datatype == TFLOAT)
    {
        ffpkye(fptr, keyname, *(float *) value, 6, comm, status);
    }
    else if (datatype == TDOUBLE)
    {
        ffpkyd(fptr, keyname, *(double *) value, 14, comm, status);
    }
    else
        *status = BAD_DATATYPE;

    return(*status);
} 
/*-------------------------------------------------------------------------*/
int ffprec(fitsfile *fptr,     /* I - FITS file pointer        */
           const char *card,   /* I - string to be written     */
           int *status)        /* IO - error status            */
/*
  write a keyword record (80 bytes long) to the end of the header
*/
{
    char tcard[81];
    size_t len, ii;
    long nblocks;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if ( (fptr->datastart - fptr->headend) == 80) /* only room for END card */
    {
        nblocks = 1;
        if (ffiblk(fptr, nblocks, 0, status) > 0) /* insert 2880-byte block */
            return(*status);  
    }

    strncpy(tcard,card,80);
    tcard[80] = '\0';

    len = strlen(tcard);
    for (ii=len; ii < 80; ii++)    /* fill card with spaces if necessary */
        tcard[ii] = ' ';

    for (ii=0; ii < 8; ii++)       /* make sure keyword name is uppercase */
        tcard[ii] = toupper(tcard[ii]);

    fftkey(tcard, status);        /* test keyword name contains legal chars */

    fftrec(tcard, status);          /* test rest of keyword for legal chars */

    ffmbyt(fptr, fptr->headend, IGNORE_EOF, status); /* move to end header */

    ffpbyt(fptr, 80, tcard, status);   /* write the 80 byte card */

    if (*status <= 0)
       fptr->headend += 80;           /* update end-of-header position */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkys( fitsfile *fptr,     /* I - FITS file pointer        */
            char *keyname,      /* I - name of keyword to write */
            char *value,        /* I - keyword value            */
            char *comm,         /* I - keyword comment          */
            int  *status)       /* IO - error status            */
/*
  Write (put) the keyword, value and comment into the FITS header.
  The value string will be truncated at 68 characters which is the
  maximum length that will fit on a single FITS keyword.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffs2c(value, valstring, status);   /* put quotes around the string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkls( fitsfile *fptr,     /* I - FITS file pointer        */
            char *keyname,      /* I - name of keyword to write */
            char *value,        /* I - keyword value            */
            char *comm,         /* I - keyword comment          */
            int  *status)       /* IO - error status            */
/*
  Write (put) the keyword, value and comment into the FITS header.
  This routine is a modified version of ffpkys which supports the
  HEASARC long string convention and can write arbitrarily long string
  keyword values.  The value is continued over multiple keywords that
  have the name COMTINUE without an equal sign in column 9 of the card.
  This routine also supports simple string keywords which are less than
  69 characters in length.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];
    char tstring[FLEN_VALUE], *cptr;
    int next, remain, vlen, nquote, nchar, contin;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    remain = strlen(value);    /* number of characters to write out */
    next = 0;                  /* pointer to next character to write */
    
    /* count the number of single quote characters are in the string */
    nquote = 0;
    cptr = strchr(value, '\'');   /* search for quote character */

    while (cptr)  /* search for quote character */
    {
        nquote++;            /*  increment no. of quote characters  */
        cptr++;              /*  increment pointer to next character */
        cptr = strchr(cptr, '\'');  /* search for another quote char */
    }

    /* each quote character is expanded to 2 quotes, so leave enough space */
    nchar = 68 - nquote;    /*  max of 68 chars fit in a FITS string value */

    contin = 0;
    while (remain > 0)
    {
        strncpy(tstring, &value[next], nchar); /* copy string to temp buff */
        tstring[nchar] = '\0';
        ffs2c(tstring, valstring, status);  /* put quotes around the string */

        if (remain > nchar)   /* if string is continued, put & as last char */
        {
            vlen = strlen(valstring);
            nchar -= 1;        /* outputting one less character now */

            if (valstring[vlen-2] != '\'')
                valstring[vlen-2] = '&';  /*  over write last char with &  */
            else
            { /* last char was a pair of single quotes, so over write both */
                valstring[vlen-3] = '&';
                valstring[vlen-1] = '\0';
            }
        }

        ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/

        if (contin)           /* This is a CONTINUEd keyword */
            strncpy(card, "CONTINUE   ", 10);  /* overwrite the name and = */
 
        ffprec(fptr, card, status);  /* write the keyword */

        contin = 1;
        remain -= nchar;
        next  += nchar;
        nchar += 1;    /* this compensates for the earlier decrement */
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffplsw( fitsfile *fptr,     /* I - FITS file pointer  */
            int  *status)       /* IO - error status       */
/*
  Write the LONGSTRN keyword and a series of related COMMENT keywords
  which document that this FITS header may contain long string keyword
  values which are continued over multiple keywords using the HEASARC
  long string keyword convention.  If the LONGSTRN keyword already exists
  then this routine simple returns without doing anything.
*/
{
    char valstring[FLEN_VALUE], comm[FLEN_COMMENT];
    int tstatus;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    tstatus = 0;
    if (ffgkys(fptr, "LONGSTRN", valstring, comm, &tstatus) == 0)
        return(*status);     /* keyword already exists, so just return */

    ffpkys(fptr, "LONGSTRN", "OGIP 1.0", 
       "The HEASARC Long String Convention may be used.", status);

    ffpcom(fptr,
    "This FITS file may contain long string keyword values that are", status);

    ffpcom(fptr,
    "continued over multiple keywords.  The HEASARC convention uses the &",
    status);

    ffpcom(fptr,
    "character at the end of each substring which is then continued", status);

    ffpcom(fptr,
    "on the next keyword which has the name CONTINUE.", status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyl( fitsfile *fptr,     /* I - FITS file pointer        */
            char *keyname,      /* I - name of keyword to write */
            int  value,         /* I - keyword value            */
            char *comm,         /* I - keyword comment          */
            int  *status)       /* IO - error status            */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Values equal to 0 will result in a False FITS keyword; any other
  non-zero value will result in a True FITS keyword.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffl2c(value, valstring, status);   /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyj( fitsfile *fptr,     /* I - FITS file pointer        */
            char *keyname,      /* I - name of keyword to write */
            long value,         /* I - keyword value            */
            char *comm,         /* I - keyword comment          */
            int  *status)       /* IO - error status            */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes an integer keyword value.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffi2c(value, valstring, status);   /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyf( fitsfile *fptr,      /* I - FITS file pointer                   */
            char  *keyname,      /* I - name of keyword to write            */
            float value,         /* I - keyword value                       */
            int   decim,         /* I - number of decimal places to display */
            char  *comm,         /* I - keyword comment                     */
            int   *status)       /* IO - error status                       */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes a fixed float keyword value.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffr2f(value, decim, valstring, status);   /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkye( fitsfile *fptr,      /* I - FITS file pointer                   */
            char  *keyname,      /* I - name of keyword to write            */
            float value,         /* I - keyword value                       */
            int   decim,         /* I - number of decimal places to display */
            char  *comm,         /* I - keyword comment                     */
            int   *status)       /* IO - error status                       */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes an exponential float keyword value.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffr2e(value, decim, valstring, status);   /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyg( fitsfile *fptr,      /* I - FITS file pointer                   */
            char  *keyname,      /* I - name of keyword to write            */
            double value,        /* I - keyword value                       */
            int   decim,         /* I - number of decimal places to display */
            char  *comm,         /* I - keyword comment                     */
            int   *status)       /* IO - error status                       */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes a fixed double keyword value.*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffd2f(value, decim, valstring, status);  /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyd( fitsfile *fptr,      /* I - FITS file pointer                   */
            char  *keyname,      /* I - name of keyword to write            */
            double value,        /* I - keyword value                       */
            int   decim,         /* I - number of decimal places to display */
            char  *comm,         /* I - keyword comment                     */
            int   *status)       /* IO - error status                       */
/*
  Write (put) the keyword, value and comment into the FITS header.
  Writes an exponential double keyword value.*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffd2e(value, decim, valstring, status);  /* convert to formatted string */
    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkyt( fitsfile *fptr,      /* I - FITS file pointer        */
            char  *keyname,      /* I - name of keyword to write */
            long  intval,        /* I - integer part of value    */
            double fraction,     /* I - fractional part of value */
            char  *comm,         /* I - keyword comment          */
            int   *status)       /* IO - error status            */
/*
  Write (put) a 'triple' precision keyword where the integer and
  fractional parts of the value are passed in separate parameters to
  increase the total amount of numerical precision.
*/
{
    char valstring[FLEN_VALUE];
    char card[FLEN_CARD];
    char fstring[20], *cptr;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (fraction > 1. || fraction < 0.)
        return(*status = BAD_F2C);

    ffi2c(intval, valstring, status);  /* convert integer to string */
    ffd2f(fraction, 16, fstring, status);  /* convert to 16 decimal string */

    cptr = strchr(fstring, '.');    /* find the decimal point */
    strcat(valstring, cptr);    /* append the fraction to the integer */

    ffmkky(keyname, valstring, comm, card);  /* construct the keyword*/
    ffprec(fptr, card, status);  /* write the keyword*/

    return(*status);
}
/*-----------------------------------------------------------------*/
int ffpcom( fitsfile *fptr,      /* I - FITS file pointer   */
            const char  *comm,   /* I - comment string      */
            int   *status)       /* IO - error status       */
/*
  Write 1 or more COMMENT keywords.  If the comment string is too
  long to fit on a single keyword (70 chars) then it will automatically
  be continued on multiple CONTINUE keywords.
*/
{
    char card[FLEN_CARD];
    int len, ii;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    len = strlen(comm);
    ii = 0;

    for (; len > 0; len -= 70)
    {
        strcpy(card, "COMMENT   ");
        strncat(card, &comm[ii], 70);
        ffprec(fptr, card, status);
        ii += 70;
    }

    return(*status);
}
/*-----------------------------------------------------------------*/
int ffphis( fitsfile *fptr,      /* I - FITS file pointer  */
            const char *history, /* I - history string     */
            int   *status)       /* IO - error status      */
/*
  Write 1 or more HISTORY keywords.  If the history string is too
  long to fit on a single keyword (70 chars) then it will automatically
  be continued on multiple HISTORY keywords.
*/
{
    char card[FLEN_CARD];
    int len, ii;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    len = strlen(history);
    ii = 0;

    for (; len > 0; len -= 70)
    {
        strcpy(card, "HISTORY   ");
        strncat(card, &history[ii], 70);
        ffprec(fptr, card, status);
        ii += 70;
    }

    return(*status);
}
/*-----------------------------------------------------------------*/
int ffpdat( fitsfile *fptr,      /* I - FITS file pointer  */
            int   *status)       /* IO - error status      */
/*
  Write the DATE keyword into the FITS header.  If the keyword already
  exists then the date will simply be updated in the existing keyword.
*/
{
    char date[9], card[FLEN_CARD];
    time_t tp;
    struct tm *ptr;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    time(&tp);
    ptr = localtime(&tp);
    strftime(date, 9, "%d/%m/%y", ptr);

    strcpy(card, "DATE    = '");
    strcat(card, date);
    strcat(card, "'           / FITS file creation date (dd/mm/yy)");

    ffucrd(fptr, "DATE", card, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkns( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            char *value[],      /* I - array of pointers to keyword values  */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes string keywords.
  The value strings will be truncated at 68 characters, and the HEASARC
  long string keyword convention is not supported by this routine.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkys(fptr, keyname, value[ii], tcomment, status);
        else
            ffpkys(fptr, keyname, value[ii], comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpknl( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            int  *value,        /* I - array of keyword values              */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes logical keywords
  Values equal to zero will be written as a False FITS keyword value; any
  other non-zero value will result in a True FITS keyword.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);

        if (repeat)
            ffpkyl(fptr, keyname, value[ii], tcomment, status);
        else
            ffpkyl(fptr, keyname, value[ii], comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpknj( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            long *value,        /* I - array of keyword values              */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Write integer keywords
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkyj(fptr, keyname, value[ii], tcomment, status);
        else
            ffpkyj(fptr, keyname, value[ii], comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpknf( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            float *value,       /* I - array of keyword values              */
            int decim,          /* I - number of decimals to display        */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed float values.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkyf(fptr, keyname, value[ii], decim, tcomment, status);
        else
            ffpkyf(fptr, keyname, value[ii], decim, comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkne( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            float *value,       /* I - array of keyword values              */
            int decim,          /* I - number of decimals to display        */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential float values.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkye(fptr, keyname, value[ii], decim, tcomment, status);
        else
            ffpkye(fptr, keyname, value[ii], decim, comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpkng( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            double *value,      /* I - array of keyword values              */
            int decim,          /* I - number of decimals to display        */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed double values.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkyg(fptr, keyname, value[ii], decim, tcomment, status);
        else
            ffpkyg(fptr, keyname, value[ii], decim, comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpknd( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyroot,      /* I - root name of keywords to write       */
            int  nstart,        /* I - starting index number                */
            int  nkey,          /* I - number of keywords to write          */
            double *value,      /* I - array of keyword values              */
            int decim,          /* I - number of decimals to display        */
            char *comm[],       /* I - array of pointers to keyword comment */
            int  *status)       /* IO - error status                        */
/*
  Write (put) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential double values.
*/
{
    char keyname[FLEN_KEYWORD], tcomment[FLEN_COMMENT];
    int ii, jj, repeat, len;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;
    len = strlen(comm[0]);

    while (len > 0  && comm[0][len - 1] == ' ')
       len--;                               /* ignore trailing blanks */

    if (comm[0][len - 1] == '&')
    {
        len = minvalue(len, FLEN_COMMENT);
        tcomment[0] = '\0';
        strncat(tcomment, comm[0], len-1); /* don't copy the final '&' char */
        repeat = 1;
    }

    for (ii=0, jj=nstart; ii < nkey; ii++, jj++)
    {
        ffkeyn(keyroot, jj, keyname, status);
        if (repeat)
            ffpkyd(fptr, keyname, value[ii], decim, tcomment, status);
        else
            ffpkyd(fptr, keyname, value[ii], decim, comm[ii], status);

        if (*status > 0)
            return(*status);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffptdm( fitsfile *fptr, /* I - FITS file pointer                        */
            int colnum,     /* I - column number                            */
            int naxis,      /* I - number of axes in the data array         */
            long naxes[],   /* I - length of each data axis                 */
            int *status)    /* IO - error status                            */
/*
  write the TDIMnnn keyword describing the dimensionality of a column
*/
{
    char keyname[FLEN_KEYWORD], tdimstr[FLEN_VALUE], comm[FLEN_COMMENT];
    char value[80];
    int ii;

    if (*status > 0)
        return(*status);

    if (colnum < 1 || colnum > 999)
        return(*status = BAD_COL_NUM);

    if (naxis < 1)
        return(*status = BAD_DIMEN);

    ffkeyn("TDIM", colnum, keyname, status);      /* construct TDIMn name */

    strcpy(tdimstr, "(");            /* start constructing the TDIM value */   

    for (ii = 0; ii < naxis; ii++)
    {
        if (ii > 0)
            strcat(tdimstr, ",");   /* append the comma separator */

        if (naxes[ii] < 0)
            return(*status = BAD_NAXES);

        sprintf(value, "%ld", naxes[ii]);
        strcat(tdimstr, value);     /* append the axis size */
    }

    strcat(tdimstr, ")" );          /* append the closing parenthesis */

    strcpy(comm, "size of the multidimensional array");
    ffpkys(fptr, keyname, tdimstr, comm, status);  /* write the keyword */
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffphps( fitsfile *fptr, /* I - FITS file pointer                        */
            int bitpix,     /* I - number of bits per data value pixel      */
            int naxis,      /* I - number of axes in the data array         */
            long naxes[],   /* I - length of each data axis                 */
            int *status)    /* IO - error status                            */
/*
  write STANDARD set of required primary header keywords
*/
{
    int simple = 1;     /* does file conform to FITS standard? 1/0  */
    long pcount = 0;    /* number of group parameters (usually 0)   */
    long gcount = 1;    /* number of random groups (usually 1 or 0) */
    int extend = 1;     /* may FITS file have extensions?           */

    ffphpr(fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status);
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffphpr( fitsfile *fptr, /* I - FITS file pointer                        */
            int simple,     /* I - does file conform to FITS standard? 1/0  */
            int bitpix,     /* I - number of bits per data value pixel      */
            int naxis,      /* I - number of axes in the data array         */
            long naxes[],   /* I - length of each data axis                 */
            long pcount,    /* I - number of group parameters (usually 0)   */
            long gcount,    /* I - number of random groups (usually 1 or 0) */
            int extend,     /* I - may FITS file have extensions?           */
            int *status)    /* IO - error status                            */
/*
  write required primary header keywords
*/
{
    int ii;
    long longbitpix;
    char name[FLEN_KEYWORD], comm[FLEN_COMMENT], message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    if (fptr->headend != fptr->headstart[fptr->curhdu] )
        return(*status = HEADER_NOT_EMPTY);

    if (fptr->curhdu == 0)
    {                /* write primary array header */
        if (simple)
            strcpy(comm, "file does conform to FITS standard");
        else
            strcpy(comm, "file does not conform to FITS standard");

        ffpkyl(fptr, "SIMPLE", simple, comm, status);
    }
    else
    {               /* write IMAGE extension header */
        strcpy(comm, "IMAGE extension");
        ffpkys(fptr, "XTENSION", "IMAGE", comm, status);
    }

    longbitpix = bitpix;

    /* test for the 2 special cases that represent unsigned integers */
    if (longbitpix == USHORT_IMG)
        longbitpix = SHORT_IMG;
    else if (longbitpix == ULONG_IMG)
        longbitpix = LONG_IMG;

    if (longbitpix != BYTE_IMG && longbitpix != SHORT_IMG && 
        longbitpix != LONG_IMG &&
        longbitpix != FLOAT_IMG && longbitpix != DOUBLE_IMG)
    {
        sprintf(message,
        "Illegal value for BITPIX keyword: %d", bitpix);
        ffpmsg(message);
        return(*status = BAD_BITPIX);
    }

    strcpy(comm, "number of bits per data pixel");
    if (ffpkyj(fptr, "BITPIX", longbitpix, comm, status) > 0)
        return(*status);

    if (naxis < 0 || naxis > 999)
    {
        sprintf(message,
        "Illegal value for NAXIS keyword: %d", naxis);
        ffpmsg(message);
        return(*status = BAD_NAXIS);
    }

    strcpy(comm, "number of data axes");
    ffpkyj(fptr, "NAXIS", naxis, comm, status);

    strcpy(comm, "length of data axis ");
    for (ii = 0; ii < naxis; ii++)
    {
        if (naxes[ii] < 0)
        {
            sprintf(message,
            "Illegal value for NAXIS%d keyword: %ld", ii + 1,  naxes[ii]);
            ffpmsg(message);
            return(*status = BAD_NAXES);
        }

        sprintf(&comm[20], "%d", ii + 1);
        ffkeyn("NAXIS", ii + 1, name, status);
        ffpkyj(fptr, name, naxes[ii], comm, status);
    }

    if (fptr->curhdu == 0)  /* the primary array */
    {
        if (extend)
        {
            /* only write EXTEND keyword if value = true */
            strcpy(comm, "FITS dataset may contain extensions");
            ffpkyl(fptr, "EXTEND", extend, comm, status);
        }

        if (pcount < 0)
            return(*status = BAD_PCOUNT);

        else if (gcount < 1)
            return(*status = BAD_GCOUNT);

        else if (pcount > 0 || gcount > 1)
        {
            /* only write these keyword if non-standard values */
            strcpy(comm, "random group records are present");
            ffpkyl(fptr, "GROUPS", 1, comm, status);

            strcpy(comm, "number of random group parameters");
            ffpkyj(fptr, "PCOUNT", pcount, comm, status);
  
            strcpy(comm, "number of random groups");
            ffpkyj(fptr, "GCOUNT", gcount, comm, status);
        }

      /* write standard block of self-documentating comments */
      ffpcom(fptr,
      "FITS (Flexible Image Transport System) format defined in Astronomy and",
      status);

      ffpcom(fptr,
      "Astrophysics Supplement Series v44/p363, v44/p371, v73/p359, v73/p365.",
      status);

      ffpcom(fptr,
      "Contact the NASA Science Office of Standards and Technology for the",
      status);

      ffpcom(fptr,
      "FITS Definition document #100 and other FITS information.", status);
    }

    else  /* an IMAGE extension */

    {   /* image extension; cannot have random groups */
        if (pcount != 0)
            *status = BAD_PCOUNT;

        else if (gcount != 1)
            *status = BAD_GCOUNT;

        else
        {
            strcpy(comm, "required keyword; must = 0");
            ffpkyj(fptr, "PCOUNT", pcount, comm, status);
  
            strcpy(comm, "required keyword; must = 1");
            ffpkyj(fptr, "GCOUNT", gcount, comm, status);
        }
    }

    /* Write the BSCALE and BZERO keywords, if an unsigned integer image */
    if (bitpix == USHORT_IMG)
    {
        strcpy(comm, "offset data range to that of unsigned short");
        ffpkyg(fptr, "BZERO", 32768., 0, comm, status);
        strcpy(comm, "default scaling factor");
        ffpkyg(fptr, "BSCALE", 1.0, 0, comm, status);
    }
    else if (bitpix == ULONG_IMG)
    {
        strcpy(comm, "offset data range to that of unsigned long");
        ffpkyg(fptr, "BZERO", 2147483648., 0, comm, status);
        strcpy(comm, "default scaling factor");
        ffpkyg(fptr, "BSCALE", 1.0, 0, comm, status);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffphtb(fitsfile *fptr,  /* I - FITS file pointer                        */
           long naxis1,     /* I - width of row in the table                */
           long naxis2,     /* I - number of rows in the table              */
           int tfields,     /* I - number of columns in the table           */
           char **ttype,    /* I - name of each column                      */
           long *tbcol,     /* I - byte offset in row to each column        */
           char **tform,    /* I - value of TFORMn keyword for each column  */
           char **tunit,    /* I - value of TUNITn keyword for each column  */
           char *extnm,   /* I - value of EXTNAME keyword, if any         */
           int *status)     /* IO - error status                            */
/*
  Put required Header keywords into the ASCII TaBle:
*/
{
    int ii;
    char tfmt[30], name[FLEN_KEYWORD], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);
    else if (fptr->headend != fptr->headstart[fptr->curhdu] )
        return(*status = HEADER_NOT_EMPTY);
    else if (naxis1 < 0)
        return(*status = NEG_WIDTH);
    else if (naxis2 < 0)
        return(*status = NEG_ROWS);
    else if (tfields < 0 || tfields > 999)
        return(*status = BAD_TFIELDS);
    
    ffpkys(fptr, "XTENSION", "TABLE", "ASCII table extension", status);
    ffpkyj(fptr, "BITPIX", 8, "8-bit ASCII characters", status);
    ffpkyj(fptr, "NAXIS", 2, "2-dimensional ASCII table", status);
    ffpkyj(fptr, "NAXIS1", naxis1, "width of table in characters", status);
    ffpkyj(fptr, "NAXIS2", naxis2, "number of rows in table", status);
    ffpkyj(fptr, "PCOUNT", 0, "no group parameters (required keyword)", status);
    ffpkyj(fptr, "GCOUNT", 1, "one data group (required keyword)", status);
    ffpkyj(fptr, "TFIELDS", tfields, "number of fields in each row", status);

    for (ii = 0; ii < tfields; ii++) /* loop over every column */
    {
        if ( *(ttype[ii]) )  /* optional TTYPEn keyword */
        {
          sprintf(comm, "label for field %3d", ii + 1);
          ffkeyn("TTYPE", ii + 1, name, status);
          ffpkys(fptr, name, ttype[ii], comm, status);
        }

        if (tbcol[ii] < 1 || tbcol[ii] > naxis1)
           *status = BAD_TBCOL;

        sprintf(comm, "beginning column of field %3d", ii + 1);
        ffkeyn("TBCOL", ii + 1, name, status);
        ffpkyj(fptr, name, tbcol[ii], comm, status);

        strcpy(tfmt, tform[ii]);  /* required TFORMn keyword */
        ffupch(tfmt);
        ffkeyn("TFORM", ii + 1, name, status);
        ffpkys(fptr, name, tfmt, "Fortran-77 format of field", status);

        if ( *(tunit[ii]) )       /* optional TUNITn keyword */
        {
          ffkeyn("TUNIT", ii + 1, name, status);
          ffpkys(fptr, name, tunit[ii], "physical unit of field", status) ;
        }

        if (*status > 0)
            break;       /* abort loop on error */
    }

    if (extnm[0])       /* optional EXTNAME keyword */
        ffpkys(fptr, "EXTNAME", extnm,
               "name of this ASCII table extension", status);

    if (*status > 0)
        ffpmsg("Failed to write ASCII table header keywords (ffphtb)");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffphbn(fitsfile *fptr,  /* I - FITS file pointer                        */
           long naxis2,     /* I - number of rows in the table              */
           int tfields,     /* I - number of columns in the table           */
           char **ttype,    /* I - name of each column                      */
           char **tform,    /* I - value of TFORMn keyword for each column  */
           char **tunit,    /* I - value of TUNITn keyword for each column  */
           char *extnm,   /* I - value of EXTNAME keyword, if any         */
           long pcount,     /* I - size of the variable length heap area    */
           int *status)     /* IO - error status                            */
/*
  Put required Header keywords into the Binary Table:
*/
{
    int ii, datatype;
    long repeat, width, naxis1;

    char tfmt[30], name[FLEN_KEYWORD], comm[FLEN_COMMENT];
    char *cptr;

    if (*status > 0)
        return(*status);
    else if (fptr->headend != fptr->headstart[fptr->curhdu] )
        return(*status = HEADER_NOT_EMPTY);
    else if (naxis2 < 0)
        return(*status = NEG_ROWS);
    else if (pcount < 0)
        return(*status = BAD_PCOUNT);
    else if (tfields < 0 || tfields > 999)
        return(*status = BAD_TFIELDS);

    ffpkys(fptr, "XTENSION", "BINTABLE", "binary table extension", status);
    ffpkyj(fptr, "BITPIX", 8, "8-bit bytes", status);
    ffpkyj(fptr, "NAXIS", 2, "2-dimensional binary table", status);

    naxis1 = 0;
    for (ii = 0; ii < tfields; ii++)  /* sum the width of each field */
    {
        ffbnfm(tform[ii], &datatype, &repeat, &width, status);

        if (datatype == TSTRING)
            naxis1 += repeat;   /* one byte per char */
        else if (datatype == TBIT)
            naxis1 += (repeat + 7) / 8;
        else if (datatype > 0)
            naxis1 += repeat * (datatype / 10);
        else   /* this is a variable length descriptor (neg. datatype) */
            naxis1 += 8;

        if (*status > 0)
            break;       /* abort loop on error */
    }

    ffpkyj(fptr, "NAXIS1", naxis1, "width of table in bytes", status);
    ffpkyj(fptr, "NAXIS2", naxis2, "number of rows in table", status);

    /*
      the initial value of PCOUNT (= size of the variable length array heap)
      should always be zero.  If any variable length data is written, then
      the value of PCOUNT will be updated when the HDU is closed
    */
    ffpkyj(fptr, "PCOUNT", 0, "size of special data area", status);
    ffpkyj(fptr, "GCOUNT", 1, "one data group (required keyword)", status);
    ffpkyj(fptr, "TFIELDS", tfields, "number of fields in each row", status);

    for (ii = 0; ii < tfields; ii++) /* loop over every column */
    {
        if ( *(ttype[ii]) )  /* optional TTYPEn keyword */
        {
          sprintf(comm, "label for field %3d", ii + 1);
          ffkeyn("TTYPE", ii + 1, name, status);
          ffpkys(fptr, name, ttype[ii], comm, status);
        }

        strcpy(tfmt, tform[ii]);  /* required TFORMn keyword */
        ffupch(tfmt);

        ffkeyn("TFORM", ii + 1, name, status);
        strcpy(comm, "data format of field");

        ffbnfm(tfmt, &datatype, &repeat, &width, status);

        if (datatype == TSTRING)
            strcat(comm, ": ASCII Character");
        else if (datatype == TBIT)
           strcat(comm, ": BIT");
        else if (datatype == TBYTE)
           strcat(comm, ": BYTE");
        else if (datatype == TLOGICAL)
           strcat(comm, ": 1-byte LOGICAL");
        else if (datatype == TSHORT)
           strcat(comm, ": 2-byte INTEGER");
        else if (datatype == TUSHORT)
           strcat(comm, ": 2-byte INTEGER");
        else if (datatype == TLONG)
           strcat(comm, ": 4-byte INTEGER");
        else if (datatype == TULONG)
           strcat(comm, ": 4-byte INTEGER");
        else if (datatype == TFLOAT)
           strcat(comm, ": 4-byte REAL");
        else if (datatype == TDOUBLE)
           strcat(comm, ": 8-byte DOUBLE");
        else if (datatype == TCOMPLEX)
           strcat(comm, ": COMPLEX");
        else if (datatype == TDBLCOMPLEX)
           strcat(comm, ": DOUBLE COMPLEX");
        else if (datatype < 0)
           strcat(comm, ": variable length array");

        if (datatype == TUSHORT) 
        {
           /* Replace the 'U' with an 'I' in the TFORMn code */
           cptr = tfmt;
           while (*cptr != 'U') 
              cptr++;

           *cptr = 'I';
           ffpkys(fptr, name, tfmt, comm, status);

           /* write the TZEROn and TSCALn keywords */
           ffkeyn("TZERO", ii + 1, name, status);
           strcpy(comm, "offset for unsigned integers");

           ffpkyg(fptr, name, 32768., 0, comm, status);

           ffkeyn("TSCAL", ii + 1, name, status);
           strcpy(comm, "data are not scaled");
           ffpkyg(fptr, name, 1., 0, comm, status);
        }
        else if (datatype == TULONG) 
        {
           /* Replace the 'V' with an 'J' in the TFORMn code */
           cptr = tfmt;
           while (*cptr != 'V') 
              cptr++;

           *cptr = 'J';
           ffpkys(fptr, name, tfmt, comm, status);

           /* write the TZEROn and TSCALn keywords */
           ffkeyn("TZERO", ii + 1, name, status);
           strcpy(comm, "offset for unsigned integers");

           ffpkyg(fptr, name, 2147483648., 0, comm, status);

           ffkeyn("TSCAL", ii + 1, name, status);
           strcpy(comm, "data are not scaled");
           ffpkyg(fptr, name, 1., 0, comm, status);
        }
        else
        {
           ffpkys(fptr, name, tfmt, comm, status);
        }

        if ( *(tunit[ii]) )       /* optional TUNITn keyword */
        {
          ffkeyn("TUNIT", ii + 1, name, status);
          ffpkys(fptr, name, tunit[ii],
             "physical unit of field", status);
        }

        if (*status > 0)
            break;       /* abort loop on error */
    }

    if (extnm[0])       /* optional EXTNAME keyword */
        ffpkys(fptr, "EXTNAME", extnm,
               "name of this binary table extension", status);

    if (*status > 0)
        ffpmsg("Failed to write binary table header keywords (ffphbn)");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffi2c(long ival,   /* I - value to be converted to a string */
          char *cval,  /* O - character string representation of the value */
          int *status) /* IO - error status */
/*
  convert  value to a null-terminated formatted string.
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    cval[0] = '\0';

    if (sprintf(cval, "%ld", ival) < 0)
    {
        ffpmsg("Error in ffi2c converting integer to string");
        *status = BAD_I2C;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffl2c(int lval,    /* I - value to be converted to a string */
          char *cval,  /* O - character string representation of the value */
          int *status) /* IO - error status ) */
/*
  convert logical value to a null-terminated formatted string.  If the
  input value == 0, then the output character is the letter F, else
  the output character is the letter T.  The output string is null terminated.
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (lval)
        strcpy(cval,"T");
    else
        strcpy(cval,"F");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffs2c(char *instr,   /* I - null terminated input string  */
          char *outstr,  /* O - null terminated quoted output string */
          int *status)   /* IO - error status */
/*
  convert an input string to a quoted string. Leading spaces 
  are significant.  FITS string keyword values must be at least 
  8 chars long so pad out string with spaces if necessary.
      Example:   km/s ==> 'km/s    '
  Single quote characters in the input string will be replace by
  two single quote characters. e.g., o'brian ==> 'o''brian'
*/
{
    size_t len, ii, jj;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    outstr[0] = '\'';      /* start output string with a quote */

    len = strlen(instr);
    if (len > 68)
        len = 68;    /* limit input string to 68 chars */

    for (ii=0, jj=1; ii < len && jj < 69; ii++, jj++)
    {
        outstr[jj] = instr[ii];  /* copy each char from input to output */
        if (instr[ii] == '\'')
        {
            jj++;
            outstr[jj]='\'';   /* duplicate any apostrophies in the input */
        }
    }

    for (; jj < 9; jj++)       /* pad string so it is at least 8 chars long */
        outstr[jj] = ' ';

    if (jj == 70)   /* only occurs if the last char of string was a quote */
        outstr[69] = '\0';
    else
    {
        outstr[jj] = '\'';         /* append closing quote character */
        outstr[jj+1] = '\0';          /* terminate the string */
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffr2f(float fval,   /* I - value to be converted to a string */
          int  decim,   /* I - number of decimal places to display */
          char *cval,   /* O - character string representation of the value */
          int  *status) /* IO - error status */
/*
  convert float value to a null-terminated F format string
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    cval[0] = '\0';

    if (decim < 0)
    {
        ffpmsg("Error in ffr2f:  no. of decimal places < 0");
        return(*status = BAD_DECIM);
    }

    if (sprintf(cval, "%.*f", decim, fval) < 0)
    {
        ffpmsg("Error in ffr2f converting float to string");
        *status = BAD_F2C;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffr2e(float fval,  /* I - value to be converted to a string */
         int decim,    /* I - number of decimal places to display */
         char *cval,   /* O - character string representation of the value */
         int *status)  /* IO - error status */
/*
  convert float value to a null-terminated exponential format string
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    cval[0] = '\0';

    if (decim < 0)
    {
        ffpmsg("Error in ffr2e:  no. of decimal places < 0");
        return(*status = BAD_DECIM);
    }

    if (sprintf(cval, "%.*E", decim, fval) < 0)
    {
        ffpmsg("Error in ffr2e converting float to string");
        *status = BAD_F2C;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffd2f(double dval,  /* I - value to be converted to a string */
          int decim,    /* I - number of decimal places to display */
          char *cval,   /* O - character string representation of the value */
          int *status)  /* IO - error status */
/*
  convert double value to a null-terminated F format string
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    cval[0] = '\0';

    if (decim < 0)
    {
        ffpmsg("Error in ffd2f:  no. of decimal places < 0");
        return(*status = BAD_DECIM);
    }

    if (sprintf(cval, "%.*f", decim, dval) < 0)
    {
        ffpmsg("Error in ffd2f converting double to string");
        *status = BAD_F2C;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffd2e(double dval,  /* I - value to be converted to a string */
          int decim,    /* I - number of decimal places to display */
          char *cval,   /* O - character string representation of the value */
          int *status)  /* IO - error status */
/*
  convert double value to a null-terminated exponential format string.
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    cval[0] = '\0';

    if (decim < 0)
    {
        ffpmsg("Error in ffd2e:  no. of decimal places < 0");
        return(*status = BAD_DECIM);
    }

    if (sprintf(cval, "%.*E", decim, dval) < 0)
    {
        ffpmsg("Error in ffd2e converting double to string");
        *status = BAD_F2C;
    }

    return(*status);
}

