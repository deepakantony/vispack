/*  This file, getkey.c, contains routines that read keywords from         */
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
#include "fitsio2.h"

/*--------------------------------------------------------------------------*/
int ffghsp(fitsfile *fptr,  /* I - FITS file pointer                     */
           int *nexist,     /* O - number of existing keywords in header */
           int *nmore,      /* O - how many more keywords will fit       */
           int *status)     /* IO - error status                         */
/*
  returns the number of existing keywords (not counting the END keyword)
  and the number of more keyword that will fit in the current header 
  without having to insert more FITS blocks.
*/
{
    *nexist = ( (fptr->headend) - (fptr->headstart[fptr->curhdu]) ) / 80;

    if (fptr->datastart == DATA_UNDEFINED)
        *nmore = -1;   /* data not written yet, so room for any keywords */
    else
        /* calculate space available between the data and the END card */
        *nmore = (fptr->datastart - fptr->headend) / 80 - 1;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffghps(fitsfile *fptr, /* I - FITS file pointer                     */
          int *nexist,     /* O - number of existing keywords in header */
          int *position,   /* O - position of next keyword to be read   */
          int *status)     /* IO - error status                         */
/*
  return the number of existing keywords and the position of the next
  keyword that will be read.
*/
{
  *nexist = ( (fptr->headend) - (fptr->headstart[fptr->curhdu]) ) / 80;
  *position = ( (fptr->nextkey) - (fptr->headstart[fptr->curhdu]) ) / 80 + 1;
  return(*status);
}
/*--------------------------------------------------------------------------*/
int ffmaky(fitsfile *fptr,    /* I - FITS file pointer                    */
          int nrec,           /* I - one-based keyword number to move to  */
          int *status)        /* IO - error status                        */
{
/*
  move pointer to the specified absolute keyword position.  E.g. this keyword 
  will then be read by the next call to ffgnky.
*/

    fptr->nextkey = fptr->headstart[fptr->curhdu] + ( (nrec - 1) * 80);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffmrky(fitsfile *fptr,    /* I - FITS file pointer                   */
          int nmove,          /* I - relative number of keywords to move */
          int *status)        /* IO - error status                       */
{
/*
  move pointer to the specified keyword position relative to the current
  position.  E.g. this keyword  will then be read by the next call to ffgnky.
*/

    int absrec;

    absrec = ( (fptr->nextkey) - (fptr->headstart[fptr->curhdu]) ) / 80
              + 1 + nmove;

    ffmaky(fptr, absrec, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgnky(fitsfile *fptr,  /* I - FITS file pointer     */
           char *card,      /* O - card string           */
           int *status)     /* IO - error status         */
/*
  read the next keyword from the header - used internally by cfitsio
*/
{
    int jj, nrec;
    long bytepos, endhead;
    char message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    card[0] = '\0';  /* make sure card is terminated, even affer read error */

/*
  Check that nextkey points to a legal keyword position.  Note that headend
  is the current end of the header, i.e., the position where a new keyword
  would be appended, however, if there are more than 1 FITS block worth of
  blank keywords at the end of the header (36 keywords per 2880 byte block)
  then the actual physical END card must be located at a starting position
  which is just 2880 bytes prior to the start of the data unit.
*/

    bytepos = fptr->nextkey;
    endhead = maxvalue( (fptr->headend), (fptr->datastart - 2880) );

    if (bytepos > endhead ||            /* nextkey must be < endhead and */
        bytepos < fptr->headstart[fptr->curhdu] )    /* > than headstart */
    {
        nrec = (bytepos - fptr->headstart[fptr->curhdu]) / 80 + 1;
        sprintf(message, "Cannot get keyword number %d.  It does not exist.",
                nrec);
        ffpmsg(message);
        return(*status = KEY_OUT_BOUNDS);
    }
      
    ffmbyt(fptr, bytepos, REPORT_EOF, status);  /* move to read pos. */

    card[80] = '\0';  /* make sure card is terminate, even if ffgbyt fails */

    if (ffgbyt(fptr, 80, card, status) <= 0) 
    {
        ffmrky(fptr, 1, status);    /* increment pointer to next keyword */

        for (jj=79; jj >= 0; jj--)  /* replace trailing blanks with nulls */
        {
            if (card[jj] == ' ')
                card[jj] = '\0';
            else
                break;
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgnxk( fitsfile *fptr,     /* I - FITS file pointer              */
            char **inclist,     /* I - list of included keyword names */
            int ninc,           /* I - number of names in inclist     */
            char **exclist,     /* I - list of excluded keyword names */
            int nexc,           /* I - number of names in exclist     */
            char *card,         /* O - first matching keyword         */
            int  *status)       /* IO - error status                  */
/*
    Return the next keyword that matches one of the names in inclist
    but does not match any of the names in exclist.  The search
    goes from the current position to the end of the header, only.
    Wild card characters may be used in the name lists ('*', '?' and '#').
*/
{
    int casesn, match, exact;
    long ii, jj;
    char keybuf[FLEN_CARD];

    card[0] = '\0';
    if (*status > 0)
        return(*status);

    casesn = FALSE;

    /* get next card, and return with an error if hit end of header */
    while( ffgcrd(fptr, "*", keybuf, status) <= 0)
    {
        /* does keyword match any names in the include list? */
        for (ii = 0; ii < ninc; ii++)
        {
            ffcmps(inclist[ii], keybuf, casesn, &match, &exact);
            if (match)
            {
                /* does keyword match any names in the exclusion list? */
                for (jj = 0; jj < nexc; jj++)
                {
                    ffcmps(exclist[jj], keybuf, casesn, &match, &exact);
                    if (match)
                        break;
                }

                if (jj >= nexc)
                {
                    /* not in exclusion list, so return this keyword */
                    strcat(card, keybuf);
                    return(*status);
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgky( fitsfile *fptr,     /* I - FITS file pointer        */
           int  datatype,      /* I - datatype of the value    */
           char *keyname,      /* I - name of keyword to read  */
           void *value,        /* O - keyword value            */
           char *comm,         /* O - keyword comment          */
           int  *status)       /* IO - error status            */
/*
  Read (gut) the keyword value and comment from the FITS header.
  Reads a keyword value with the datatype specified by the 2nd argument.
*/
{
    long longval;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (datatype == TSTRING)
    {
        ffgkys(fptr, keyname, (char *) value, comm, status);
    }
    else if (datatype == TBYTE)
    {
        ffgkyj(fptr, keyname, &longval, comm, status);
        *(unsigned char *) value = longval;
    }
    else if (datatype == TUSHORT)
    {
        ffgkyj(fptr, keyname, &longval, comm, status);
        *(unsigned short *) value = longval;
    }
    else if (datatype == TSHORT)
    {
        ffgkyj(fptr, keyname, &longval, comm, status);
        *(short *) value = longval;
    }
    else if (datatype == TINT)
    {
        ffgkyj(fptr, keyname, &longval, comm, status);
        *(int *) value = longval;
    }
    else if (datatype == TLOGICAL)
    {
        ffgkyl(fptr, keyname, (int *) value, comm, status);
    }
    else if (datatype == TULONG)
    {
        ffgkyj(fptr, keyname, &longval, comm, status);
        *(unsigned long *) value = longval;
    }
    else if (datatype == TLONG)
    {
        ffgkyj(fptr, keyname, (long *) value, comm, status);
    }
    else if (datatype == TFLOAT)
    {
        ffgkye(fptr, keyname, (float *) value, comm, status);
    }
    else if (datatype == TDOUBLE)
    {
        ffgkyd(fptr, keyname, (double *) value, comm, status);
    }
    else
        *status = BAD_DATATYPE;

    return(*status);
} 
/*--------------------------------------------------------------------------*/
int ffgkey( fitsfile *fptr,     /* I - FITS file pointer        */
            char *keyname,      /* I - name of keyword to read  */
            char *keyval,       /* O - keyword value            */
            char *comm,         /* O - keyword comment          */
            int  *status)       /* IO - error status            */
/*
  Read (get) the named keyword, returning the keyword value and comment.
  The value is just the literal string of characters in the value field
  of the keyword.  In the case of a string valued keyword, the returned
  value includes the leading and closing quote characters.  The value may be
  up to 70 characters long, and the comment may be up to 72 characters long.
  If the keyword has no value (no equal sign in column 9) then a null value
  is returned.
*/
{
    char card[FLEN_CARD];

    keyval[0] = '\0';
    if (comm)
       comm[0] = '\0';

    if (ffgcrd(fptr, keyname, card, status) > 0)    /* get the 80-byte card */
        return(*status);

    ffpsvc(card, keyval, comm, status);      /* parse the value and comment */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgrec( fitsfile *fptr,     /* I - FITS file pointer          */
            int nrec,           /* I - number of keyword to read  */
            char *card,         /* O - keyword card               */
            int  *status)       /* IO - error status              */
/*
  Read (get) the nrec-th keyword, returning the entire keyword card up to
  80 characters long.  The first keyword in the header has nrec = 1, not 0.
  The returned card value is null terminated with any trailing blank 
  characters removed.  If nrec = 0, then this routine simply moves the
  current header pointer to the top of the header.
*/
{
    char sbuff[FLEN_CARD];

    if (*status > 0)
        return(*status);

    if (nrec == 0)
    {
        ffmaky(fptr, 1, status);  /* simply move to beginning of header */
        card[0] = '\0';           /* and return null card */
    }
    else if (nrec > 0)
    {
        ffmaky(fptr, nrec, status);
        ffgnky(fptr, card, status);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcrd( fitsfile *fptr,     /* I - FITS file pointer        */
            char *name,         /* I - name of keyword to read  */
            char *card,         /* O - keyword card             */
            int  *status)       /* IO - error status            */
/*
  Read (get) the named keyword, returning the entire keyword card up to
  80 characters long.  The first keyword in the header has nrec = 1, not 0.
  The returned card value is null terminated with any trailing blank 
  characters removed.

  If the input name contains wild cards ('?' matches any single char
  and '*' matches any sequence of chars, # matches any string of decimal
  digits) then the search ends once the end of header is reached and does 
  not automatically resume from the top of the header.
*/
{
    int nkeys, nextkey, ntodo, namelen, ii, jj, wild, match, exact;
    char keyname[10], ctemp;

    if (*status > 0)
        return(*status);

    /* does input name contain wild card chars?  ('?',  '*', or '#') */
    if (strchr(name,'?') || strchr(name,'*') || strchr(name,'#'))
        wild = 1;
    else
        wild = 0;

    strncpy(keyname, name, 8);
    keyname[8] = '\0';  /* make sure string is terminated */

    namelen=strlen(keyname);

    for (ii=0; ii < namelen; ii++)       
        keyname[ii] = toupper(name[ii]);    /*  make sure upper case  */

    for (ii=namelen; ii < 8; ii++)
        keyname[ii] = ' '; /*  pad name with blanks to 8 characters  */

    ffghps(fptr, &nkeys, &nextkey, status); /* get no. keywords and position */

    ntodo = nkeys - nextkey + 1;  /* first, read from next keyword to end */
    for (jj=0; jj < 2; jj++)
    {
        for (ii = 0; ii < ntodo; ii++)
        {
            ffgnky(fptr, card, status);     /*  get next keyword */
            if (wild)
            {
                ctemp = card[8];   /* save 9th char in temporary variable */
                card[8] = '\0';    /* terminate the keyword name */
                ffcmps(keyname, card, 1, &match, &exact);
                if (match)
                {
                    card[8] = ctemp;  /* restore the 9th char */
                    return(*status); /* found a matching keyword */
                }
            }
            else if (strncmp(keyname, card, 8) == 0)
                    return(*status);  /* found the matching keyword */
        }

        if (wild || jj == 1)
            break;  /* stop at end of header if template contains wildcards */

        ffmaky(fptr, 1, status);  /* reset pointer to beginning of header */
        ntodo = nextkey - 1;      /* number of keyword to read */ 
    }

    return(*status = KEY_NO_EXIST);  /* couldn't find the keyword */
}
/*--------------------------------------------------------------------------*/
int ffgkys( fitsfile *fptr,     /* I - FITS file pointer         */
            char *keyname,      /* I - name of keyword to read   */
            char *value,        /* O - keyword value             */
            char *comm,         /* O - keyword comment           */
            int  *status)       /* IO - error status             */
/*
  Get KeYword with a String value:
  Read (get) a simple string valued keyword.  The returned value may be up to 
  68 chars long ( + 1 null terminator char).  The routine does not support the
  HEASARC convention for continuing long string values over multiple keywords.
  The ffgkls routine may be used to read long continued strings. The returned
  comment string may be up to 69 characters long (including null terminator).
*/
{
    char valstring[FLEN_VALUE];

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */
    ffc2s(valstring, value, status);   /* remove quotes from string */
 
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkls( fitsfile *fptr,     /* I - FITS file pointer         */
            char *keyname,      /* I - name of keyword to read   */
            char **value,       /* O - pointer to keyword value  */
            char *comm,         /* O - keyword comment           */
            int  *status)       /* IO - error status             */
/*
  Get Keyword with possible Long String value:
  Read (get) the named keyword, returning the value and comment.
  The returned value string may be arbitrarily long (by using the HEASARC
  convention for continuing long string values over multiple keywords) so
  this routine allocates the required memory for the returned string value.
  It is up to the calling routine to free the memory once it is finished
  with the value string.  The returned comment string may be up to 69
  characters long.
*/
{
    char valstring[FLEN_VALUE];
    int contin;
    size_t len;

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */

    /* allocate space, minus the 2 quote chars, plus 1 for null */
    *value = (char *) malloc(strlen(valstring)-1);

    ffc2s(valstring, *value, status);   /* convert string to value */
    len = strlen(*value);

    /* If last character is a & then value may be continued on next keyword */
    contin = 1;
    while (contin)  
    {
        if (*(*value+len-1) == '&')  /*  is last char an anpersand?  */
        {
            ffgcnt(fptr, valstring, status);
            if (valstring)    /* a null valstring indicates no continuation */
            {
               *(*value+len-1) = '\0';         /* erase the trailing & char */
               len += strlen(valstring) - 1;
               *value = (char *) realloc(*value, len); /* increase str size */
               strcat(*value, valstring);     /* append the continued chars */
            }
            else
                contin = 0;
        }
        else
            contin = 0;
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcnt( fitsfile *fptr,     /* I - FITS file pointer         */
            char *value,        /* O - continued string value    */
            int  *status)       /* IO - error status             */
/*
  Attempt to read the next keyword, returning the string value
  if it is a continuation of the previous string keyword value.
  This uses the HEASARC convention for continuing long string values
  over multiple keywords.  Each continued string is terminated with a
  backslash character, and the continuation follows on the next keyword
  which must have the name CONTINUE without an equal sign in column 9
  of the card.  If the next card is not a continuation, then the returned
  value string will be null.
*/
{
    int tstatus;
    char card[FLEN_CARD], strval[FLEN_VALUE], comm[FLEN_COMMENT];

    tstatus = 0;
    value[0] = '\0';

    if (ffgnky(fptr, card, &tstatus) > 0)  /*  read next keyword  */
        return(*status);                   /*  hit end of header  */

    if (strncmp(card, "CONTINUE  ", 10) == 0)  /* a continuation card? */
    {
        strncpy(card, "D2345678=  ", 10); /* overwrite a dummy keyword name */
        ffpsvc(card, strval, comm, &tstatus);  /*  get the string value  */
        ffc2s(strval, value, &tstatus);    /* remove the surrounding quotes */

        if (tstatus)       /*  return null if error status was returned  */
           value[0] = '\0';
    }
    else
        ffmrky(fptr, -1, status);  /* reset the keyword pointer */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkyl( fitsfile *fptr,     /* I - FITS file pointer         */
            char *keyname,      /* I - name of keyword to read   */
            int  *value,        /* O - keyword value             */
            char *comm,         /* O - keyword comment           */
            int  *status)       /* IO - error status             */
/*
  Read (get) the named keyword, returning the value and comment.
  The returned value = 1 if the keyword is true, else = 0 if false.
  The comment may be up to 69 characters long.
*/
{
    char valstring[FLEN_VALUE];

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */
    ffc2l(valstring, value, status);   /* convert string to value */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkyj( fitsfile *fptr,     /* I - FITS file pointer         */
            char *keyname,      /* I - name of keyword to read   */
            long *value,        /* O - keyword value             */
            char *comm,         /* O - keyword comment           */
            int  *status)       /* IO - error status             */
/*
  Read (get) the named keyword, returning the value and comment.
  The value will be implicitly converted to a (long) integer if it not
  already of this datatype.  The comment may be up to 69 characters long.
*/
{
    char valstring[FLEN_VALUE];

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */
    ffc2i(valstring, value, status);   /* convert string to value */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkye( fitsfile *fptr,     /* I - FITS file pointer         */
            char  *keyname,     /* I - name of keyword to read   */
            float *value,       /* O - keyword value             */
            char  *comm,        /* O - keyword comment           */
            int   *status)      /* IO - error status             */
/*
  Read (get) the named keyword, returning the value and comment.
  The value will be implicitly converted to a float if it not
  already of this datatype.  The comment may be up to 69 characters long.
*/
{
    char valstring[FLEN_VALUE];

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */
    ffc2r(valstring, value, status);   /* convert string to value */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkyd( fitsfile *fptr,      /* I - FITS file pointer         */
            char   *keyname,     /* I - name of keyword to read   */
            double *value,       /* O - keyword value             */
            char   *comm,        /* O - keyword comment           */
            int    *status)      /* IO - error status             */
/*
  Read (get) the named keyword, returning the value and comment.
  The value will be implicitly converted to a double if it not
  already of this datatype.  The comment may be up to 69 characters long.
*/
{
    char valstring[FLEN_VALUE];

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */
    ffc2d(valstring, value, status);   /* convert string to value */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkyt( fitsfile *fptr,      /* I - FITS file pointer                 */
            char   *keyname,     /* I - name of keyword to read           */
            long   *ivalue,      /* O - integer part of keyword value     */
            double *fraction,    /* O - fractional part of keyword value  */
            char   *comm,        /* O - keyword comment                   */
            int    *status)      /* IO - error status                     */
/*
  Read (get) the named keyword, returning the value and comment.
  The integer and fractional parts of the value are returned in separate
  variables, to allow more numerical precision to be passed.  This
  effectively passes a 'triple' precision value, with a 4-byte integer
  and an 8-byte fraction.  The comment may be up to 69 characters long.
*/
{
    char valstring[FLEN_VALUE];
    char *loc;

    ffgkey(fptr, keyname, valstring, comm, status);  /* read the keyword */

    /*  read the entire value string as a double, to get the integer part */
    ffc2d(valstring, fraction, status);

    *ivalue = (long) *fraction;

    *fraction = *fraction - *ivalue;

    /* see if we need to read the fractional part again with more precision */
    /* look for decimal point, without an exponential E or D character */

    loc = strchr(valstring, '.');
    if (loc)
    {
        if (!strchr(valstring, 'E') && !strchr(valstring, 'D'))
            ffc2d(loc, fraction, status);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkyn( fitsfile *fptr,      /* I - FITS file pointer             */
            int    nkey,         /* I - number of the keyword to read */
            char   *keyname,     /* O - name of the keyword           */
            char   *value,       /* O - keyword value                 */
            char   *comm,        /* O - keyword comment               */
            int    *status)      /* IO - error status                 */
/*
  Read (get) the nkey-th keyword returning the keyword name, value and comment.
  The value is just the literal string of characters in the value field
  of the keyword.  In the case of a string valued keyword, the returned
  value includes the leading and closing quote characters.  The value may be
  up to 70 characters long, and the comment may be up to 72 characters long.
  If the keyword has no value (no equal sign in column 9) then a null value
  is returned.  If comm = NULL, then do not return the comment string.
*/
{
    char card[FLEN_CARD], sbuff[FLEN_CARD];
    int ii;

    keyname[0] = '\0';
    value[0] = '\0';
    if (comm)
        comm[0] = '\0';

    if (ffgrec(fptr, nkey, card, status) > 0 )  /* get the 80-byte card */
        return(*status);

    strncpy(keyname,card,8);  /* first 8 characters = the keyword name */
    keyname[8] = '\0';

    for (ii=7; ii >= 0; ii--)  /* replace trailing blanks with nulls */
    {
        if (keyname[ii] == ' ')
            keyname[ii] = '\0';
        else
            break;
    }

    if (ffpsvc(card, value, comm, status) > 0)   /* parse value and comment */
        return(*status);

    if (fftkey(keyname, status) > 0)  /* test keyword name; catches no END */
    {
     sprintf(sbuff,"Name of keyword no. %d contains illegal character(s): %s",
              nkey, keyname);
     ffpmsg(sbuff);

     if (nkey % 36 == 0)  /* test if at beginning of 36-card FITS record */
            ffpmsg("  (This may indicate a missing END keyword).");
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkns( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyname,      /* I - root name of keywords to read        */
            int  nstart,        /* I - starting index number                */
            int  nmax,          /* I - maximum number of keywords to return */
            char *value[],      /* O - array of pointers to keyword values  */
            int  *nfound,       /* O - number of values that were returned  */
            int  *status)       /* IO - error status                        */
/*
  Read (get) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NMAX -1) inclusive.  
  This routine does NOT support the HEASARC long string convention.
*/
{
    int nend, lenroot, ii, nkeys, mkeys, tstatus;
    long ival;
    char keyroot[FLEN_KEYWORD], keyindex[8], card[FLEN_CARD];
    char svalue[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = '\0';
    strncat(keyroot, keyname, 8);
     
    lenroot = strlen(keyroot);
    if (lenroot == 0 || lenroot > 7)     /*  root must be 1 - 7 chars long  */
        return(*status);

    for (ii=0; ii < lenroot; ii++)           /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);

    ffghps(fptr, &nkeys, &mkeys, status);  /*  get the number of keywords  */

    for (ii=3; ii <= nkeys; ii++)  
    {
       if (ffgrec(fptr, ii, card, status) > 0)     /*  get next keyword  */
           return(*status);

       if (strncmp(keyroot, card, lenroot) == 0)  /* see if keyword matches */
       {
           keyindex[0] = '\0';
           strncat(keyindex, &card[lenroot], 8-lenroot);  /*  copy suffix */

           tstatus = 0;
           if (ffc2ii(keyindex, &ival, &tstatus) <= 0)     /*  test suffix  */
           {
               if (ival <= nend && ival >= nstart)
               {
                 ffpsvc(card, svalue, comm, status);   /*  parse the value */

                 ffc2s(svalue, value[ival-nstart], status); /* convert */
                 if (ival > *nfound)
                       *nfound = ival;  /* record the max index found */ 
               }
           }
       }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgknl( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyname,      /* I - root name of keywords to read        */
            int  nstart,        /* I - starting index number                */
            int  nmax,          /* I - maximum number of keywords to return */
            int  *value,        /* O - array of keyword values              */
            int  *nfound,       /* O - number of values that were returned  */
            int  *status)       /* IO - error status                        */
/*
  Read (get) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NMAX -1) inclusive.  
  The returned value = 1 if the keyword is true, else = 0 if false.
*/
{
    int nend, lenroot, ii, nkeys, mkeys, tstatus;
    long ival;
    char keyroot[FLEN_KEYWORD], keyindex[8], card[FLEN_CARD];
    char svalue[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = '\0';
    strncat(keyroot, keyname, 8);

    lenroot = strlen(keyroot);
    if (lenroot == 0 || lenroot > 7)     /*  root must be 1 - 7 chars long  */
        return(*status);

    for (ii=0; ii < lenroot; ii++)           /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);

    ffghps(fptr, &nkeys, &mkeys, status);  /*  get the number of keywords  */

    ffmaky(fptr, 3, status);  /* move to 3rd keyword (skip 1st 2 keywords) */

    for (ii=3; ii <= nkeys; ii++)  
    {
       if (ffgnky(fptr, card, status) > 0)     /*  get next keyword  */
           return(*status);

       if (strncmp(keyroot, card, lenroot) == 0)  /* see if keyword matches */
       {
           keyindex[0] = '\0';
           strncat(keyindex, &card[lenroot], 8-lenroot);  /*  copy suffix */

           tstatus = 0;
           if (ffc2ii(keyindex, &ival, &tstatus) <= 0)    /*  test suffix  */
               if (ival <= nend && ival >= nstart)
               {
                 ffpsvc(card, svalue, comm, status);   /*  parse the value */
                 ffc2l(svalue, &value[ival-nstart], status); /* convert*/
                 if (ival > *nfound)
                       *nfound = ival;   /* record the max index found */
               }
       }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgknj( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyname,      /* I - root name of keywords to read        */
            int  nstart,        /* I - starting index number                */
            int  nmax,          /* I - maximum number of keywords to return */
            long *value,        /* O - array of keyword values              */
            int  *nfound,       /* O - number of values that were returned  */
            int  *status)       /* IO - error status                        */
/*
  Read (get) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NMAX -1) inclusive.  
*/
{
    int nend, lenroot, ii, nkeys, mkeys, tstatus;
    long ival;
    char keyroot[FLEN_KEYWORD], keyindex[8], card[FLEN_CARD];
    char svalue[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = '\0';
    strncat(keyroot, keyname, 8);

    lenroot = strlen(keyroot);
    if (lenroot == 0 || lenroot > 7)     /* root must be 1 - 7 chars long */
        return(*status);

    for (ii=0; ii < lenroot; ii++)           /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);

    ffghps(fptr, &nkeys, &mkeys, status);  /*  get the number of keywords  */

    ffmaky(fptr, 3, status);  /* move to 3rd keyword (skip 1st 2 keywords) */

    for (ii=3; ii <= nkeys; ii++)  
    {
       if (ffgnky(fptr, card, status) > 0)     /*  get next keyword  */
           return(*status);

       if (strncmp(keyroot, card, lenroot) == 0)  /* see if keyword matches */
       {
           keyindex[0] = '\0';
           strncat(keyindex, &card[lenroot], 8-lenroot);  /*  copy suffix */

           tstatus = 0;
           if (ffc2ii(keyindex, &ival, &tstatus) <= 0)     /*  test suffix  */
               if (ival <= nend && ival >= nstart)
               {
                 ffpsvc(card, svalue, comm, status);   /*  parse the value */
                 ffc2i(svalue, &value[ival-nstart], status);  /* convert */
                 if (ival > *nfound)
                     *nfound = ival;           /* record the max index found */ 
              }
       }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgkne( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyname,      /* I - root name of keywords to read        */
            int  nstart,        /* I - starting index number                */
            int  nmax,          /* I - maximum number of keywords to return */
            float *value,       /* O - array of keyword values              */
            int  *nfound,       /* O - number of values that were returned  */
            int  *status)       /* IO - error status                        */
/*
  Read (get) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NMAX -1) inclusive.  
*/
{
    int nend, lenroot, ii, nkeys, mkeys, tstatus;
    long ival;
    char keyroot[FLEN_KEYWORD], keyindex[8], card[FLEN_CARD];
    char svalue[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = '\0';
    strncat(keyroot, keyname, 8);

    lenroot = strlen(keyroot);
    if (lenroot == 0 || lenroot > 7)     /*  root must be 1 - 7 chars long  */
        return(*status);

    for (ii=0; ii < lenroot; ii++)           /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);

    ffghps(fptr, &nkeys, &mkeys, status);  /*  get the number of keywords  */

    ffmaky(fptr, 3, status);  /* move to 3rd keyword (skip 1st 2 keywords) */

    for (ii=3; ii <= nkeys; ii++)  
    {
       if (ffgnky(fptr, card, status) > 0)     /*  get next keyword  */
           return(*status);

       if (strncmp(keyroot, card, lenroot) == 0)  /* see if keyword matches */
       {
           keyindex[0] = '\0';
           strncat(keyindex, &card[lenroot], 8-lenroot);  /*  copy suffix */

           tstatus = 0;
           if (ffc2ii(keyindex, &ival, &tstatus) <= 0)     /*  test suffix  */
               if (ival <= nend && ival >= nstart)
               {
                 ffpsvc(card, svalue, comm, status);   /*  parse the value */
                 ffc2r(svalue, &value[ival-nstart], status); /* convert */
                 if (ival > *nfound)
                     *nfound = ival;          /* record the max index found */
              }
       }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgknd( fitsfile *fptr,     /* I - FITS file pointer                    */
            char *keyname,      /* I - root name of keywords to read        */
            int  nstart,        /* I - starting index number                */
            int  nmax,          /* I - maximum number of keywords to return */
            double *value,      /* O - array of keyword values              */
            int  *nfound,       /* O - number of values that were returned  */
            int  *status)       /* IO - error status                        */
/*
  Read (get) an indexed array of keywords with index numbers between
  NSTART and (NSTART + NMAX -1) inclusive.  
*/
{
    int nend, lenroot, ii, nkeys, mkeys, tstatus;
    long ival;
    char keyroot[FLEN_KEYWORD], keyindex[8], card[FLEN_CARD];
    char svalue[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = '\0';
    strncat(keyroot, keyname, 8);

    lenroot = strlen(keyroot);
    if (lenroot == 0 || lenroot > 7)     /*  root must be 1 - 7 chars long  */
        return(*status);

    for (ii=0; ii < lenroot; ii++)           /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);

    ffghps(fptr, &nkeys, &mkeys, status);  /*  get the number of keywords  */

    ffmaky(fptr, 3, status);  /* move to 3rd keyword (skip 1st 2 keywords) */

    for (ii=3; ii <= nkeys; ii++)  
    {
       if (ffgnky(fptr, card, status) > 0)     /*  get next keyword  */
           return(*status);

       if (strncmp(keyroot, card, lenroot) == 0)   /* see if keyword matches */
       {
           keyindex[0] = '\0';
           strncat(keyindex, &card[lenroot], 8-lenroot);  /*  copy suffix */

           tstatus = 0;
           if (ffc2ii(keyindex, &ival, &tstatus) <= 0)      /*  test suffix */ 
               if (ival <= nend && ival >= nstart) /* is index within range? */
               {
                 ffpsvc(card, svalue, comm, status);   /*  parse the value */
                 ffc2d(svalue, &value[ival-nstart], status); /* convert */
                 if (ival > *nfound)
                     *nfound = ival;           /* record the max index found */
              }
       }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgtdm(fitsfile *fptr,  /* I - FITS file pointer                        */
           int colnum,      /* I - number of the column to read             */
           int maxdim,      /* I - maximum no. of dimensions to read;       */
           int *naxis,      /* O - number of axes in the data array         */
           long naxes[],    /* O - length of each data axis                 */
           int *status)     /* IO - error status                            */
/*
  parse the TDIMnnn keyword to get the dimensionality of a column
*/
{
    int tstatus;
    long dimsize;
    char keyname[FLEN_KEYWORD], tdimstr[FLEN_VALUE], comm[FLEN_COMMENT];
    char *loc, *lastloc;
    tcolumn *colptr;

    if (*status > 0)
        return(*status);

    if (colnum < 1 || colnum > fptr->tfield)
        return(*status = BAD_COL_NUM);

    colptr = fptr->tableptr;   /* set pointer to the first column */
    colptr += (colnum - 1);    /* increment to the correct column */

    ffkeyn("TDIM", colnum, keyname, status);      /* construct keyword name */
    tstatus = 0;
    ffgkys(fptr, keyname, tdimstr, comm, &tstatus); /* try reading keyword */

    if (tstatus)   /* TDIMnnn keyword doesn't exist? */
    {
        *naxis = 1;                   /* default = 1 dimensional */
        naxes[0] = colptr->trepeat;  /* default length = repeat count */
    }
    else
    {
        *naxis = 0;

        loc = strchr(tdimstr, '(' );  /* find the opening quote */
        if (!loc)
            return(*status = BAD_TDIM);

        while (loc)
        {
            loc++;
            dimsize = strtol(loc, &loc, 10);  /* read size of next dimension */
            if (*naxis < maxdim)
                naxes[*naxis] = dimsize;

            (*naxis)++;
            lastloc = loc;
            loc = strchr(loc, ',');  /* look for comma separating next dimension */
        }

        loc = strchr(lastloc, ')' );  /* check for the closing quote */
        if (!loc)
            return(*status = BAD_TDIM);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffghpr(fitsfile *fptr,  /* I - FITS file pointer                        */
           int maxdim,      /* I - maximum no. of dimensions to read;       */
           int *simple,     /* O - does file conform to FITS standard? 1/0  */
           int *bitpix,     /* O - number of bits per data value pixel      */
           int *naxis,      /* O - number of axes in the data array         */
           long naxes[],    /* O - length of each data axis                 */
           long *pcount,    /* O - number of group parameters (usually 0)   */
           long *gcount,    /* O - number of random groups (usually 1 or 0) */
           int *extend,     /* O - may FITS file haave extensions?          */
           int *status)     /* IO - error status                            */
/*
  Get keywords from the Header of the PRimary array:
  Check that the keywords conform to the FITS standard and return the
  parameters which determine the size and structure of the primary array
  or IMAGE extension.
*/
{
    int idummy;
    long ldummy;
    double ddummy;

    ffgphd(fptr, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend,
          &ddummy, &ddummy, &ldummy, &idummy, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffghtb(fitsfile *fptr,  /* I - FITS file pointer                        */
           int maxfield,    /* I - maximum no. of columns to read;          */
           long *naxis1,    /* O - length of table row in bytes             */
           long *naxis2,    /* O - number of rows in the table              */
           int *tfields,    /* O - number of columns in the table           */
           char **ttype,    /* O - name of each column                      */
           long *tbcol,     /* O - byte offset in row to each column        */
           char **tform,    /* O - value of TFORMn keyword for each column  */
           char **tunit,    /* O - value of TUNITn keyword for each column  */
           char *extnm,   /* O - value of EXTNAME keyword, if any         */
           int *status)     /* IO - error status                            */
/*
  Get keywords from the Header of the ASCII TaBle:
  Check that the keywords conform to the FITS standard and return the
  parameters which describe the table.
*/
{
    int ii, maxf, nfound, tstatus;
    long pcount, fields;
    char comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    if (ffgttb(fptr, naxis1, naxis2, &pcount, &fields, status) > 0)
        return(*status);

    *tfields = fields;

    if (maxfield < 0)
        maxf = *tfields;
    else
        maxf = minvalue(maxfield, *tfields);

    if (maxf > 0)
    {
        for (ii = 0; ii < maxf; ii++)
        {   /* initialize optional keyword values */
            *ttype[ii] = '\0';   
            *tunit[ii] = '\0';
        }

        ffgkns(fptr, "TTYPE", 1, maxf, ttype, &nfound, status);
        ffgkns(fptr, "TUNIT", 1, maxf, tunit, &nfound, status);

        if (*status > 0)
            return(*status);

        ffgknj(fptr, "TBCOL", 1, maxf, tbcol, &nfound, status);

        if (*status > 0 || nfound != maxf)
        {
        ffpmsg(
        "Required TBCOL keyword(s) not found in ASCII table header (ffghtb).");
        return(*status = NO_TBCOL);
        }

        ffgkns(fptr, "TFORM", 1, maxf, tform, &nfound, status);

        if (*status > 0 || nfound != maxf)
        {
        ffpmsg(
        "Required TFORM keyword(s) not found in ASCII table header (ffghtb).");
        return(*status = NO_TFORM);
        }

        extnm[0] = '\0';

        tstatus = *status;
        ffgkys(fptr, "EXTNAME", extnm, comm, status);

        if (*status == KEY_NO_EXIST)
            *status = tstatus;  /* keyword not required, so ignore error */

    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffghbn(fitsfile *fptr,  /* I - FITS file pointer                        */
           int maxfield,    /* I - maximum no. of columns to read;          */
           long *naxis2,    /* O - number of rows in the table              */
           int *tfields,    /* O - number of columns in the table           */
           char **ttype,    /* O - name of each column                      */
           char **tform,    /* O - TFORMn value for each column             */
           char **tunit,    /* O - TUNITn value for each column             */
           char *extnm,   /* O - value of EXTNAME keyword, if any         */
           long *pcount,    /* O - value of PCOUNT keyword                  */
           int *status)     /* IO - error status                            */
/*
  Get keywords from the Header of the ASCII TaBle:
  Check that the keywords conform to the FITS standard and return the
  parameters which describe the table.
*/
{
    int ii, maxf, nfound, tstatus;
    long naxis1, fields;
    char comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    if (ffgttb(fptr, &naxis1, naxis2, pcount, &fields, status) > 0)
        return(*status);

    *tfields = fields;

    if (maxfield < 0)
        maxf = *tfields;
    else
        maxf = minvalue(maxfield, *tfields);

    if (maxf > 0)
    {
        for (ii = 0; ii < maxf; ii++)
        {   /* initialize optional keyword values */
            *ttype[ii] = '\0';   
            *tunit[ii] = '\0';
        }

        ffgkns(fptr, "TTYPE", 1, maxf, ttype, &nfound, status);
        ffgkns(fptr, "TUNIT", 1, maxf, tunit, &nfound, status);

        if (*status > 0)
            return(*status);

        ffgkns(fptr, "TFORM", 1, maxf, tform, &nfound, status);

        if (*status > 0 || nfound != maxf)
        {
        ffpmsg(
        "Required TFORM keyword(s) not found in binary table header (ffghbn).");
        return(*status = NO_TFORM);
        }

        extnm[0] = '\0';

        tstatus = *status;
        ffgkys(fptr, "EXTNAME", extnm, comm, status);

        if (*status == KEY_NO_EXIST)
            *status = tstatus;  /* keyword not required, so ignore error */

    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgphd(fitsfile *fptr,  /* I - FITS file pointer                        */
           int maxdim,      /* I - maximum no. of dimensions to read;       */
           int *simple,     /* O - does file conform to FITS standard? 1/0  */
           int *bitpix,     /* O - number of bits per data value pixel      */
           int *naxis,      /* O - number of axes in the data array         */
           long naxes[],    /* O - length of each data axis                 */
           long *pcount,    /* O - number of group parameters (usually 0)   */
           long *gcount,    /* O - number of random groups (usually 1 or 0) */
           int *extend,     /* O - may FITS file haave extensions?          */
           double *bscale,  /* O - array pixel linear scaling factor        */
           double *bzero,   /* O - array pixel linear scaling zero point    */
           long *blank,     /* O - value used to represent undefined pixels */
           int *nspace,     /* O - number of blank keywords prior to END    */
           int *status)     /* IO - error status                            */
{
/*
  Get the Primary HeaDer parameters.  Check that the keywords conform to
  the FITS standard and return the parameters which determine the size and
  structure of the primary array or IMAGE extension.
*/
    int unknown, found_end, tstatus, ii, nextkey;
    long longbitpix, longnaxis, axislen;
    char message[FLEN_ERRMSG], keyword[FLEN_KEYWORD];
    char name[FLEN_KEYWORD], value[FLEN_VALUE], comm[FLEN_COMMENT];

    if (*status > 0)
        return(*status);

    *simple = 1;
    unknown = 0;

    /*--------------------------------------------------------------------*/
    /*  Get 1st keyword of HDU and test whether it is SIMPLE or XTENSION  */
    /*--------------------------------------------------------------------*/
    ffgkyn(fptr, 1, name, value, comm, status);

    if (fptr->curhdu == 0)  /* Is this the beginning of the FITS file? */
    {
        if (!strcmp(name, "SIMPLE"))
        {
            if (value[0] == 'F')
                *simple=0;          /* not a simple FITS file */
            else if (value[0] != 'T')
                return(*status = BAD_SIMPLE);
        }

        else
        {
            sprintf(message,
                   "First keyword of the file is not SIMPLE: %s", name);
            ffpmsg(message);
            return(*status = NO_SIMPLE);
        }
    }

    else    /* not beginning of the file, so presumably an IMAGE extension */
    {
        if (!strcmp(name, "XTENSION"))
        {
            if (strcmp(value, "\'IMAGE   \'" ) &&
                strcmp(value, "\'IUEIMAGE\'") )
            {
                unknown = 1;  /* unknown type of extension; press on anyway */
                sprintf(message,
                   "This is not an IMAGE extension: %s", value);
                ffpmsg(message);
            }
        }

        else  /* error: 1st keyword of extension != XTENSION */
        {
            sprintf(message,
            "First keyword of the extension is not XTENSION: %s", name);
            ffpmsg(message);
            return(*status = NO_XTENSION);
        }
    }

    /*----------------------------------------------------------------*/
    /*  Get 2nd keyword;  test whether it is BITPIX with legal value  */
    /*----------------------------------------------------------------*/
    ffgkyn(fptr, 2, name, value, comm, status);  /* BITPIX = 2nd keyword */

    if (strcmp(name, "BITPIX"))
    {
        sprintf(message,
        "Second keyword of the extension is not BITPIX: %s", name);
        ffpmsg(message);
        return(*status = NO_BITPIX);
    }

    if (ffc2ii(value,  &longbitpix, status) > 0)
    {
        sprintf(message,
        "Value of BITPIX keyword is not an integer: %s", value);
        ffpmsg(message);
        return(*status = BAD_BITPIX);
    }
    else if (longbitpix != BYTE_IMG && longbitpix != SHORT_IMG &&
             longbitpix != LONG_IMG &&
             longbitpix != FLOAT_IMG && longbitpix != DOUBLE_IMG)
    {
        sprintf(message,
        "Illegal value for BITPIX keyword: %s", value);
        ffpmsg(message);
        return(*status = BAD_BITPIX);
    }
    *bitpix = longbitpix;  /* do explicit type conversion */


    /*---------------------------------------------------------------*/
    /*  Get 3rd keyword;  test whether it is NAXIS with legal value  */
    /*---------------------------------------------------------------*/
    ffgtkn(fptr, 3, "NAXIS",  &longnaxis, status);

    if (*status == BAD_ORDER)
        return(*status = NO_NAXIS);
    else if (*status == NOT_POS_INT)
        return(*status = BAD_NAXIS);
    else
        *naxis = longnaxis;  /* do explicit type conversion */


    /*---------------------------------------------------------*/
    /*  Get the next NAXISn keywords and test for legal values */
    /*---------------------------------------------------------*/
    for (ii=0, nextkey=4; ii < *naxis; ii++, nextkey++)
    {
        ffkeyn("NAXIS", ii+1, keyword, status);
        ffgtkn(fptr, 4+ii, keyword, &axislen, status);

        if (*status == BAD_ORDER)
            return(*status = NO_NAXES);
        else if (*status == NOT_POS_INT)
            return(*status = BAD_NAXES);
        else if (ii < maxdim)
            naxes[ii] = axislen;
    }


    /*---------------------------------------------------------*/
    /*  now look for other keywords of interest:               */
    /*  BSCALE, BZERO, BLANK, PCOUNT, GCOUNT, EXTEND, and END  */
    /*---------------------------------------------------------*/

    /*  initialize default values in case keyword is not present */
    *bscale = 1.0;
    *bzero  = 0.0;
    *pcount = 0;
    *gcount = 1;
    *extend = 0;
    *blank = NO_NULL; /* by default, no null value for BITPIX = 8, 16, 32 */

    *nspace = 0;
    found_end = 0;
    tstatus = *status;

    for (; !found_end; nextkey++)  
    {
      /* get next keyword */
      if (ffgkyn(fptr, nextkey, name, value, comm, status) > 0)
      {
        if (*status == KEY_OUT_BOUNDS)
        {
          found_end = 1;  /* simply hit the end of the header */
          *status = tstatus;  /* reset error status */
        }
        else          
        {
          ffpmsg("Failed to find the END keyword in header (ffgphd).");
        }
      }

      else   /* got the next keyword without error */
      {
        if (!strcmp(name, "BSCALE"))
        {
            if (ffc2dd(value, bscale, status) > 0) /* convert to double */
            {
                sprintf(message,
                "Error reading BSCALE keyword value as a double: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "BZERO"))
        {
            if (ffc2dd(value, bzero, status) > 0) /* convert to double */
            {
                sprintf(message,
                "Error reading BZERO keyword value as a double: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "BLANK"))
        {
            if (ffc2ii(value, blank, status) > 0) /* convert to long */
            {
                sprintf(message,
                "Error reading BLANK keyword value as an integer: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "PCOUNT"))
        {
            if (ffc2ii(value, pcount, status) > 0) /* convert to long */
            {
                sprintf(message,
                "Error reading PCOUNT keyword value as an integer: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "GCOUNT"))
        {
            if (ffc2ii(value, gcount, status) > 0) /* convert to long */
            {
                sprintf(message,
                "Error reading GCOUNT keyword value as an integer: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "EXTEND"))
        {
            if (ffc2ll(value, extend, status) > 0) /* convert to int */
            {
                sprintf(message,
                "Error reading EXTEND keyword value as a logical: %s", value);
                ffpmsg(message);
            }
        }

        else if (!strcmp(name, "END"))
            found_end = 1;

        else if (!name[0] && !value[0] && !comm[0])
            *nspace = *nspace + 1;  /* this is a blank card in the header */

        else
            *nspace = 0;  /* reset count of blank keywords immediately
                            before the END keyword to zero   */
      }

      if (*status > 0)  /* exit on error after writing error message */
      {
        if (fptr->curhdu == 0)
            ffpmsg(
            "Failed to read the required primary array header keywords.");
        else
            ffpmsg(
            "Failed to read the required image extension header keywords.");

        return(*status);
      }
    }


    if (unknown)
       *status = NOT_IMAGE;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgttb(fitsfile *fptr,      /* I - FITS file pointer*/
           long *rowlen,        /* O - length of a table row, in bytes */
           long *nrows,         /* O - number of rows in the table */
           long *pcount,        /* O - value of PCOUNT keyword */
           long *tfields,       /* O - number of fields in the table */
           int *status)         /* IO - error status    */
{
/*
  Get and Test TaBle;
  Test that this is a legal ASCII or binary table and get some keyword values.
  We assume that the calling routine has already tested the 1st keyword
  of the extension to ensure that this is really a table extension.
*/
    if (*status > 0)
        return(*status);

    if (fftkyn(fptr, 2, "BITPIX", "8", status) == BAD_ORDER) /* 2nd keyword */
        return(*status = NO_BITPIX);  /* keyword not BITPIX */
    else if (*status == NOT_POS_INT)
        return(*status = BAD_BITPIX); /* value != 8 */

    if (fftkyn(fptr, 3, "NAXIS", "2", status) == BAD_ORDER) /* 3rd keyword */
        return(*status = NO_NAXIS);  /* keyword not NAXIS */
    else if (*status == NOT_POS_INT)
        return(*status = BAD_NAXIS); /* value != 2 */

    if (ffgtkn(fptr, 4, "NAXIS1", rowlen, status) == BAD_ORDER) /* 4th key */
        return(*status = NO_NAXES);  /* keyword not NAXIS1 */
    else if (*status == NOT_POS_INT)
        return(*status == BAD_NAXES); /* bad NAXIS1 value */

    if (ffgtkn(fptr, 5, "NAXIS2", nrows, status) == BAD_ORDER) /* 5th key */
        return(*status = NO_NAXES);  /* keyword not NAXIS2 */
    else if (*status == NOT_POS_INT)
        return(*status == BAD_NAXES); /* bad NAXIS2 value */

    if (ffgtkn(fptr, 6, "PCOUNT", pcount, status) == BAD_ORDER) /* 6th key */
        return(*status = NO_PCOUNT);  /* keyword not PCOUNT */
    else if (*status == NOT_POS_INT)
        return(*status = BAD_PCOUNT); /* bad PCOUNT value */

    if (fftkyn(fptr, 7, "GCOUNT", "1", status) == BAD_ORDER) /* 7th keyword */
        return(*status = NO_GCOUNT);  /* keyword not GCOUNT */
    else if (*status == NOT_POS_INT)
        return(*status = BAD_GCOUNT); /* value != 1 */

    if (ffgtkn(fptr, 8, "TFIELDS", tfields, status) == BAD_ORDER) /* 8th key*/
        return(*status = NO_TFIELDS);  /* keyword not TFIELDS */
    else if (*status == NOT_POS_INT || *tfields > 999)
        return(*status == BAD_TFIELDS); /* bad TFIELDS value */


    if (*status > 0)
       ffpmsg(
       "Error reading required keywords in the table header (FTGTTB).");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgtkn(fitsfile *fptr,  /* I - FITS file pointer              */
           int numkey,      /* I - number of the keyword to read  */
           char *name,      /* I - expected name of the keyword   */
           long *value,     /* O - integer value of the keyword   */
           int *status)     /* IO - error status                  */
{
/*
  test that keyword number NUMKEY has the expected name and get the
  integer value of the keyword.  Return an error if the keyword
  name does not match the input name, or if the value of the
  keyword is not a positive integer.
*/
    char keyname[FLEN_KEYWORD], valuestring[FLEN_VALUE];
    char comm[FLEN_COMMENT], message[FLEN_ERRMSG];
   
    if (*status > 0)
        return(*status);
    
    keyname[0] = '\0';
    valuestring[0] = '\0';

    if (ffgkyn(fptr, numkey, keyname, valuestring, comm, status) <= 0)
    {
        if (strcmp(keyname, name) )
            *status = BAD_ORDER;  /* incorrect keyword name */

        ffc2ii(valuestring, value, status);  /* convert to integer */

        if (*status > 0 || *value < 0 )
           *status = NOT_POS_INT;
    }

    if (*status > 0)
    {
        sprintf(message,
          "ffgtkn found unexpected keyword or value for keyword no. %d.",
          numkey);
        ffpmsg(message);

        sprintf(message,
          " Expected positive integer keyword %s, but instead", name);
        ffpmsg(message);

        sprintf(message,
          " found keyword %s with value %s", keyname, valuestring);
        ffpmsg(message);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fftkyn(fitsfile *fptr,  /* I - FITS file pointer              */
           int numkey,      /* I - number of the keyword to read  */
           char *name,      /* I - expected name of the keyword   */
           char *value,     /* I - expected value of the keyword  */
           int *status)     /* IO - error status                  */
{
/*
  test that keyword number NUMKEY has the expected name and the
  expected value string.
*/
    char keyname[FLEN_KEYWORD], valuestring[FLEN_VALUE];
    char comm[FLEN_COMMENT], message[FLEN_ERRMSG];
   
    if (*status > 0)
        return(*status);
    
    keyname[0] = '\0';
    valuestring[0] = '\0';

    if (ffgkyn(fptr, numkey, keyname, valuestring, comm, status) <= 0)
    {
        if (strcmp(keyname, name) )
            *status = BAD_ORDER;  /* incorrect keyword name */

        if (strcmp(value, valuestring) )
            *status = NOT_POS_INT;  /* incorrect keyword value */
    }

    if (*status > 0)
    {
        sprintf(message,
          "fftkyn found unexpected keyword or value for keyword no. %d.",
          numkey);
        ffpmsg(message);

        sprintf(message,
          " Expected keyword %s with value %s, but", name, value);
        ffpmsg(message);

        sprintf(message,
          " found keyword %s with value %s", keyname, valuestring);
        ffpmsg(message);
    }

    return(*status);
}

