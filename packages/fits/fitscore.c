/*  This file, fitscore.c, contains the core set of FITSIO routines.       */

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
#include <math.h>
#include <ctype.h>
#include "fitsio2.h"
/*--------------------------------------------------------------------------*/
float ffvers(float *version)  /* IO - version number */
/*
  return the current version number of the FITSIO software
*/
{
 *version = 1.23;  /* 24 Apr 1997 */

 /*   *version = 1.22;   18 Apr 1997 */
 /*   *version = 1.21;   26 Mar 1997 */
 /*   *version = 1.2;    29 Jan 1997 */
 /*   *version = 1.11;   04 Dec 1996 */
 /*   *version = 1.101;  13 Nov 1996 */
 /*   *version = 1.1;     6 Nov 1996 */
 /*   *version = 1.04;   17 Sep 1996 */
 /*   *version = 1.03;   20 Aug 1996 */
 /*   *version = 1.02;   15 Aug 1996 */
 /*   *version = 1.01;   12 Aug 1996 */

    return(*version);
}
/*--------------------------------------------------------------------------*/
void ffgerr(int status,     /* I - error status value */
            char *errtext)  /* O - error message (max 30 char long + null) */
/*
  Return a short descriptive error message that corresponds to the input
  error status value.  The message may be up to 30 characters long, plus
  the terminating null character.
*/
{
  if (status < 300)
  {
    switch (status) {

    case 0:
       strcpy(errtext, "OK - no error");
       break;
    case 101:
       strcpy(errtext, "same input and output files");
       break;
    case 103:
       strcpy(errtext, "attempt to open too many files");
       break;
    case 104:
       strcpy(errtext, "could not open the named file");
       break;
    case 105:
       strcpy(errtext, "couldn't create the named file");
       break;
    case 106:
       strcpy(errtext, "error writing to FITS file");
       break;
    case 107:
       strcpy(errtext, "tried to move past end of file");
       break;
    case 108:
       strcpy(errtext, "error reading from FITS file");
       break;
    case 110:
       strcpy(errtext, "could not close the file");
       break;
    case 111:
       strcpy(errtext, "array dimensions too big");
       break;
    case 112:
       strcpy(errtext, "cannot write to readonly file");
       break;
    case 201:
       strcpy(errtext, "header already has keywords");
       break;
    case 202:
       strcpy(errtext, "keyword not found in header");
       break;
    case 203:
       strcpy(errtext, "keyword number out of bounds");
       break;
    case 204:
       strcpy(errtext, "keyword value field is blank");
       break;
    case 205:
       strcpy(errtext, "string missing closing quote");
       break;
    case 207:
       strcpy(errtext, "illegal character in keyword");
       break;
    case 208:
       strcpy(errtext, "required keywords out of order");
       break;
    case 209:
       strcpy(errtext, "keyword value not positive int");
       break;
    case 210:
       strcpy(errtext, "END keyword not found");
       break;
    case 211:
       strcpy(errtext, "illegal BITPIX keyword value");
       break;
    case 212:
       strcpy(errtext, "illegal NAXIS keyword value");
       break;
    case 213:
       strcpy(errtext, "illegal NAXISn keyword value");
       break;
    case 214:
       strcpy(errtext, "illegal PCOUNT keyword value");
       break;
    case 215:
       strcpy(errtext, "illegal GCOUNT keyword value");
       break;
    case 216:
       strcpy(errtext, "illegal TFIELDS keyword value");
       break;
    case 217:
       strcpy(errtext, "negative table row size");
       break;
    case 218:
       strcpy(errtext, "negative number of rows");
       break;
    case 219:
       strcpy(errtext, "named column not found");
       break;
    case 220:
       strcpy(errtext, "illegal SIMPLE keyword value");
       break;
    case 221:
       strcpy(errtext, "first keyword not SIMPLE");
       break;
    case 222:
       strcpy(errtext, "second keyword not BITPIX");
       break;
    case 223:
       strcpy(errtext, "third keyword not NAXIS");
       break;
    case 224:
       strcpy(errtext, "missing NAXISn keywords");
       break;
    case 225:
       strcpy(errtext, "first keyword not XTENSION");
       break;
    case 226:
       strcpy(errtext, "CHDU not an ASCII table");
       break;
    case 227:
       strcpy(errtext, "CHDU not a binary table");
       break;
    case 228:
       strcpy(errtext, "PCOUNT keyword not found");
       break;
    case 229:
       strcpy(errtext, "GCOUNT keyword not found");
       break;
    case 230:
       strcpy(errtext, "TFIELDS keyword not found");
       break;
    case 231:
       strcpy(errtext, "missing TBCOLn keyword");
       break;
    case 232:
       strcpy(errtext, "missing TFORMn keyword");
       break;
    case 233:
       strcpy(errtext, "CHDU not an IMAGE extension");
       break;
    case 234:
       strcpy(errtext, "illegal TBCOLn keyword value");
       break;
    case 235:
       strcpy(errtext, "CHDU not a table extension");
       break;
    case 236:
       strcpy(errtext, "column exceeds width of table");
       break;
    case 237:
       strcpy(errtext, "more than 1 matching col. name");
       break;
    case 241:
       strcpy(errtext, "row width not = field widths");
       break;
    case 251:
       strcpy(errtext, "unknown FITS extension type");
       break;
    case 252:
       strcpy(errtext, "unknown FITS record type");
       break;
    case 253:
       strcpy(errtext, "END keyword is not blank");
       break;
    case 261:
       strcpy(errtext, "illegal TFORM format code");
       break;
    case 262:
       strcpy(errtext, "unknown TFORM datatype code");
       break;
    case 263:
       strcpy(errtext, "illegal TDIMn keyword value");
       break;
    }
  }
  else
  {
    switch(status) {

    case 301:
       strcpy(errtext, "illegal HDU number");
       break;
    case 302:
       strcpy(errtext, "column number < 1 or > tfields");
       break;
    case 304:
       strcpy(errtext, "negative byte address");
       break;
    case 306:
       strcpy(errtext, "negative number of elements");
       break;
    case 307:
       strcpy(errtext, "bad first row number");
       break;
    case 308:
       strcpy(errtext, "bad first element number");
       break;
    case 309:
       strcpy(errtext, "not an ASCII (A) column");
       break;
    case 310:
       strcpy(errtext, "not a logical (L) column");
       break;
    case 311:
       strcpy(errtext, "bad ASCII table datatype");
       break;
    case 312:
       strcpy(errtext, "bad binary table datatype");
       break;
    case 314:
       strcpy(errtext, "null value not defined");
       break;
    case 317:
       strcpy(errtext, "not a variable length column");
       break;
    case 320:
       strcpy(errtext, "illegal number of dimensions");
       break;
    case 322:
       strcpy(errtext, "BSCALE or TSCALn = 0.");
       break;
    case 323:
       strcpy(errtext, "illegal axis length < 1");
       break;
    case 401:
       strcpy(errtext, "bad int to string conversion");
       break;
    case 402:
       strcpy(errtext, "bad float to string conversion");
       break;
    case 403:
       strcpy(errtext, "keyword value not integer");
       break;
    case 404:
       strcpy(errtext, "keyword value not logical");
       break;
    case 405:
       strcpy(errtext, "keyword value not floating pt");
       break;
    case 406:
       strcpy(errtext, "keyword value not double");
       break;
    case 407:
       strcpy(errtext, "bad string to int conversion");
       break;
    case 408:
       strcpy(errtext, "bad string to float conversion");
       break;
    case 409:
       strcpy(errtext, "bad string to double convert");
       break;
    case 410:
       strcpy(errtext, "illegal datatype code value");
       break;
    case 411:
       strcpy(errtext, "illegal no. of decimals");
       break;
    case 412:
       strcpy(errtext, "datatype conversion overflow");
       break;
    case 501:
       strcpy(errtext, "WCS angle too large");
       break;
    case 502:
       strcpy(errtext, "bad WCS coordinate");
       break;
    case 503:
       strcpy(errtext, "error in WCS calculation");
       break;
    case 504:
       strcpy(errtext, "bad WCS projection type");
       break;
    case 505:
       strcpy(errtext, "WCS keywords not found");
       break;
    }
  }
  return;
}
/*--------------------------------------------------------------------------*/
void ffpmsg(const char *err_message)
/*
  put message on to error stack
*/
{
    ffxmsg(1, (char *)err_message);
    return;
}
/*--------------------------------------------------------------------------*/
int ffgmsg(char *err_message)
/*
  get oldest message from error stack
*/
{
    ffxmsg(-1, err_message);
    return(*err_message);
}
/*--------------------------------------------------------------------------*/
void ffcmsg(void)
/*
  erase all messages in the error stack
*/
{
    char *dummy;

    ffxmsg(0, dummy);
    return;
}
/*--------------------------------------------------------------------------*/
void ffxmsg( int action,
            char *errmsg)
/*
  general routine to get, put, or clear the error message stack 
*/
{
    int ii;
    size_t len;
#define errmsgsiz 50
    static char *txtbuff[errmsgsiz];
    static nummsg = 0;

    if (action == -1)  /* return and remove oldest message from stack */ 
    {
      if (nummsg > 0)
      {
        strcpy(errmsg, txtbuff[0]);   /* copy oldest message to output */

        free(txtbuff[0]);  /* free the memory for this msg */
           
        nummsg--;  
        for (ii = 0; ii < nummsg; ii++)
             txtbuff[ii] = txtbuff[ii + 1]; /* shift remaining pointers */
      }
      else
          errmsg[0] = 0;  /*  no messages in the stack */
    }
    else if (action == 1)  /* add new message to stack */
    {
      if (nummsg == errmsgsiz)
      {
        free(txtbuff[0]);  /* buffer full; delete oldest msg */
        nummsg--;
        for (ii = 0; ii < nummsg; ii++)
             txtbuff[ii] = txtbuff[ii + 1];   /* shift remaining pointers */
      }

      len = minvalue(strlen(errmsg), 80);

      txtbuff[nummsg] = (char *) malloc(len + 1);

      if (!txtbuff[nummsg])
      {
        printf("\nmalloc failed in the ffpmsg routine of cfitsio.\n");
        printf("%s\n", errmsg);
      }
      else
      {
        *txtbuff[nummsg] = '\0';  /* initialize a null string */
        strncat(txtbuff[nummsg], errmsg, len);
        nummsg++;
      }
    }
    else if (action == 0)
    {
      for (ii = 0; ii < nummsg; ii ++)
        free(txtbuff[ii]);

      nummsg = 0;
    }
    return;
}
/*--------------------------------------------------------------------------*/
int fftkey(char *keyword,    /* I -  keyword name */
           int *status)      /* IO - error status */
/*
  Test that the keyword name conforms to the FITS standard.  Must contain
  only capital letters, digits, minus or underscore chars.  Trailing spaces
  are allowed.
*/
{
    size_t maxchr, ii;
    int spaces=0;
    char msg[81];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    maxchr=strlen(keyword);
    if (maxchr > 8)
        maxchr = 8;

    for (ii = 0; ii < maxchr; ii++)
    {
        if ( (keyword[ii] >= 'A' && keyword[ii] <= 'Z') ||
             (keyword[ii] >= '0' && keyword[ii] <= '9') ||
              keyword[ii] == '-' || keyword[ii] == '_'   )
              {
                if (spaces)
                {
                 sprintf(msg, "Keyword name contains embedded space(s): %.8s",
                        keyword);
                    ffpmsg(msg);
                    return(*status = BAD_KEYCHAR);        
                }
              }
        else if (keyword[ii] == ' ')
            spaces = 1;

        else     
        {
            sprintf(msg, "Character %d in this keyword is illegal: %.8s",
                    ii+1, keyword);
            ffpmsg(msg);

            /* explicitly flag the 2 most common cases */
            if (keyword[ii] == 0) 
                ffpmsg(" (This a NULL (0) character).");                
            else if (keyword[ii] == 9)
                ffpmsg(" (This an ASCII TAB (9) character).");                

            return(*status = BAD_KEYCHAR);        

        }                
    }
    return(*status);        
}
/*--------------------------------------------------------------------------*/
int fftrec(char *card,       /* I -  keyword card to test */
           int *status)      /* IO - error status */
/*
  Test that the keyword card conforms to the FITS standard.  Must contain
  only printable ASCII characters;
*/
{
    size_t ii, maxchr;
    char msg[81];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    maxchr = strlen(card);

    for (ii = 8; ii < maxchr; ii++)
    {
        if (!isprint(card[ii]))
        {
            sprintf(msg, "Character %d in this keyword record is illegal:",
              (int) (ii+1) );
            ffpmsg(msg);

            strncpy(msg, card, 80);
            msg[80] = '\0';
            ffpmsg(msg);
            return(*status = BAD_KEYCHAR);        
        }
    }
    return(*status);        
}
/*--------------------------------------------------------------------------*/
void ffupch(char *string)
/*
  convert string to upper case, in place.
*/
{
    size_t len, ii;

    len = strlen(string);
    for (ii = 0; ii < len; ii++)
        string[ii] = toupper(string[ii]);
    return;
}
/*--------------------------------------------------------------------------*/
void ffmkky(char *keyname,   /* I - keyword name    */
            char *value,     /* I - keyword value   */
            char *comm,      /* I - keyword comment */
            char *card)      /* O - constructed keyword card */
/*
  Make a complete FITS 80-byte keyword card from the input name, value and
  comment strings. Output card is null terminated without any trailing blanks.
*/
{
    size_t len, ii;

    strncpy(card, keyname, 8);   /* copy keyword name to buffer */
    card[8] = '\0'; 
   
    len = strlen(card);
    for (ii = len; ii < 8; ii++)
        card[ii] = ' ';          /* pad keyword name with spaces */

    card[8] = '=';              /* append '= ' in columns 9-10 */
    card[9] = ' ';
    card[10] = '\0';                 /* terminate the partial string */

    len = strlen(value);        
    if (len > 0)
    {
        if (value[0] == '\'')  /* is this a quoted string value? */
        {
            strncat(card, value, 70);       /* append the value string */
            len = strlen(card);

            if (comm[0] != 0)
            {
                for (ii = len; ii < 30; ii++)
                  strcat(card, " "); /* add spaces so field ends in col 30 */
            }
        }
        else
        {
            for (ii = len; ii < 20; ii++)
                strcat(card, " ");  /* add spaces so field ends in col 30 */

            strncat(card, value, 70);       /* append the value string */
        }

        len = strlen(card);
        if ((len < 77) && ( strlen(comm) > 0) )  /* room for a comment? */
        {
            strcat(card, " / ");   /* append comment separator */
            strncat(card, comm, 77 - len); /* append comment (what fits) */
        } 
    }
    else
    {
        card[8] = ' ';   /* keywords with no value have no equal sign */ 
        strncat(card, comm, 70);   /* append comment (whatever fits) */
    }
}
/*--------------------------------------------------------------------------*/
int ffmkey(fitsfile *fptr,    /* I - FITS file pointer  */
           char *card,        /* I - card string value  */
           int *status)       /* IO - error status      */
/*
  replace the previously read card (i.e. starting 80 bytes before the
  fptr->nextkey position) with the contents of the input card.
*/
{
    char tcard[81];
    size_t len, ii;

    strncpy(tcard,card,80);
    tcard[80] = '\0';

    len = strlen(tcard);
    for (ii=len; ii < 80; ii++)    /* fill card with spaces if necessary */
        tcard[ii] = ' ';

    for (ii=0; ii < 8; ii++)       /* make sure keyword name is uppercase */
        tcard[ii] = toupper(tcard[ii]);

    fftkey(tcard, status);        /* test keyword name contains legal chars */
    fftrec(tcard, status);        /* test rest of keyword for legal chars   */

    /* move position of keyword to be over written */
    ffmbyt(fptr, (fptr->nextkey) - 80, REPORT_EOF, status); 
    ffpbyt(fptr, 80, tcard, status);   /* write the 80 byte card */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffkeyn(char *keyroot,   /* I - root string for keyword name */
           int value,       /* I - index number to be appended to root name */
           char *keyname,   /* O - output root + index keyword name */
           int *status)     /* IO - error status  */
/*
  Construct a keyword name string by appending the index number to the root.
  e.g., if root = "TTYPE" and value = 12 then keyname = "TTYPE12".
*/
{
    char suffix[4];
    size_t rootlen;

    keyname[0] = '\0';            /* initialize output name to null */
    rootlen = strlen(keyroot);

    if (rootlen == 0 || rootlen > 7 || value < 1 || value > 999)
       return(*status = 206);

    sprintf(suffix, "%d", value); /* construct keyword suffix */

    if ( strlen(suffix) + rootlen > 8)
       return(*status = 206);

    strcpy(keyname, keyroot);   /* copy root string to name string */
    strcat(keyname, suffix);    /* append suffix to the root */
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffnkey(int value,       /* I - index number to be appended to root name */
           char *keyroot,   /* I - root string for keyword name */
           char *keyname,   /* O - output root + index keyword name */
           int *status)     /* IO - error status  */
/*
  Construct a keyword name string by appending the root string to the index
  number. e.g., if root = "TTYPE" and value = 12 then keyname = "12TTYPE".
*/
{
    size_t rootlen;

    keyname[0] = '\0';            /* initialize output name to null */
    rootlen = strlen(keyroot);

    if (rootlen == 0 || rootlen > 7 || value < 1 || value > 999)
       return(*status = 206);

    sprintf(keyname, "%d", value); /* construct keyword prefix */

    if (rootlen +  strlen(keyname) > 8)
       return(*status = 206);

    strcat(keyname, keyroot);  /* append root to the prefix */
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpsvc(char *card,    /* I - FITS header card (nominally 80 bytes long) */
           char *value,   /* O - value string parsed from the card */
           char *comm,    /* O - comment string parsed from the card */
           int *status)   /* IO - error status   */
/*
  ParSe the Value and Comment strings from the input header card string.
  If the card contains a quoted string value, the returned value string
  includes the enclosing quote characters.  If comm = NULL, don't return
  the comment string.
*/
{
    int jj;
    size_t ii, cardlen, nblank;
    char errmsg[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    value[0] = '\0';
    if (comm)
        comm[0] = '\0';

    cardlen = strlen(card);
    
    if (cardlen < 9  ||
        strncmp(card, "COMMENT ", 8) == 0 ||  /* keywords with no value */
        strncmp(card, "HISTORY ", 8) == 0 ||
        strncmp(card, "END     ", 8) == 0 ||
        strncmp(card, "        ", 8) == 0 ||
        strncmp(&card[8],      "= ", 2) != 0  ) /* no '= ' in cols 9-10 */
    {
        /*  no value and the comment extends from cols 9 - 80  */
        if (comm != NULL)
        {
          if (cardlen > 8)
             strcpy(comm, &card[8]);

          jj=strlen(comm);
          for (jj--; jj >= 0; jj--)  /* replace trailing blanks with nulls */
          {
            if (comm[jj] == ' ')
                comm[jj] = '\0';
            else
                break;
          }
        }
        return(*status);
    }

    nblank = strspn(&card[10], " ");  /*  find number of leading blanks */

    if (nblank + 10 == cardlen)
    {
        strcpy(errmsg,"The keyword ");
        strncat(errmsg, card, 8);
        strcat(errmsg, " has no value string after the equal sign:");
        ffpmsg(errmsg);
        ffpmsg(card);
        return(*status = NO_VALUE);
    }

    ii = nblank + 10;
    if (card[ii] == '\'' )  /* is this a quoted string value? */
    {
        value[0] = card[ii];
        for (jj=1, ii++; ii < cardlen; ii++, jj++)
        {
            if (card[ii] == '\'')  /*  is this the closing quote?  */
            {
                if (card[ii+1] == '\'')  /* 2 successive quotes? */ 
                {
                   value[jj] = card[ii];
                   ii++;  
                   jj++;
                }
                else
                {
                    value[jj] = card[ii];
                    break;   /* found the closing quote, so exit this loop  */
                }
            }
            value[jj] = card[ii];  /* copy the next character to the output */
        }

        if (ii == cardlen)
        {
            value[jj] = '\0';  /*  terminate the bad value string  */
            ffpmsg("This keyword string value has no closing quote:");
            ffpmsg(card);
            return(*status = 205);
        }
        else
        {
            value[jj+1] = '\0';  /*  terminate the good value string  */
            ii++;   /*  point to the character following the value  */
        }
    }

    else   /*  an integer, floating point, or logical FITS value string  */
    {
        nblank = strcspn(&card[ii], " /");  /* find the end of the token */
        strncpy(value, &card[ii], nblank);
        value[nblank] = '\0';
        ii = ii + nblank;
    }

    /*  now find the comment string, if any  */
    if (comm)
    {
      comm[0] = '\0';

      nblank = strspn(&card[ii], " ");  /*  find next non-space character  */
      ii = ii + nblank;

      if (ii < 80)
      {
        if (card[ii] == '/')   /*  ignore the slash separator  */
        {
            ii++;
            if (card[ii] == ' ')  /*  also ignore the following space  */
                ii++;
        }
        strcat(comm, &card[ii]);  /*  copy the remaining characters  */

        jj=strlen(comm);
        for (jj--; jj >= 0; jj--)  /* replace trailing blanks with nulls */
        {
            if (comm[jj] == ' ')
                comm[jj] = '\0';
            else
                break;
        }
      }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgthd(char *tmplt, /* I - input header template string */
           char *card,  /* O - returned FITS header record */
           int *hdtype, /* O - how to interpreter the returned card string */ 
            /*
              -2 = modify the name of a keyword; the old keyword name
                   is returned starting at address chars[0]; the new name
                   is returned starting at address char[40] (to be consistent
                   with the Fortran version).  Both names are null terminated. 
              -1 = card contains the name of a keyword that is to be deleted
               0 = append this keyword if it doesn't already exist, or 
                   modify the value if the keyword already exists.
               1 = append this comment keyword ('HISTORY', 
                   'COMMENT', or blank keyword name) 
               2  =  this is the END keyword; do not write it to the header
            */
           int *status)   /* IO - error status   */
/*
  'Get Template HeaDer'
  parse a template header line and create a formated
  character string which is suitable for appending to a FITS header 
*/
{
    char keyname[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
    char *tok, *suffix;
    int len, vlen, more;
    double dval;

    if (*status > 0)
        return(*status);

    card[0]   = '\0';
    *hdtype   = 0;

    if (!strncmp(tmplt, "        ", 8) )
    {
        /* if first 8 chars of template are blank, then this is a comment */
        strncat(card, tmplt, 80);
        *hdtype = 1;
        return(*status);
    }

    tok = tmplt;   /* point to start of template string */
 
    keyname[0] = '\0';
    value[0]   = '\0';
    comment[0] = '\0';

    len = strspn(tok, " ");  /* no. of spaces before keyword */
    tok += len;

    if (tok[0] == '-')  /* is there a leading minus sign? */
    {
        /* first token is name of keyword to be deleted or renamed */
        *hdtype = -1;
        tok++;
        len = strspn(tok, " ");  /* no. of spaces before keyword */
        tok += len;
        if (len < 8)  /* not a blank name? */
        {
          len = strcspn(tok, " =");  /* length of name */
          if (len > 8)
            return(*status = BAD_KEYCHAR);

          strncat(card, tok, len);
          ffupch(card);
          if (fftkey(card, status) > 0)
              return(*status);      /* illegal chars in name */

          tok += len;
        }

        /* second token, if present, is the new name for the keyword */

        len = strspn(tok, " ");  /* no. of spaces before next token */
        tok += len;

        if (tok[0] == '\0' || tok[0] == '=')
            return(*status);  /* no second token */

        *hdtype = -2;
        len = strcspn(tok, " ");  /* length of new name */
        if (len > 8)
          return(*status = BAD_KEYCHAR);

        /* copy the new name to card + 40;  This is awkward, */
        /* but is consistent with the way the Fortran FITSIO works */
        strncpy(&card[40], tok, len);
        ffupch(&card[40]);
        fftkey(&card[40], status);
    }
    else  /* no negative sign at beginning of template */
    {
      /* get the keyword name token */

      len = strcspn(tok, " =");  /* length of keyword name */
      if (len > 8)
        return(*status = BAD_KEYCHAR);

      strncat(keyname, tok, len);
      ffupch(keyname);
      if (fftkey(keyname, status) > 0)
          return(*status);

      if (!strcmp(keyname, "END") )
      {
         strcpy(card, "END");
         *hdtype = 2;
         return(*status);
      }

      tok += len;
      len = strspn(tok, " =");  /* space between name and value */
      
      if (strcmp(keyname, "COMMENT") && strcmp(keyname, "HISTORY") )
      {
        /* Get the Value token */
        tok += len;

        if (*tok == '\'') /* is value enclosed in quotes? */
        {
          more = TRUE;
          while (more)
          {
            tok++;  /* temporarily move past the quote char */
            len = strcspn(tok, "'");  /* length of quoted string */
            tok--;
            strncat(value, tok, len + 2); 
 
            tok += len + 1;
            if (tok[0] != '\'')   /* check there is a closing quote */
              return(*status = NO_QUOTE);

            tok++;
            if (tok[0] != '\'')  /* 2 quote chars = literal quote */
              more = FALSE;
          }
        }
        else   /* not a quoted string value */
        {
          len = strcspn(tok, " /"); /* length of value string */

          strncat(value, tok, len);
          if (tok[0] != 'T' && tok[0] != 'F') /* not a logical value */
          {
            dval = strtod(value, &suffix); /* try to read value as number */

            if (*suffix != '\0' && *suffix != ' ' && *suffix != '/')
            { 
              /* value is not a number; must enclose it in quotes */
              strcpy(value, "'");
              strncat(value, tok, len);
              strcat(value, "'");

              /* the following useless statement stops the compiler warning */
              /* that dval is not used anywhere */
              if (dval == 0.)
                 len += (int) dval; 
            }
          }
          tok += len;
        }

        len = strspn(tok, " /"); /* no. of spaces between value and comment */
        tok += len;
      }
      else
      {
        *hdtype = 1;   /* simply append COMMENT and HISTORY keywords */
        strcpy(card, keyname);
        strncat(card, tok, 73);
        return(*status);
      }

      vlen = strlen(value);
      if (vlen > 0 && vlen < 10 && value[0] == '\'')
      {
          /* pad quoted string with blanks so it is at least 8 chars long */
          value[vlen-1] = '\0';
          strncat(value, "        ", 10 - vlen);
          strcat(&value[9], "'");
      }

      /* get the comment string */
      strncpy(comment, tok, 70);

      /* construct the complete FITS header card */
      ffmkky(keyname, value, comment, card);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffasfm(char *tform,    /* I - format code from the TFORMn keyword */
           int *datacode,  /* O - numerical datatype code */
           long *width,    /* O - width of the field, in chars */
           int *decimals,  /* O - number of decimal places (F, E, D format) */
           int *status)    /* IO - error status      */
{
/*
  parse the ASCII table TFORM column format to determine the data
  type, the field width, and number of decimal places (if relevant)
*/
    int ii;
    long longval;
    float fwidth;
    char *form, temp[FLEN_VALUE], message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    *datacode = 0;
    *width = 0;
    *decimals = 0;

    ii = 0;
    while (tform[ii] != 0 && tform[ii] == ' ') /* find first non-blank char */
         ii++;

    strcpy(temp, &tform[ii]); /* copy format string */
    ffupch(temp);     /* make sure it is in upper case */
    form = temp;      /* point to start of format string */


    if (form[0] == 0)
    {
        ffpmsg("Error: ASCII table TFORM code is blank");
        return(*status = BAD_TFORM);
    }

    /*-----------------------------------------------*/
    /*       determine default datatype code         */
    /*-----------------------------------------------*/
    if (form[0] == 'A')
        *datacode = TSTRING;
    else if (form[0] == 'I')
        *datacode = TLONG;
    else if (form[0] == 'F')
        *datacode = TFLOAT;
    else if (form[0] == 'E')
        *datacode = TFLOAT;
    else if (form[0] == 'D')
        *datacode = TDOUBLE;
    else
    {
        sprintf(message,
        "Illegal ASCII table TFORMn datatype: \'%s\'", tform);
        ffpmsg(message);
        return(*status = BAD_TFORM_DTYPE);
    }

    form++;  /* point to the start of field width */

    if (*datacode == TSTRING || *datacode == TLONG)
    { 
        /*-----------------------------------------------*/
        /*              A or I data formats:             */
        /*-----------------------------------------------*/

        if (ffc2ii(form, width, status) <= 0)  /* read the width field */
        {
            if (*width <= 0)
            {
                *width = 0;
                *status = BAD_TFORM;
            }
            else
            {                
                /* set to shorter precision if I4 or less */
                if (*width <= 4 && *datacode == TLONG)
                    *datacode = TSHORT;
            }
        }
    }
    else
    {  
        /*-----------------------------------------------*/
        /*              F, E or D data formats:          */
        /*-----------------------------------------------*/

        if (ffc2rr(form, &fwidth, status) <= 0) /* read ww.dd width field */
        {
           if (fwidth <= 0.)
            *status = BAD_TFORM;
          else
          {
            *width = (long) fwidth;  /* convert from float to long */

            if (*width > 7 && *(form-1) == 'F')
                *datacode = TDOUBLE;  /* type double if >7 digits */

            if (*width < 10)
                form = form + 1; /* skip 1 digit  */
            else
                form = form + 2; /* skip 2 digits */

            if (form[0] == '.') /* should be a decimal point here */
            {
                form++;  /*  point to start of decimals field */

                if (ffc2ii(form, &longval, status) <= 0) /* read decimals */
                {
                    *decimals = longval;  /* long to short convertion */

                    if (*decimals >= *width)  /* width < no. of decimals */
                        *status = BAD_TFORM; 
                }
            }
          }
        }
    }
    if (*status > 0)
    {
        *status = BAD_TFORM;
        sprintf(message,"Illegal ASCII table TFORMn code: \'%s\'", tform);
        ffpmsg(message);
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffbnfm(char *tform,     /* I - format code from the TFORMn keyword */
           int *datacode,   /* O - numerical datatype code */
           long *repeat,    /* O - repeat count of the field  */
           long *width,     /* O - width of the field, in chars */
           int *status)     /* IO - error status      */
{
/*
  parse the binary table TFORM column format to determine the data
  type, repeat count, and the field width (if it is an ASCII (A) field)
*/
    size_t ii, nchar;
    int variable, iread;
    char *form, temp[FLEN_VALUE], message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    *datacode = 0;
    *repeat = 0;
    *width = 0;

    nchar = strlen(tform);

    for (ii = 0; ii < nchar; ii++)
    {
        if (tform[ii] != ' ')     /* find first non-space char */
            break;
    }

    if (ii == nchar)
    {
        ffpmsg("Error: binary table TFORM code is blank (ffbnfm).");
        return(*status = BAD_TFORM);
    }

    strcpy(temp, &tform[ii]); /* copy format string */
    ffupch(temp);     /* make sure it is in upper case */
    form = temp;      /* point to start of format string */

    /*-----------------------------------------------*/
    /*       get the repeat count                    */
    /*-----------------------------------------------*/

    ii = 0;
    while(isdigit(form[ii]))
        ii++;   /* look for leading digits in the field */

    if (ii == 0)
        *repeat = 1;  /* no explicit repeat count */
    else
        sscanf(form,"%ld", repeat);  /* read repeat count */

    /*-----------------------------------------------*/
    /*             determine datatype code           */
    /*-----------------------------------------------*/

    form = form + ii;  /* skip over the repeat field */

    if (form[0] == 'P')
    {
        variable = 1;  /* this is a variable length column */
        *repeat = 1;   /* disregard any other repeat value */
        form++;        /* move to the next data type code char */
    }
    else
        variable = 0;

    if (form[0] == 'U')  /* internal code to signify unsigned integer */
    { 
        *datacode = TUSHORT;
        *width = 2;
    }
    else if (form[0] == 'I')
    {
        *datacode = TSHORT;
        *width = 2;
    }
    else if (form[0] == 'V') /* internal code to signify unsigned integer */
    {
        *datacode = TULONG;
        *width = 4;
    }
    else if (form[0] == 'J')
    {
        *datacode = TLONG;
        *width = 4;
    }
    else if (form[0] == 'E')
    {
        *datacode = TFLOAT;
        *width = 4;
    }
    else if (form[0] == 'D')
    {
        *datacode = TDOUBLE;
        *width = 8;
    }
    else if (form[0] == 'A')
    {
        *datacode = TSTRING;

        /*
          the following code is used to support the non-standard
          datatype of the form rAw where r = total width of the field
          and w = width of fixed-length substrings within the field.
        */
        iread = 0;
        if (form[1] != 0)
            iread = sscanf(&form[1],"%ld", width);

        if (iread != 1)
            *width = *repeat;  
    }
    else if (form[0] == 'L')
    {
        *datacode = TLOGICAL;
        *width = 1;
    }
    else if (form[0] == 'X')
    {
        *datacode = TBIT;
        *width = 1;
    }
    else if (form[0] == 'B')
    {
        *datacode = TBYTE;
        *width = 1;
    }
    else if (form[0] == 'C')
    {
        *datacode = TCOMPLEX;
        *width = 8;
    }
    else if (form[0] == 'M')
    {
        *datacode = TDBLCOMPLEX;
        *width = 16;
    }
    else
    {
        sprintf(message,
        "Illegal binary table TFORMn datatype: \'%s\' ", tform);
        ffpmsg(message);
        return(*status = BAD_TFORM_DTYPE);
    }

    if (variable)
        *datacode = *datacode * (-1); /* flag variable cols w/ neg type code */

    return(*status);
}
/*--------------------------------------------------------------------------*/
void ffcfmt(char *tform,    /* value of an ASCII table TFORMn keyword */
            char *cform)    /* equivalent format code in C language syntax */
/*
  convert the FITS format string for an ASCII Table extension column into the
  equivalent C format string that can be used in a printf statement.
*/
{
    int ii;

    cform[0] = '\0';
    ii = 0;
    while (tform[ii] != 0 && tform[ii] == ' ') /* find first non-blank char */
         ii++;

    if (tform[ii] == 0)
        return;    /* input format string was blank */

    cform[0] = '%';  /* start the format string */

    strcpy(&cform[1], &tform[ii + 1]); /* append the width and decimal code */


    if (tform[ii] == 'A')
        strcat(cform, "s");
    else if (tform[ii] == 'I')
        strcat(cform, ".0lf");  /*  0 precision to suppress decimal point */
    if (tform[ii] == 'F')
        strcat(cform, "lf");
    if (tform[ii] == 'E')
        strcat(cform, "lE");
    if (tform[ii] == 'D')
        strcat(cform, "lE");

    return;
}
/*--------------------------------------------------------------------------*/
int ffgcno( fitsfile *fptr,  /* I - FITS file pionter                       */
            int  casesen,    /* I - case sensitive string comparison? 0=no  */
            char *templt,    /* I - input name of column (w/wildcards)      */
            int  *colnum,    /* O - number of the named column; 1=first col */
            int  *status)    /* IO - error status                           */
/*
  Determine the column number corresponding to an input column name.
  The first column of the table = column 1;  
  This supports the * and ? wild cards in the input template.
*/
{
    char colname[FLEN_VALUE];  /*  temporary string to hold column name  */

    ffgcnn(fptr, casesen, templt, colname, colnum, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcnn( fitsfile *fptr,  /* I - FITS file pointer                       */
            int  casesen,    /* I - case sensitive string comparison? 0=no  */
            char *templt,    /* I - input name of column (w/wildcards)      */
            char *colname,   /* O - full column name up to 68 + 1 chars long*/
            int  *colnum,    /* O - number of the named column; 1=first col */
            int  *status)    /* IO - error status                           */
/*
  Return the full column name and column number of the next column whose
  TTYPEn keyword value matches the input template string.
  The template may contain the * and ? wildcards.  Status = 237 is
  returned if the match is not unique.  If so, one may call this routine
  again with input status=237  to get the next match.  A status value of
  219 is returned when there are no more matching columns.
*/
{
    char errmsg[FLEN_ERRMSG];
    static int startcol;
    int tstatus, ii, founde, foundw, match, exact, unique;
    long ivalue;
    tcolumn *colptr;

    if (*status <= 0)
    {
        startcol = 0;   /* start search with first column */
        tstatus = 0;
    }
    else if (*status == COL_NOT_UNIQUE) /* start search from previous spot */
    {
        tstatus = COL_NOT_UNIQUE;
        *status = 0;
    }
    else
        return(*status);  /* bad input status value */

    colname[0] = 0;    /* initialize null return */
    *colnum = 0;

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)   /* rescan header to get col struct */
            return(*status);

    colptr = fptr->tableptr;   /* pointer to first column */
    colptr += (startcol);      /* offset to starting column */

    founde = FALSE;   /* initialize 'found exact match' flag */
    foundw = FALSE;   /* initialize 'found wildcard match' flag */

    for (ii = startcol; ii < fptr->tfield; ii++, colptr++)
    {
        ffcmps(templt, colptr->ttype, casesen, &match, &exact);
        if (match)
        {
            if (founde && exact)
            {
                /* warning: this is the second exact match we've found     */
                /*reset pointer to first match so next search starts there */
               startcol = *colnum;
               return(*status = COL_NOT_UNIQUE);
            }
            else if (founde)   /* a wildcard match */
            {
                /* already found exact match so ignore this non-exact match */
            }
            else if (exact)
            {
                /* this is the first exact match we have found, so save it. */
                strcpy(colname, colptr->ttype);
                *colnum = ii + 1;
                founde = TRUE;
            }
            else if (foundw)
            {
                /* we have already found a wild card match, so not unique */
                /* continue searching for other matches                   */
                unique = FALSE;
            }
            else
            {
               /* this is the first wild card match we've found. save it */
               strcpy(colname, colptr->ttype);
               *colnum = ii + 1;
               startcol = *colnum;
               foundw = TRUE;
               unique = TRUE;
            }
        }
    }

    /* OK, we've checked all the names now see if we got any matches */
    if (founde)
    {
        if (tstatus == COL_NOT_UNIQUE)  /* we did find 1 exact match but */
            *status = COL_NOT_UNIQUE;   /* there was a previous match too */
    }
    else if (foundw)
    {
        /* found one or more wildcard matches; report error if not unique */
       if (!unique || tstatus == COL_NOT_UNIQUE)
           *status = COL_NOT_UNIQUE;
    }
    else
    {
        /* didn't find a match; check if template is a positive integer */
        ffc2ii(templt, &ivalue, &tstatus);
        if (tstatus ==  0 && ivalue <= fptr->tfield && ivalue > 0)
        {
            *colnum = ivalue;

            colptr = fptr->tableptr;   /* pointer to first column */
            colptr += (ivalue - 1);    /* offset to correct column */
            strcpy(colname, colptr->ttype);
        }
        else
        {
            *status = COL_NOT_FOUND;
            if (tstatus != COL_NOT_UNIQUE)
            {
              sprintf(errmsg, "ffgcnn could not find column: %.45s", templt);
              ffpmsg(errmsg);
            }
        }
    }
    
    startcol = *colnum;  /* save pointer for next time */
    return(*status);
}
/*--------------------------------------------------------------------------*/
void ffcmps(char *templt,   /* I - input template (may have wildcards)      */
            char *colname,  /* I - full column name up to 68 + 1 chars long */
            int  casesen,   /* I - case sensitive string comparison? 1=yes  */
            int  *match,    /* O - do template and colname match? 1=yes     */
            int  *exact)    /* O - do strings exactly match, or wildcards   */
/*
  compare the template to the string and test if they match.
  The strings are limited to 68 characters or less (the max. length
  of a FITS string keyword value.  This routine reports whether
  the two strings match and whether the match is exact or
  involves wildcards.

  This algorithm is very similar to the way unix filename wildcards
  work except that this first treats a wild card as a literal character
  when looking for a match.  If there is no literal match, then
  it interpretes it as a wild card.  So the template 'AB*DE'
  is considered to be an exact rather than a wild card match to
  the string 'AB*DE'.  The '#' wild card in the template string will 
  match any consecutive string of decimal digits in the colname.
  
*/
{
    int found, t1, s1;
    char temp[FLEN_VALUE], col[FLEN_VALUE];

    *match = FALSE;
    *exact = TRUE;

    strncpy(temp, templt, FLEN_VALUE); /* copy strings to work area */
    strncpy(col, colname, FLEN_VALUE);
    temp[FLEN_VALUE -1] = '\0';  /* make sure strings are terminated */
    col[FLEN_VALUE -1]  = '\0';

    if (!casesen)
    {             /* convert both strings to uppercase before comparison */
        ffupch(temp);
        ffupch(col);
    }

    if (!strcmp(temp, col) )
    {
        *match = TRUE;     /* strings exactly match */
        return;
    }

    *exact = FALSE;    /* strings don't exactly match */

    t1 = 0;   /* start comparison with 1st char of each string */
    s1 = 0;

    while(1)  /* compare corresponding chars in each string */
    {
      if (temp[t1] == '\0' || col[s1] == '\0')
      { 
         /* completely scanned one or both strings so they match */
         *match = TRUE;
         return;
      }

      if (temp[t1] == col[s1] || (temp[t1] == '?' && col[s1] != ' ') )
      {
        s1++;  /* corresponding chars in the 2 strings match */
        t1++;  /* increment both pointers and loop back again */
      }
      else if (temp[t1] == '#' && isdigit(col[s1]) )
      {
        s1++;  /* corresponding chars in the 2 strings match */
        t1++;  /* increment both pointers */

        /* find the end of the string of digits */
        while (isdigit(col[s1]) ) 
            s1++;        
      }
      else if (temp[t1] == '*')
      {    
        /* get next char from template and look for it in the col name */
        t1++;
        if (temp[t1] == '\0' || temp[t1] == ' ')
        {
          /* reached end of template so strings match */
          *match = TRUE;
          return;
        }

        found = FALSE;
        while (col[s1] && !found)
        {
          if (temp[t1] == col[s1])
          {
            t1++;  /* found matching characters; incre both pointers */
            s1++;  /* and loop back to compare next chars */
            found = TRUE;
          }
          else
            s1++;  /* increment the column name pointer and try again */
        }

        if (!found)
          return;  /* hit end of column name and failed to find a match */
      }
      else
        return;   /* strings don't match */
    }
}
/*--------------------------------------------------------------------------*/
int ffgtcl( fitsfile *fptr,  /* I - FITS file pointer                       */
            int  colnum,     /* I - column number                           */
            int *typecode,   /* O - datatype code (21 = short, etc)         */
            long *repeat,    /* O - repeat count of field                   */
            long *width,     /* O - if ASCII, width of field or unit string */
            int  *status)    /* IO - error status                           */
/*
  Get Type of table column. 
  Returns the datatype code of the column, as well as the vector
  repeat count and (if it is an ASCII character column) the
  width of the field or a unit string within the field.  This supports the
  TFORMn = 'rAw' syntax for specifying arrays of substrings, so
  if TFORMn = '60A12' then repeat = 60 and width = 12.
*/
{
    tcolumn *colptr;
    int decims;

    if (*status > 0)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    if (colnum < 1 || colnum > fptr->tfield)
        return(*status = BAD_COL_NUM);

    colptr = fptr->tableptr;   /* pointer to first column */
    colptr += (colnum - 1);    /* offset to correct column */

    if (fptr->hdutype == ASCII_TBL)
    {
       ffasfm(colptr->tform, typecode, width, &decims, status);
       *repeat = 1;
    }
    else
    {
      *typecode = colptr->tdatatype;
      *width = colptr->twidth;
      *repeat = colptr->trepeat;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgacl( fitsfile *fptr,   /* I - FITS file pointer                      */
            int  colnum,      /* I - column number                          */
            char *ttype,      /* O - TTYPEn keyword value                   */
            long *tbcol,      /* O - TBCOLn keyword value                   */
            char *tunit,      /* O - TUNITn keyword value                   */
            char *tform,      /* O - TFORMn keyword value                   */
            double *tscal,    /* O - TSCALn keyword value                   */
            double *tzero,    /* O - TZEROn keyword value                   */
            char *tnull,      /* O - TNULLn keyword value                   */
            char *tdisp,      /* O - TDISPn keyword value                   */
            int  *status)     /* IO - error status                          */
/*
  get ASCII column keyword values
*/
{
    char name[FLEN_KEYWORD], comm[FLEN_COMMENT];
    tcolumn *colptr;
    int tstatus;

    if (*status > 0)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    if (colnum < 1 || colnum > fptr->tfield)
        return(*status = BAD_COL_NUM);

    /* get what we can from the column structure */

    colptr = fptr->tableptr;   /* pointer to first column */
    colptr += (colnum -1);     /* offset to correct column */

    strcpy(ttype, colptr->ttype);
    *tbcol = (colptr->tbcol) + 1;  /* first col is 1, not 0 */
    strcpy(tform, colptr->tform);
    *tscal = colptr->tscale;
    *tzero = colptr->tzero;
    strcpy(tnull, colptr->strnull);

    /* read keywords to get additional parameters */

    ffkeyn("TUNIT", colnum, name, status);
    tstatus = 0;
    *tunit = '\0';
    ffgkys(fptr, name, tunit, comm, &tstatus);

    ffkeyn("TDISP", colnum, name, status);
    tstatus = 0;
    *tdisp = '\0';
    ffgkys(fptr, name, tdisp, comm, &tstatus);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgbcl( fitsfile *fptr,   /* I - FITS file pointer                      */
            int  colnum,      /* I - column number                          */
            char *ttype,      /* O - TTYPEn keyword value                   */
            char *tunit,      /* O - TUNITn keyword value                   */
            char *dtype,      /* O - datatype char: I, J, E, D, etc.        */
            long *repeat,     /* O - vector column repeat count             */
            double *tscal,    /* O - TSCALn keyword value                   */
            double *tzero,    /* O - TZEROn keyword value                   */
            long *tnull,      /* O - TNULLn keyword value integer cols only */
            char *tdisp,      /* O - TDISPn keyword value                   */
            int  *status)     /* IO - error status                          */
/*
  get BINTABLE column keyword values
*/
{
    char name[FLEN_KEYWORD], comm[FLEN_COMMENT];
    tcolumn *colptr;
    int tstatus;

    if (*status > 0)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    if (colnum < 1 || colnum > fptr->tfield)
        return(*status = BAD_COL_NUM);

    /* get what we can from the column structure */

    colptr = fptr->tableptr;   /* pointer to first column */
    colptr += (colnum -1);     /* offset to correct column */

    strcpy(ttype, colptr->ttype);
    if (colptr->tdatatype < 0)  /* add the "P" prefix for */
        strcpy(dtype, "P");     /* variable length columns */
    else
        dtype[0] = 0;

    if      (abs(colptr->tdatatype) == TBIT)
        strcat(dtype, "X");
    else if (abs(colptr->tdatatype) == TBYTE)
        strcat(dtype, "B");
    else if (abs(colptr->tdatatype) == TLOGICAL)
        strcat(dtype, "L");
    else if (abs(colptr->tdatatype) == TSTRING)
        strcat(dtype, "A");
    else if (abs(colptr->tdatatype) == TSHORT)
        strcat(dtype, "I");
    else if (abs(colptr->tdatatype) == TLONG)
        strcat(dtype, "J");
    else if (abs(colptr->tdatatype) == TFLOAT)
        strcat(dtype, "E");
    else if (abs(colptr->tdatatype) == TDOUBLE)
        strcat(dtype, "D");
    else if (abs(colptr->tdatatype) == TCOMPLEX)
        strcat(dtype, "C");
    else if (abs(colptr->tdatatype) == TDBLCOMPLEX)
        strcat(dtype, "M");

    *repeat = colptr->trepeat;
    *tscal  = colptr->tscale;
    *tzero  = colptr->tzero;
    *tnull  = colptr->tnull;

    /* read keywords to get additional parameters */

    ffkeyn("TUNIT", colnum, name, status);
    tstatus = 0;
    *tunit = '\0';
    ffgkys(fptr, name, tunit, comm, &tstatus);

    ffkeyn("TDISP", colnum, name, status);
    tstatus = 0;
    *tdisp = '\0';
    ffgkys(fptr, name, tdisp, comm, &tstatus);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffghdn(fitsfile *fptr,   /* I - FITS file pointer                      */
            int *chdunum)    /* O - number of the CHDU; 1 = primary array  */
/*
  Return the number of the Current HDU in the FITS file.  The primary array
  is HDU number 1.  Note that this is one of the few cfitsio routines that
  does not return the error status value as the value of the function.
*/
{
    *chdunum = (fptr->curhdu) + 1;       

    return(*chdunum);
}
/*--------------------------------------------------------------------------*/
void ffghad(fitsfile *fptr,   /* I - FITS file pointer                     */
            long *chduaddr,   /* O - byte offset to beginning of CHDU      */
            long *nextaddr)   /* O - byte offset to beginning of next HDU  */
/*
  Return the address (= byte offset) in the FITS file to the beginning of
  the current HDU and the offset to the beginning of the next HDU following
  the CHDU.
*/
{
    *chduaddr = fptr->headstart[fptr->curhdu];       
    *nextaddr = fptr->headstart[(fptr->curhdu) + 1];       

    return;
}
/*--------------------------------------------------------------------------*/
int ffrhdu(fitsfile *fptr,    /* I - FITS file pointer */
           int *hdutype,      /* O - type of HDU       */
           int *status)       /* IO - error status     */
/*
  read the required keywords of the CHDU and initialize the corresponding
  structure elements that describe the format of the HDU
*/
{
    int ii, tstatus;
    char card[FLEN_CARD];
    char name[FLEN_KEYWORD], value[FLEN_VALUE], comm[FLEN_COMMENT];
    char xtension[FLEN_VALUE];

    if (*status > 0)
        return(*status);

    if (ffgrec(fptr, 1, card, status) > 0 )  /* get the 80-byte card */
    {
        ffpmsg("Cannot read first keyword in header (ffrhdu).");
        return(*status);
    }

    strncpy(name,card,8);  /* first 8 characters = the keyword name */
    name[8] = '\0';

    for (ii=7; ii >= 0; ii--)  /* replace trailing blanks with nulls */
    {
        if (name[ii] == ' ')
            name[ii] = '\0';
        else
            break;
    }
    if (ffpsvc(card, value, comm, status) > 0)   /* parse value and comment */
    {
        ffpmsg("Cannot read value of first  keyword in header (ffrhdu):");
        ffpmsg(card);
        return(*status);
    }

    if (!strcmp(name, "SIMPLE"))        /* this is the primary array */
    {
       ffpinit(fptr, status);           /* initialize the primary array */
       *hdutype = 0;
    }

    else if (!strcmp(name, "XTENSION"))   /* this is an XTENSION keyword */
    {
        if (ffc2s(value, xtension, status) > 0)  /* get the value string */
        {
            ffpmsg("Bad value string for XTENSION keyword:");
            ffpmsg(value);
            return(*status);
        }

        if (!strcmp(xtension, "TABLE"))
        {
            ffainit(fptr, status);       /* initialize the ASCII table */
            *hdutype = 1;
        }

        else if (!strcmp(xtension, "BINTABLE") ||
                 !strcmp(xtension, "A3DTABLE") ||
                 !strcmp(xtension, "3DTABLE") )
        {
            ffbinit(fptr, status);       /* initialize the binary table */
            *hdutype = 2;
        }

        else
        {
            tstatus = 0;
            ffpinit(fptr, &tstatus);       /* probably an IMAGE extension */

            if (tstatus == UNKNOWN_EXT)
                *hdutype = -1;       /* don't recognize this extension type */
            else
            {
                *status = tstatus;
                *hdutype = 0;
            }
        }
    }

    else     /*  not the start of a new extension */
    {
        if (card[0] == 0  ||
            card[0] == 10)     /* some editors append this character to EOF */
        {           
            *status = END_OF_FILE;
        }
        else
        {
          *status = UNKNOWN_REC;  /* found unknown type of record */
          ffpmsg
        ("Extension doesn't start with SIMPLE or XTENSION keyword. (ffrhdu)");
        }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpinit(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)          /* IO - error status     */
/*
  initialize the parameters defining the structure of the primary array
  or an Image extension 
*/
{
    int groups, tstatus, simple, bitpix, naxis, extend, nspace;
    int ttype, bytlen, ii;
    long naxes[999], pcount, gcount, npix, blank;
    double bscale, bzero;
    char comm[FLEN_COMMENT];
    tcolumn *colptr;

    if (*status > 0)
        return(*status);

    fptr->hdutype = IMAGE_HDU; /* primary array or IMAGE extension  */
    fptr->headend = 2000000000;  /* temporarily set huge header size  */

    groups = 0;
    tstatus = *status;

    /* get all the descriptive info about this HDU */
    ffgphd(fptr, 999, &simple, &bitpix, &naxis, naxes, &pcount, &gcount, 
           &extend, &bscale, &bzero, &blank, &nspace, status);

    if (*status == NOT_IMAGE)
        *status = tstatus;    /* ignore 'unknown extension type' error */
    else if (*status > 0)
        return(*status);

    /*
       the logical end of the header is 80 bytes before the current position, 
       minus any trailing blank keywords just before the END keyword.
    */
    fptr->headend = fptr->nextkey - (80 * (nspace + 1));

    /* the data unit begins at the beginning of the next logical block */
    fptr->datastart = ( (fptr->nextkey - 80) / 2880 + 1) * 2880;

    if (naxis > 0 && naxes[0] == 0)  /* test for 'random groups' */
    {
        tstatus = 0;
        if (ffgkyl(fptr, "GROUPS", &groups, comm, &tstatus))
            groups = 0;          /* GROUPS keyword not found */
    }

    if (bitpix == BYTE_IMG)   /* test  bitpix and set the datatype code */
    {
        ttype=TBYTE;
        bytlen=1;
    }
    else if (bitpix == SHORT_IMG)
    {
        ttype=TSHORT;
        bytlen=2;
    }
    else if (bitpix == LONG_IMG)
    {
        ttype=TLONG;
        bytlen=4;
    }
    else if (bitpix == FLOAT_IMG)
    {
        ttype=TFLOAT;
        bytlen=4;
    }
    else if (bitpix == DOUBLE_IMG)
    {
        ttype=TDOUBLE;
        bytlen=8;
    }
        
    /*   calculate the size of the primary array  */
    if (naxis == 0)
        npix = 0;
    else
    {
        if (groups)
            npix = 1;  /* NAXIS1 = 0 is a special flag for 'random groups' */
        else
            npix = naxes[0];

        for (ii=1; ii < naxis; ii++)
             npix = npix*naxes[ii];   /* calc number of pixels in the array */
    }

    /*
       now we know everything about the array; just fill in the parameters:
       the next HDU begins in the next logical block after the data
    */
    fptr->headstart[ fptr->curhdu + 1] =
         fptr->datastart + 
         ( (pcount + npix) * bytlen * gcount + 2879) / 2880 * 2880;

    /*
      initialize the fictitious heap starting address (immediately following
      the array data) and a zero length heap.  This is used to find the
      end of the data when checking the fill values in the last block. 
    */
    fptr->heapstart = (pcount + npix) * bytlen * gcount;
    fptr->heapsize = 0;

    if (naxis == 0)
    {
        fptr->rowlength = 0;    /* rows have zero length */
        fptr->tfield = 0;       /* table has no fields   */
        fptr->tableptr = 0;     /* set a null table structure pointer */
    }
    else
    {
      /*
        The primary array is actually interpreted as a binary table.  There
        are two columns: the first column contains the group parameters if any.
        The second column contains the primary array of data as a single vector
        column element. In the case of 'random grouped' format, each group
        is stored in a separate row of the table.
      */
        fptr->rowlength = (pcount + npix) * bytlen; /* total size of image */
        fptr->tfield = 2;  /* 2 fields: group parameters and the image */

        if (fptr->tableptr)
           free(fptr->tableptr); /* free memory for the old CHDU structure */

        colptr = (tcolumn *) calloc(2, sizeof(tcolumn) ) ;

        if (!colptr)
        {
          ffpmsg
          ("malloc failed to get memory for FITS array descriptors (ffpinit)");
          fptr->tableptr = 0;  /* set a null table structure pointer */
          return(*status = ARRAY_TOO_BIG);
        }

        /* copy the table structure address to the fitsfile structure */
        fptr->tableptr = colptr; 

        /* the first column represents the group parameters, if any */
        colptr->tbcol = 0;
        colptr->tdatatype = ttype;
        colptr->twidth = bytlen;
        colptr->trepeat = pcount;
        colptr->tscale = 1.;
        colptr->tzero = 0.;
        colptr->tnull = blank;

        colptr++;  /* increment pointer to the second column */

        /* the second column represents the image array */
        colptr->tbcol = pcount * bytlen; /* col starts after the group parms */
        colptr->tdatatype = ttype; 
        colptr->twidth = bytlen;
        colptr->trepeat = npix;
        colptr->tscale = bscale;
        colptr->tzero = bzero;
        colptr->tnull = blank;
    }

    /* reset next keyword pointer to the start of the header */
    fptr->nextkey = fptr->headstart[ fptr->curhdu ];

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffainit(fitsfile *fptr,      /* I - FITS file pointer */
            int *status)         /* IO - error status     */
{
/*
  initialize the parameters defining the structure of an ASCII table 
*/
    int  ii, nspace, tbcoln;
    long nrows, rowlen, pcount, tfield;
    tcolumn *colptr;
    char name[FLEN_KEYWORD], value[FLEN_VALUE], comm[FLEN_COMMENT];
    char message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    fptr->hdutype = ASCII_TBL;  /* set that this is an ASCII table */
    fptr->headend = 2000000000;       /* temporarily set huge header size  */

    /* get table parameters and test that the header is a valid: */
    if (ffgttb(fptr, &rowlen, &nrows, &pcount, &tfield, status) > 0)  
       return(*status);

    if (pcount != 0)
    {
       ffpmsg("PCOUNT keyword not equal to 0 in ASCII table (ffainit).");
       return(*status = BAD_PCOUNT);
    }

    fptr->rowlength = rowlen; /* store length of a row */
    fptr->tfield = tfield;    /* store number of table fields in each row */

    if (fptr->tableptr)
       free(fptr->tableptr); /* free memory for the old CHDU structure */

    /* mem for column structures ; space is initialized = 0  */
    colptr = (tcolumn *) calloc(tfield, sizeof(tcolumn) );
    if (!colptr)
        {
          ffpmsg
          ("malloc failed to get memory for FITS table descriptors (ffainit)");
          fptr->tableptr = 0;  /* set a null table structure pointer */
          return(*status = ARRAY_TOO_BIG);
        }

    /* copy the table structure address to the fitsfile structure */
    fptr->tableptr = colptr; 

    /*  initialize the table field parameters */
    for (ii = 0; ii < tfield; ii++, colptr++)
    {
        colptr->ttype[0] = '\0';  /* null column name */
        colptr->tscale = 1.;
        colptr->tzero  = 0.;
        colptr->strnull[0] = ASCII_NULL_UNDEFINED;  /* null value undefined */
        colptr->tbcol = -1;          /* initialize to illegal value */
        colptr->tdatatype = -9999;   /* initialize to illegal value */
    }

    /*
      Initialize the fictitious heap starting address (immediately following
      the table data) and a zero length heap.  This is used to find the
      end of the table data when checking the fill values in the last block. 
      There is no special data following an ASCII table.
    */
    fptr->heapstart = rowlen * nrows;
    fptr->heapsize = 0;

    /* now search for the table column keywords and the END keyword */

    for (nspace = 0, ii = 8; 1; ii++)  /* infinite loop  */
    {
        ffgkyn(fptr, ii, name, value, comm, status);

        if (*status == END_OF_FILE)
        {
            ffpmsg("END keyword not found in ASCII table header (ffainit).");
            return(*status = NO_END);
        }
        else if (*status > 0)
            return(*status);

        else if (name[0] == 'T')   /* keyword starts with 'T' ? */
            ffgtbp(fptr, name, value, status); /* test if column keyword */

        else if (!strcmp(name, "END"))  /* is this the END keyword? */
            break;

        if (!name[0] && !value[0] && !comm[0])  /* a blank keyword? */
            nspace++;

        else
            nspace = 0;
    }

    /* test that all the required keywords were found and have legal values */
    colptr = fptr->tableptr;
    for (ii = 0; ii < tfield; ii++, colptr++)
    {
        tbcoln = colptr->tbcol;  /* the starting column number (zero based) */

        if (colptr->tdatatype == -9999)
        {
            ffkeyn("TFORM", ii+1, name, status);  /* construct keyword name */
            sprintf(message,"Required %s keyword not found (ffainit).", name);
            ffpmsg(message);
            return(*status = NO_TFORM);
        }

        else if (tbcoln == -1)
        {
            ffkeyn("TBCOL", ii+1, name, status); /* construct keyword name */
            sprintf(message,"Required %s keyword not found (ffainit).", name);
            ffpmsg(message);
            return(*status = NO_TBCOL);
        }

        else if (fptr->rowlength != 0 && 
                (tbcoln < 0 || tbcoln >= fptr->rowlength ) )
        {
            ffkeyn("TBCOL", ii+1, name, status);  /* construct keyword name */
            sprintf(message,"Value of %s keyword out of range: %d (ffainit).",
            name, tbcoln);
            ffpmsg(message);
            return(*status = BAD_TBCOL);
        }

        else if (fptr->rowlength != 0 && 
                 tbcoln + colptr->twidth > fptr->rowlength )
        {
            sprintf(message,"Column %d is too wide to fit in table (ffainit)",
            ii+1);
            ffpmsg(message);
            sprintf(message, " TFORM = %s and NAXIS1 = %ld",
                    colptr->tform, fptr->rowlength);
            ffpmsg(message);
            return(*status = COL_TOO_WIDE);
        }
    }

    /*
      now we know everything about the table; just fill in the parameters:
      the 'END' record is 80 bytes before the current position, minus
      any trailing blank keywords just before the END keyword.
    */
    fptr->headend = fptr->nextkey - (80 * (nspace + 1));
 
    /* the data unit begins at the beginning of the next logical block */
    fptr->datastart = ( (fptr->nextkey - 80) / 2880 + 1) * 2880;

    /* the next HDU begins in the next logical block after the data  */
    fptr->headstart[ fptr->curhdu + 1] = fptr->datastart +
         ( (rowlen * nrows + 2879) / 2880 * 2880 );

    /* reset next keyword pointer to the start of the header */
    fptr->nextkey = fptr->headstart[ fptr->curhdu ];

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffbinit(fitsfile *fptr,     /* I - FITS file pointer */
            int *status)        /* IO - error status     */
{
/*
  initialize the parameters defining the structure of a binary table 
*/
    int  ii, nspace;
    long nrows, rowlen, tfield, pcount, totalwidth, strrepeat;
    tcolumn *colptr;
    char name[FLEN_KEYWORD], value[FLEN_VALUE], comm[FLEN_COMMENT];
    char message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);

    fptr->hdutype = BINARY_TBL;  /* set that this is a binary table */
    fptr->headend = 2000000000;       /* temporarily set huge header size  */

    /* get table parameters and test that the header is valid: */
    if (ffgttb(fptr, &rowlen, &nrows, &pcount, &tfield, status) > 0)
       return(*status);

    fptr->rowlength = rowlen; /* store length of a row */
    fptr->tfield = tfield;   /* store number of table fields in each row */

    if (fptr->tableptr)
       free(fptr->tableptr); /* free memory for the old CHDU structure */

   /* mem for column structures ; space is initialized = 0  */
    colptr = (tcolumn *) calloc(tfield, sizeof(tcolumn) );
    if (!colptr)
    {
        ffpmsg
        ("malloc failed to get memory for FITS table descriptors (ffbinit)");
        fptr->tableptr = 0;  /* set a null table structure pointer */
        return(*status = ARRAY_TOO_BIG);
    }

    /* copy the table structure address to the fitsfile structure */
    fptr->tableptr = colptr; 

    /*  initialize the table field parameters */
    for (ii = 0; ii < tfield; ii++, colptr++)
    {
        colptr->ttype[0] = '\0';  /* null column name */
        colptr->tscale = 1.;
        colptr->tzero  = 0.;
        colptr->tnull  = NULL_UNDEFINED; /* (integer) null value undefined */
        colptr->tdatatype = -9999;   /* initialize to illegal value */
        colptr->trepeat = 1;
        colptr->strnull[0] = '\0'; /* for ASCII string columns (TFORM = rA) */
    }

    /*
      Initialize the heap starting address (immediately following
      the table data) and the size of the heap.  This is used to find the
      end of the table data when checking the fill values in the last block. 
    */
    fptr->heapstart = rowlen * nrows;
    fptr->heapsize = pcount;

    /* now search for the table column keywords and the END keyword */

    for (nspace = 0, ii = 8; 1; ii++)  /* infinite loop  */
    {
        ffgkyn(fptr, ii, name, value, comm, status);

        if (*status == END_OF_FILE)
        {
            ffpmsg("END keyword not found in binary table header (ffbinit).");
            return(*status = NO_END);
        }
        else if (*status > 0)
            return(*status);

        else if (name[0] == 'T')   /* keyword starts with 'T' ? */
            ffgtbp(fptr, name, value, status); /* test if column keyword */

        else if (!strcmp(name, "END"))  /* is this the END keyword? */
            break;


        if (!name[0] && !value[0] && !comm[0])  /* a blank keyword? */
            nspace++;

        else
            nspace = 0; /* reset number of consecutive spaces before END */
    }

    /* test that all the required keywords were found and have legal values */
    colptr = fptr->tableptr;  /* set pointer to first column */

    for (ii = 0; ii < tfield; ii++, colptr++)
    {
        if (colptr->tdatatype == -9999)
        {
            ffkeyn("TFORM", ii+1, name, status);  /* construct keyword name */
            sprintf(message,"Required %s keyword not found (ffbinit).", name);
            ffpmsg(message);
            return(*status = NO_TFORM);
        }
    }

    /*
      now we know everything about the table; just fill in the parameters:
      the 'END' record is 80 bytes before the current position, minus
      any trailing blank keywords just before the END keyword.
    */
    fptr->headend = fptr->nextkey - (80 * (nspace + 1));
 
    /* the data unit begins at the beginning of the next logical block */
    fptr->datastart = ( (fptr->nextkey - 80) / 2880 + 1) * 2880;

    /* the next HDU begins in the next logical block after the data  */
    fptr->headstart[ fptr->curhdu + 1] = fptr->datastart +
         ( (rowlen * nrows + pcount + 2879) / 2880 * 2880 );

    /* determine the byte offset to the beginning of each column */
    ffgtbc(fptr, &totalwidth, status);

    if (totalwidth != rowlen)
    {
        sprintf(message,
        "NAXIS1 = %ld is not equal to the sum of column widths: %ld", 
        rowlen, totalwidth);
        ffpmsg(message);
        *status = BAD_ROW_WIDTH;
    }

    /* reset next keyword pointer to the start of the header */
    fptr->nextkey = fptr->headstart[ fptr->curhdu ];

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgabc(int tfields,     /* I - number of columns in the table           */
           char **tform,    /* I - value of TFORMn keyword for each column  */
           int space,       /* I - number of spaces to leave between cols   */
           long *rowlen,    /* O - total width of a table row               */
           long *tbcol,     /* O - starting byte in row for each column     */
           int *status)     /* IO - error status                            */
/*
  calculate the starting byte offset of each column of an ASCII table
  and the total length of a row, in bytes.  The input space value determines
  how many blank spaces to leave between each column (1 is recommended).
*/
{
    int ii, datacode, decims;
    long width;

    if (*status > 0)
        return(*status);

    *rowlen=0;
    tbcol[0] = 1;

    for (ii = 0; ii < tfields; ii++)
    {
        tbcol[ii] = *rowlen + 1;    /* starting byte in row of column */

        ffasfm(tform[ii], &datacode, &width, &decims, status);

        *rowlen += (width + space);  /* total length of row */
    }

    *rowlen -= space;  /*  don't add space after the last field */

    return (*status);
}
/*--------------------------------------------------------------------------*/
int ffgtbc(fitsfile *fptr,    /* I - FITS file pointer          */
           long *totalwidth,  /* O - total width of a table row */
           int *status)       /* IO - error status              */
{
/*
  calculate the starting byte offset of each column of a binary table.
  Use the values of the datatype code and repeat counts in the
  column structure. Return the total length of a row, in bytes.
*/
    int tfields, ii;
    long nbytes;
    tcolumn *colptr;

    if (*status > 0)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    tfields = fptr->tfield;
    colptr = fptr->tableptr;  /* point to first column structure */

    *totalwidth = 0;

    for (ii = 0; ii < tfields; ii++, colptr++)
    {
        colptr->tbcol = *totalwidth;  /* byte offset in row to this column */

        if (colptr->tdatatype == TSTRING)
        {
            nbytes = colptr->trepeat;   /* one byte per char */
        }
        else if (colptr->tdatatype == TBIT)
        {
            nbytes = (colptr->trepeat + 7) / 8;
        }
        else if (colptr->tdatatype > 0)
        {
            nbytes = colptr->trepeat * (colptr->tdatatype / 10);
        }
        else   /* this is a variable length descriptor (neg. tdatatype) */
            nbytes = 8;

        *totalwidth = *totalwidth + nbytes;
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgtbp(fitsfile *fptr,     /* I - FITS file pointer   */
           char *name,         /* I - name of the keyword */
           char *value,        /* I - value string of the keyword */
           int *status)        /* IO - error status       */
{
/*
  Get TaBle Parameter.  The input keyword name begins with the letter T.
  Test if the keyword is one of the table column definition keywords
  of an ASCII or binary table. If so, decode it and update the value 
  in the structure.
*/
    int tstatus, datacode, decimals;
    long width, repeat, nfield, ivalue;
    double dvalue;
    char tvalue[FLEN_VALUE];
    char message[FLEN_ERRMSG];
    tcolumn *colptr;

    if (*status > 0)
        return(*status);

    tstatus = 0;

    if (!strncmp(name, "TTYPE", 5) ||
        !strncmp(name, "TFORM", 5) ||
        !strncmp(name, "TBCOL", 5) ||
        !strncmp(name, "TSCAL", 5) ||
        !strncmp(name, "TZERO", 5) ||
        !strncmp(name, "TNULL", 5) ||
        !strncmp(name, "THEAP", 5) )
    {

    if (!strncmp(name, "THEAP", 5))
    {
        if (fptr->hdutype == ASCII_TBL)  /* ASCII table */
            return(*status);  /* ASCII tables don't have a heap */ 

        if (ffc2ii(value, &ivalue, status) > 0) 
        {
            sprintf(message,
            "Error reading value of %s as an integer: %s", name, value);
            ffpmsg(message);
            return(*status);
        }
        fptr->heapstart = ivalue; /* the starting byte of the heap */
        return(*status);
    }

    if( ffc2ii(&name[5], &nfield, &tstatus) > 0) /* read index no. */
        return(*status);    /* must not be an indexed keyword */

    if (nfield < 1 || nfield > fptr->tfield )  /* index out of range */
        return(*status);

    colptr = fptr->tableptr;        /* get pointer to columns */
    colptr = colptr + nfield - 1;   /* point to the correct column */

    if (!strncmp(name, "TTYPE", 5))
    {
        if (ffc2s(value, tvalue, &tstatus) > 0)  /* remove quotes */
            return(*status);

        strcpy(colptr->ttype, tvalue);  /* copy col name to structure */
    }

    else if (!strncmp(name, "TFORM", 5))
    {
        if (ffc2s(value, tvalue, &tstatus) > 0)  /* remove quotes */
            return(*status);

        strncpy(colptr->tform, tvalue, 9);  /* copy TFORM to structure */
        colptr->tform[9] = '\0';            /* make sure it is terminated */

        if (fptr->hdutype == ASCII_TBL)  /* ASCII table */
        {
          if (ffasfm(tvalue, &datacode, &width, &decimals, status) > 0)
              return(*status);  /* bad format code */

          colptr->tdatatype = TSTRING; /* store datatype code */
          colptr->trepeat = 1;      /* field repeat count == 1 */
          colptr->twidth = width;   /* the width of the field, in bytes */
        }
        else  /* binary table */
        {
          if (ffbnfm(tvalue, &datacode, &repeat, &width, status) > 0)
              return(*status);  /* bad format code */

          colptr->tdatatype = datacode; /* store datatype code */
          colptr->trepeat = repeat;     /* field repeat count  */
          colptr->twidth = width;   /*  width of a unit value in chars */
        }
    }

    else if (!strncmp(name, "TBCOL", 5))
    {
        if (fptr->hdutype == BINARY_TBL)
            return(*status);  /* binary tables don't have TBCOL keywords */

        if (ffc2ii(value, &ivalue, status) > 0)
        {
            sprintf(message,
            "Error reading value of %s as an integer: %s", name, value);
            ffpmsg(message);
            return(*status);
        }
        colptr->tbcol = ivalue - 1; /* convert to zero base */
    }
    else if (!strncmp(name, "TSCAL", 5))
    {
        if (ffc2dd(value, &dvalue, status) > 0)
        {
            sprintf(message,
            "Error reading value of %s as a double: %s", name, value);
            ffpmsg(message);
            return(*status);
        }
        colptr->tscale = dvalue;
    }
    else if (!strncmp(name, "TZERO", 5))
    {
        if (ffc2dd(value, &dvalue, status) > 0)
        {
            sprintf(message,
            "Error reading value of %s as a double: %s", name, value);
            ffpmsg(message);
            return(*status);
        }
        colptr->tzero = dvalue;
    }
    else if (!strncmp(name, "TNULL", 5))
    {
        if (fptr->hdutype == ASCII_TBL)  /* ASCII table */
        {
            if (ffc2s(value, tvalue, &tstatus) > 0)  /* remove quotes */
                return(*status);

            strncpy(colptr->strnull, tvalue, 17);  /* copy TNULL string */
            colptr->strnull[17] = '\0';  /* terminate the strnull field */

        }
        else  /* binary table */
        {
            if (ffc2ii(value, &ivalue, status) > 0) 
            {
                sprintf(message,
                "Error reading value of %s as an integer: %s", name, value);
                ffpmsg(message);
                return(*status);
            }
            colptr->tnull = ivalue; /* null value for integer column */
        }
    }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcpr( fitsfile *fptr, /* I - FITS file pointer                        */
        int colnum,     /* I - column number (1 = 1st column of table)      */
        long firstrow,  /* I - first row (1 = 1st row of table)             */
        long firstelem, /* I - first element within vector (1 = 1st)        */
        long nelem,     /* I - number of elements to read or write          */
        int writemode,  /* I - = 1 if writing data, = 0 if reading data     */
        double *scale,  /* O - FITS scaling factor (TSCALn keyword value)   */
        double *zero,   /* O - FITS scaling zero pt (TZEROn keyword value)  */
        char *tform,    /* O - ASCII column format: value of TFORMn keyword */
        long *twidth,   /* O - width of ASCII column (characters)           */
        int *tcode,     /* O - column datatype code: I*4=41, R*4=42, etc    */
        int *maxelem,   /* O - max number of elements that fit in buffer    */
        long *startpos, /* O - offset in file to starting row & column      */
        long *elemnum,  /* O - starting element number ( 0 = 1st element)   */
        long *incre,    /* O - byte offset between elements within a row    */
        long *repeat,   /* O - number of elements in a row (vector column)  */
        long *rowlen,   /* O - length of a row, in bytes                    */
        int  *hdutype,  /* O - HDU type: 0, 1, 2 = primary, table, bintable */
        long *tnull,    /* O - null value for integer columns               */
        char *snull,    /* O - null value for ASCII table columns           */
        int *status)    /* IO - error status                                */
/*
  Get Column PaRameters, and test starting row and element numbers for 
  validity.  
*/
{
    int nulpos;
    long datastart, tbcol;
    char message[81];
    tcolumn *colptr;

    /* rescan header if data structure is undefined */
    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               
            return(*status);

    /* Do sanity check of input parameters */
    if (firstrow < 1)
    {
        sprintf(message, "Starting row number is out of range: %ld",
                firstrow);
        ffpmsg(message);
        return(*status = BAD_ROW_NUM);
    }
    else if (fptr->hdutype != ASCII_TBL && firstelem < 1)
    {
        sprintf(message, "Starting element number is out of range: %ld",
                firstelem);
        ffpmsg(message);
        return(*status = BAD_ELEM_NUM);
    }
    else if (nelem < 0)
    {
        sprintf(message, "Negative no. of elements to read or write: %ld",
                nelem);
        ffpmsg(message);
        return(*status = NEG_BYTES);
    }
    else if (colnum < 1 || colnum > fptr->tfield)
    {
        sprintf(message, "Specified column number is out of range: %d",
                colnum);
        ffpmsg(message);
        return(*status = BAD_COL_NUM);
    }

    /*  copy relevant parameters from the structure */

    *hdutype = fptr->hdutype;      /* image, ASCII table, or BINTABLE  */
    *rowlen   = fptr->rowlength;   /* width of the table, in bytes     */
    datastart = fptr->datastart;   /* offset in file to start of table */

    colptr  = fptr->tableptr;    /* point to first column */
    colptr += (colnum - 1);      /* offset to correct column structure */

    *scale    = colptr->tscale;  /* value scaling factor;    default = 1.0 */
    *zero     = colptr->tzero;   /* value scaling zeropoint; default = 0.0 */
    *tnull    = colptr->tnull;   /* null value for integer columns         */
    tbcol     = colptr->tbcol;   /* offset to start of column within row   */
    *twidth   = colptr->twidth;  /* width of a single datum, in bytes      */
    *incre    = colptr->twidth;  /* increment between datums, in bytes     */
    *tcode    = colptr->tdatatype;
    *repeat   = colptr->trepeat;

    strcpy(tform, colptr->tform);    /* value of TFORMn keyword            */
    strcpy(snull, colptr->strnull);  /* null value for ASCII table columns */

    if (*hdutype == ASCII_TBL && snull[0] == '\0')
    {
     /* In ASCII tables, a null value is equivalent to all spaces */

       strcpy(snull, "                 ");   /* maximum of 17 spaces */
       nulpos = minvalue(17, *twidth);         /* truncate to width of column */
       snull[nulpos] = '\0';
    }

    /* Special case: interprete 'X' column as 'B' */
    if (abs(*tcode) == TBIT)
    {
        *tcode  = *tcode / TBIT * TBYTE;
        *repeat = (*repeat + 7) / 8;
    }

    /* Special case: support the 'rAw' format in BINTABLEs */
    if (*hdutype == BINARY_TBL && *tcode == TSTRING)
       *repeat = *repeat / *twidth;  /* repeat = # of unit strings in field */

    if (*hdutype == ASCII_TBL)
        *elemnum = 0;   /* ASCII tables don't have vector elements */
    else
        *elemnum = firstelem - 1;

    /* interprete complex and double complex as pairs of floats or doubles */
    if (abs(*tcode) >= TCOMPLEX)
    {
        if (*tcode > 0)
          *tcode = (*tcode + 1) / 2;
        else
          *tcode = (*tcode - 1) / 2;

        *repeat  = *repeat * 2;
        *twidth  = *twidth / 2;
        *incre   = *incre  / 2;
    }

    /* calculate no. of pixels that fit in buffer */
    /* allow for case where longs or floats are 8 bytes long */
    if (abs(*tcode) == TLONG)
       *maxelem = DBUFFSIZE / sizeof(long);
    else if (abs(*tcode) == TFLOAT)
       *maxelem = DBUFFSIZE / sizeof(float);
    else if (abs(*tcode) == TDOUBLE)
       *maxelem = DBUFFSIZE / sizeof(double);
    else if (abs(*tcode) == TSTRING)
       *maxelem = (DBUFFSIZE - 1)/ *twidth; /* leave room for final \0 */
    else
       *maxelem = DBUFFSIZE / *twidth; 

    /* calc starting byte position to 1st element of col  */
    /*  (this does not apply to variable length columns)  */
    *startpos = datastart + ((firstrow - 1) * *rowlen) + tbcol;

    if (*hdutype == IMAGE_HDU && writemode) /*  Primary VISArray or IMAGE */
    { /*
        For primary arrays, set the repeat count greater than the total
        number of pixels to be written.  This prevents an out-of-range
        error message in cases where the final image array size is not
        yet known or defined.
      */
        *repeat = *elemnum + nelem; 
    }
    else if (*tcode > 0)     /*  Fixed length table column  */
    {
        if (*elemnum >= *repeat)
        {
            sprintf(message, "Starting element number is out of range: %ld",
                    firstelem);
            ffpmsg(message);
            return(*status = BAD_ELEM_NUM);
        }
        else if (*repeat == 1 && nelem > 1)
        { /*
            When accessing a scalar column, fool the calling routine into
            thinking that this is a vector column with very big elements.
            This allows multiple values (up to the maxelem number of elements
            that will fit in the buffer) to be read or written with a single
            routine call, which increases the efficiency.
          */           
            *incre = *rowlen;
            *repeat = nelem;
        }
    }
    else    /*  Variable length Binary Table column */
    {
      *tcode *= (-1);  

      if (writemode)     /* return next empty heap address for writing */
      {
        *repeat = nelem + *elemnum; /* total no. of elements in the field */

        /*  calculate starting position (for writing new data) in the heap */
        *startpos = datastart + fptr->heapstart + fptr->heapsize;

        /*  write the descriptor into the fixed length part of table */
        ffpdes(fptr, colnum, firstrow, *repeat, fptr->heapsize, status);

        /* increment the address to the next empty heap position */
        fptr->heapsize += (*repeat * (*incre)); 
      }
      else    /*  get the read start position in the heap */
      {
        ffgdes(fptr, colnum, firstrow, repeat, startpos, status);

        if (colptr->tdatatype == -TBIT)
            *repeat = (*repeat + 7) / 8;  /* convert from bits to bytes */

        if (*elemnum >= *repeat)
        {
            sprintf(message, "Starting element number is out of range: %ld",
                    firstelem);
            ffpmsg(message);
            return(*status = BAD_ELEM_NUM);
        }

        *startpos = *startpos + datastart + fptr->heapstart;
      }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgdes(fitsfile *fptr, /* I - FITS file pointer                         */
           int colnum,     /* I - column number (1 = 1st column of table)   */
           long rownum,    /* I - row number (1 = 1st row of table)         */
           long *length,   /* O - number of elements in the row             */
           long *heapaddr, /* O - heap pointer to the data                  */
           int *status)    /* IO - error status                             */
/*
  get (read) the variable length vector descriptor from the table.
*/
{
    long bytepos, descript[2];
    tcolumn *colptr;

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    colptr = fptr->tableptr;  /* point to first column structure */
    colptr += (colnum - 1);   /* offset to the correct column */

    if (colptr->tdatatype >= 0)
        *status = NOT_VARI_LEN;

    else
    {
        bytepos = fptr->datastart + (rownum - 1) * fptr->rowlength +
                  colptr->tbcol;

        ffgi4b(fptr, bytepos, 2, 4, descript, status); /* read the descriptor */

        *length = descript[0];   /* 1st word of descriptor is the length  */
        *heapaddr = descript[1]; /* 2nd word of descriptor is the address */
     }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpdes(fitsfile *fptr, /* I - FITS file pointer                         */
           int colnum,     /* I - column number (1 = 1st column of table)   */
           long rownum,    /* I - row number (1 = 1st row of table)         */
           long length,    /* I - number of elements in the row             */
           long heapaddr,  /* I - heap pointer to the data                  */
           int *status)    /* IO - error status                             */
/*
  put (write) the variable length vector descriptor to the table.
*/
{
    long bytepos, descript[2];
    tcolumn *colptr;

    if (*status > 0)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    colptr = fptr->tableptr;  /* point to first column structure */
    colptr += (colnum - 1);   /* offset to the correct column */

    if (colptr->tdatatype >= 0)
        *status = NOT_VARI_LEN;

    else
    {
        bytepos = fptr->datastart + (rownum - 1) * fptr->rowlength +
                  colptr->tbcol;

        ffmbyt(fptr, bytepos, IGNORE_EOF, status); /* move to element */

        descript[0] = length;   /* 1st word of descriptor is the length  */
        descript[1] = heapaddr; /* 2nd word of descriptor is the address */
 
        ffpi4b(fptr, 2, 4, descript, status); /* write the descriptor */
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffchdu(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)         /* IO - error status     */
{
/*
  close the current HDU.  If we have write access to the file, then:
    - write the END keyword and pad header with blanks if necessary
    - check the data fill values, and rewrite them if not correct
*/
    long pcount;
    char comm[FLEN_COMMENT], message[FLEN_ERRMSG], valstring[FLEN_VALUE];
    char card[FLEN_CARD];

    if (fptr->writemode == 1)  /* write access to the file? */
    {
        /* if data has been written to variable length columns in a  */
        /* binary table, then we may need to update the PCOUNT value */
        if (fptr->heapsize > 0)
        {
          ffgkyj(fptr, "PCOUNT", &pcount, comm, status);
          if (fptr->heapsize > pcount)
          {
            /* would be simpler to just call ffmkyj here, but this */
            /* would force linking in all the modkey & putkey routines */
            sprintf(valstring, "%ld", fptr->heapsize);
            ffmkky("PCOUNT", valstring, comm, card);
            ffmkey(fptr, card, status);
          }

          ffuptf(fptr, status);  /* update the variable length TFORM values */
        }

        ffrdef(fptr, status);  /* scan header to redefine structure */
        ffpdfl(fptr, status);  /* insure correct data file values */
    }

    free(fptr->tableptr);      /* free memory for the CHDU structure */
    fptr->tableptr = 0;

    if (*status > 0)
    {
        sprintf(message,
        "Error while closing HDU number %d (ffchdu).", fptr->curhdu);
        ffpmsg(message);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffuptf(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)         /* IO - error status     */
/*
  Update the value of the TFORM keywords for the variable length array
  columns to make sure they all have the form 1Px(len) or Px(len) where
  'len' is the maximum length of the vector in the table (e.g., '1PE(400)')
*/
{
    int ii;
    long tflds, naxis2, maxlen, jj, length, addr;
    char comment[FLEN_COMMENT], keyname[FLEN_KEYWORD];
    char tform[FLEN_VALUE], newform[FLEN_VALUE], lenval[40];
    char card[FLEN_CARD];
    char message[FLEN_ERRMSG];

    ffgkyj(fptr, "TFIELDS", &tflds, comment, status);
    ffgkyj(fptr, "NAXIS2", &naxis2, comment, status);

    for (ii = 1; ii <= tflds; ii++)        /* loop over all the columns */
    {
      ffkeyn("TFORM", ii, keyname, status);          /* construct name */
      if (ffgkys(fptr, keyname, tform, comment, status) > 0)
      {
        sprintf(message,
        "Error while updating variable length vector TFORMn values (ffuptf).");
        ffpmsg(message);
        return(*status);
      }

      /* is this a variable array length column ? */
      if (tform[0] == 'P' || tform[1] == 'P')
      {
        if (strlen(tform) < 5)  /* is maxlen field missing? */
        {
          /* get the max length */
          maxlen = 0;
          for (jj=1; jj <= naxis2; jj++)
          {
            ffgdes(fptr, ii, jj, &length, &addr, status);
            maxlen = maxvalue(maxlen, length);
          }

          /* construct the new keyword value */
          strcpy(newform, "'");
          strcat(newform, tform);
          sprintf(lenval, "(%d)", maxlen);
          strcat(newform,lenval);
          while(strlen(newform) < 9)
             strcat(newform," ");   /* append spaces 'till length = 8 */
          strcat(newform,"'" );     /* append closing parenthesis */

          /* would be simpler to just call ffmkyj here, but this */
          /* would force linking in all the modkey & putkey routines */
          ffmkky(keyname, newform, comment, card);  /* make new card */
          ffmkey(fptr, card, status);   /* replace last read keyword */
       }
      }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffrdef(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)         /* IO - error status     */
/*
  ReDEFine the structure of a data unit.  This routine re-reads
  the CHDU header keywords to determine the structure and length of the
  current data unit.  This redefines the start of the next HDU.
*/
{
    int dummy;

    if (*status <= 0 && fptr->writemode == 1) /* write access to the file? */
        if (ffwend(fptr, status) <= 0)     /* rewrite END keyword and fill */
            ffrhdu(fptr, &dummy, status);  /* re-scan the header keywords  */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffhdef(fitsfile *fptr,      /* I - FITS file pointer                    */
           int morekeys,        /* I - reserve space for this many keywords */
           int *status)         /* IO - error status                        */
/*
  based on the number of keywords which have already been written,
  plus the number of keywords to reserve space for, we then can
  define where the data unit should start (it must start at the
  beginning of a 2880-byte logical block).

  This routine will only have any effect if the starting location of the
  data unit following the header is not already defined.  In any case,
  it is always possible to add more keywords to the header even if the
  data has already been written.  It is just more efficient to reserve
  the space in advance.
*/
{
    if (*status > 0 || morekeys < 1)
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
    {
      ffrdef(fptr, status);
      fptr->datastart = ((fptr->headend + (morekeys * 80)) / 2880 + 1) * 2880;
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffwend(fitsfile *fptr,       /* I - FITS file pointer */
            int *status)         /* IO - error status     */
/*
  write the END card and following fill (space chars) in the current header
*/
{
    int ii, tstatus;
    long endpos, nspace;
    char blankkey[FLEN_CARD], endkey[FLEN_CARD], keyrec[FLEN_CARD];

    if (*status > 0)
        return(*status);

    endpos = fptr->headend;

    /*  calc the data starting position if not currently defined */
    if (fptr->datastart == DATA_UNDEFINED)
        fptr->datastart = ( endpos / 2880 + 1 ) * 2880;

    /* calculate the number of blank keyword slots in the header */
    nspace = ( fptr->datastart - endpos ) / 80;

    /* construct a blank and END keyword (80 spaces )  */
    strcpy(blankkey, "                                        ");
    strcat(blankkey, "                                        ");
    strcpy(endkey, "END                                     ");
    strcat(endkey, "                                        ");
  
    /* check if header is already correctly terminated with END and fill */
    tstatus=0;
    ffmbyt(fptr, endpos, REPORT_EOF, &tstatus); /* move to header end */
    for (ii=0; ii < nspace; ii++)
    {
        ffgbyt(fptr, 80, keyrec, &tstatus);  /* get next keyword */
        if (strncmp(keyrec, blankkey, 80) && strncmp(keyrec, endkey, 80))
            break;
    }

    if (ii == nspace && !tstatus)
    {
        /* check if the END keyword exists at the correct position */
        endpos=maxvalue( endpos, ( fptr->datastart - 2880 ) );
        ffmbyt(fptr, endpos, REPORT_EOF, &tstatus);  /* move to END position */
        ffgbyt(fptr, 80, keyrec, &tstatus); /* read the END keyword */
        if ( !strncmp(keyrec, endkey, 80) && !tstatus)
            return(*status);    /* END card was already correct */
    }

    /* header was not correctly terminated, so write the END and blank fill */
    ffmbyt(fptr, endpos, IGNORE_EOF, status); /* move to header end */
    for (ii=0; ii < nspace; ii++)
        ffpbyt(fptr, 80, blankkey, status);  /* write the blank keywords */

    /*
    The END keyword must either be placed immediately after the last
    keyword that was written (as indicated by the headend value), or
    must be in the first 80 bytes of the 2880-byte FITS record immediately 
    preceeding the data unit, whichever is further in the file. The
    latter will occur if space has been reserved for more header keywords
    which have not yet been written.
    */

    endpos=maxvalue( endpos, ( fptr->datastart - 2880 ) );
    ffmbyt(fptr, endpos, REPORT_EOF, status);  /* move to END position */

    ffpbyt(fptr, 80, endkey, status); /*  write the END keyword to header */

    if (*status > 0)
        ffpmsg("Error while writing END card (ffwend).");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpdfl(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)         /* IO - error status     */
/*
  Write the Data Unit Fill values if they are not already correct.
  The fill values are used to fill out the last 2880 byte block of the HDU.
  Fill the data unit with zeros or blanks depending on the type of HDU
  from the end of the data to the end of the current FITS 2880 byte block
*/
{
    char chfill, fill[2880];
    long fillstart;
    int nfill, tstatus, ii;

    if (*status > 0)
        return(*status);

    if (fptr->heapstart == 0)
        return(*status);      /* null data unit, so there is no fill */

    fillstart = fptr->datastart + fptr->heapstart + fptr->heapsize;

    nfill = (fillstart + 2879) / 2880 * 2880 - fillstart;

    if (fptr->hdutype == ASCII_TBL)
        chfill = 32;         /* ASCII tables are filled with spaces */
    else
        chfill = 0;          /* all other extensions are filled with zeros */

    tstatus = 0;

    if (!nfill)  /* no fill bytes; just check that entire table exists */
    {
        fillstart--;
        nfill = 1;
        ffmbyt(fptr, fillstart, REPORT_EOF, &tstatus); /* move to last byte */
        ffgbyt(fptr, nfill, fill, &tstatus);           /* get the last byte */

        if (tstatus == 0)
            return(*status);  /* no EOF error, so everything is OK */
    }
    else
    {
        ffmbyt(fptr, fillstart, REPORT_EOF, &tstatus); /* move to fill area */
        ffgbyt(fptr, nfill, fill, &tstatus);           /* get the fill bytes */

        if (tstatus == 0)
        {
            for (ii = 0; ii < nfill; ii++)
            {
                if (fill[ii] != chfill)
                    break;
            }

            if (ii == nfill)
                return(*status);   /* all the fill values were correct */
        }
    }

    /* fill values are incorrect or have not been written, so write them */

    memset(fill, chfill, nfill);  /* fill the buffer with the fill value */

    ffmbyt(fptr, fillstart, IGNORE_EOF, status); /* move to fill area */
    ffpbyt(fptr, nfill, fill, status); /* write the fill bytes */

    if (*status > 0)
        ffpmsg("Error writing Data Unit fill bytes (ffpdfl).");

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffcrhd(fitsfile *fptr,      /* I - FITS file pointer */
           int *status)         /* IO - error status     */
/*
  CReate Header Data unit:
  Create, initialize, and move the i/o pointer to a new extension
  after the last known existing extension of the FITS file.
*/
{
    long bytepos;

    if (*status > 0)
        return(*status);

    else if (fptr->maxhdu == MAXHDU)
        *status = BAD_HDU_NUM;       /* too many HDUs in file */

    else if (ffchdu(fptr, status) <= 0)  /* close the current HDU */
    {
      bytepos = fptr->headstart[ fptr->maxhdu + 1 ];  /* last known HDU */
      ffmbyt(fptr, bytepos, IGNORE_EOF, status);  /* move file ptr to it */
      fptr->maxhdu++;       /* increment the known number of HDUs */
      fptr->curhdu = fptr->maxhdu;      /* set current HDU location */
      fptr->nextkey = bytepos;          /* next keyword = start of header */
      fptr->headend = bytepos;          /* end of header */
      fptr->datastart = DATA_UNDEFINED; /* start of data unit undefined */
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffdblk(fitsfile *fptr,      /* I - FITS file pointer                    */
           long nblocks,        /* I - number of 2880-byte blocks to delete */
           int *status)         /* IO - error status                        */
/*
  Delete the specified number of 2880-byte blocks from the end
  of the CHDU by shifting all following extensions up this this
  number of blocks.
*/
{
    char buffer[2880];
    int tstatus, ii;
    long readpos, writepos;

    if (*status > 0)
        return(*status);

    tstatus = 0;

    /* pointers to the read and write positions */
    readpos = fptr->headstart[(fptr->curhdu) + 1];
    writepos = readpos - (nblocks * 2880);

    while ( !ffmbyt(fptr, readpos, REPORT_EOF, &tstatus) &&
            !ffgbyt(fptr, 2880L, buffer, &tstatus) )
    {
        ffmbyt(fptr, writepos, REPORT_EOF, status);
        ffpbyt(fptr, 2880L, buffer, status);

        if (*status > 0)
        {
           ffpmsg("Error deleting FITS blocks (ffdblk)");
           return(*status);
        }
        readpos  += 2880;  /* increment to next block to transfer */
        writepos += 2880;
    }

    /* now fill the last nblock blocks with zeros */
    memset(buffer, 0, 2880);

    ffmbyt(fptr, writepos, REPORT_EOF, status);

    for (ii = 0; ii < nblocks; ii++)
        ffpbyt(fptr, 2880L, buffer, status);

    /* recalculate the starting location of all subsequent HDUs */
    for (ii = fptr->curhdu; ii <= fptr->maxhdu; ii++)
         fptr->headstart[ii + 1] -= (2880 * nblocks);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffmahd(fitsfile *fptr,      /* I - FITS file pointer             */
           int hdunum,          /* I - number of the HDU to move to  */
           int *exttype,        /* O - type of extension, 0, 1, or 2 */
           int *status)         /* IO - error status                 */
/*
  Move to Absolute Header Data unit.  Move to the specified HDU
  and read the header to initialize the table structure.  Note that extnum 
  is one based, so the primary array is extnum = 1.
*/
{
    int moveto, tstatus;
    char message[FLEN_ERRMSG];

    if (*status > 0)
        return(*status);
    else if (hdunum < 1 || hdunum >= MAXHDU )
        return(*status = BAD_HDU_NUM);

    while( (fptr->curhdu) + 1 != hdunum) /* at the correct HDU? */
    {
        /* move directly to the extension if we know that it exists,
        otherwise move to the highest known extension.  */
        
        moveto = minvalue(hdunum - 1, (fptr->maxhdu) + 1);

        if (fptr->headstart[moveto] < fptr->logfilesize )  /* test if HDU exists */
        {
            if (ffchdu(fptr, status) <= 0)  /* close out the current HDU */
            {
                if (ffgext(fptr, moveto, exttype, status) > 0)
                {   /* failed to get the requested extension */
                    tstatus = 0;
                    ffrhdu(fptr, exttype, &tstatus); /* restore the CHDU */
                }
            }
        }
        else
            *status = END_OF_FILE;

        if (*status > 0)
        {
            sprintf(message,
            "Failed to move to HDU number %d (ffmahd).", hdunum);
            ffpmsg(message);
            return(*status);
        }
    }

    *exttype = fptr->hdutype; /* return the type of HDU */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffmrhd(fitsfile *fptr,      /* I - FITS file pointer                    */
           int hdumov,          /* I - rel. no. of HDUs to move by (+ or -) */ 
           int *exttype,        /* O - type of extension, 0, 1, or 2        */
           int *status)         /* IO - error status                        */
/*
  Move a Relative number of Header Data units.  Offset to the specified
  extension and read the header to initialize the HDU structure. 
*/
{
    int extnum;

    if (*status > 0)
        return(*status);

    extnum = fptr->curhdu + 1 + hdumov;  /* the absolute HDU number */
    ffmahd(fptr, extnum, exttype, status);  /* move to the HDU */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgext(fitsfile *fptr,      /* I - FITS file pointer                */
           int hdunum,          /* I - no. of HDU to move get (0 based) */ 
           int *exttype,        /* O - type of extension, 0, 1, or 2    */
           int *status)         /* IO - error status                    */
/*
  Get Extension.  Move to the specified extension and initialize the
  HDU structure.
*/
{
    int xcurhdu, xmaxhdu;
    long xheadend;

    if (*status > 0)
        return(*status);

    if (ffmbyt(fptr, fptr->headstart[hdunum], REPORT_EOF, status) <= 0)
    {
        /* temporarily save current values, in case of error */
        xcurhdu = fptr->curhdu;
        xmaxhdu = fptr->maxhdu;
        xheadend = fptr->headend;

        /* set new parameter values */
        fptr->curhdu = hdunum;
        fptr->maxhdu = maxvalue(fptr->maxhdu, hdunum);
        fptr->headend = 2000000000; /* set header end to huge value for now */

        if (ffrhdu(fptr, exttype, status) > 0)
        {   /* failed to get the new HDU, so restore previous values */
            fptr->curhdu = xcurhdu;
            fptr->maxhdu = xmaxhdu;
            fptr->headend = xheadend;
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffiblk(fitsfile *fptr,      /* I - FITS file pointer               */
           long nblock,         /* I - no. of blocks to insert         */ 
           int headdata,        /* O - insert where? 0=header, 1=data  */
           int *status)         /* IO - error status                   */
/*
   insert 2880-byte blocks at the end of the current header or data unit
*/
{
    int tstatus;
    long ii, jj, insertpt;
    char charfill;
    char buff1[2880], buff2[2880];
    char *inbuff, *outbuff, *tmpbuff;

    if (*status > 0)
        return(*status);

    tstatus = *status;

    if (headdata == 0 || fptr->hdutype == ASCII_TBL)
        charfill = 32;  /* headers and ASCII tables have space (32) fill */
    else
        charfill = 0;   /* images and binary tables have zero fill */

    for (jj = 0; jj < nblock; jj++)  /* insert one block at a time */
    {

      if (headdata == 0)  
        insertpt = fptr->datastart;  /* insert just before data, or */
      else
        insertpt = fptr->headstart[fptr->curhdu + 1]; /* before next HDU */


      inbuff = buff1;   /* set pointers to input and output buffers */
      outbuff = buff2;

      memset(outbuff, charfill, 2880); /* initialize buffer with fill */

      ffmbyt(fptr, insertpt, REPORT_EOF, status);  /* move to 1st point */
      ffgbyt(fptr, 2880, inbuff, status);  /* read first block of bytes */

      while (*status <= 0)
      {
        ffmbyt(fptr, insertpt, REPORT_EOF, status);  /* insert point */
        ffpbyt(fptr, 2880, outbuff, status);  /* write the output buffer */

        if (*status > 0)
            return(*status);

        tmpbuff = inbuff;   /* swap input and output pointers */
        inbuff = outbuff;
        outbuff = tmpbuff;
        insertpt += 2880;  /* increment insert point by 1 block */

        ffmbyt(fptr, insertpt, REPORT_EOF, status);  /* move to next block */
        ffgbyt(fptr, 2880, inbuff, status);  /* read block of bytes */
      }

      *status = tstatus;  /* reset status value */
      ffmbyt(fptr, insertpt, IGNORE_EOF, status); /* move back to insert pt */
      ffpbyt(fptr, 2880, outbuff, status);  /* write the final block */

      if (headdata == 0)
        fptr->datastart += 2880;  /* update data start address */

      for (ii = fptr->curhdu; ii <= fptr->maxhdu; ii++)
         fptr->headstart[ii + 1] += 2880; /* update following HDU addresses */
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffdtyp(char *cval,  /* I - formatted string representation of the value */
           char *dtype, /* O - datatype code: C, L, F or I */
          int *status)  /* IO - error status */
/*
  determine implicit datatype of input string.
  This assumes that the string conforms to the FITS standard
  for keyword values, so may not detect all invalid formats.
*/
{

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (cval[0] == '\'')
        *dtype = 'C';          /* character string starts with a quote */
    else if (cval[0] == 'T' || cval[0] == 'F')
        *dtype = 'L';          /* logical = T or F character */
    else if (strchr(cval,'.'))
        *dtype = 'F';          /* float contains a decimal point */
    else
        *dtype = 'I';          /* if none of the above must be an integer */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2x(char *cval,   /* I - formatted string representation of the value */
          char *dtype,  /* O - datatype code: C, L, F or I  */

    /* Only one of the following will be defined, depending on datatype */
          long *ival,    /* O - integer value       */
          int *lval,     /* O - logical value       */
          char *sval,    /* O - string value        */
          double *dval,  /* O - double value        */

          int *status)   /* IO - error status */
/*
  high level routine to convert formatted character string to its
  intrinsic data type
*/
{
    ffdtyp(cval, dtype, status);     /* determine the datatype */

    if (*dtype == 'I')
        ffc2ii(cval, ival, status);
    else if (*dtype == 'F')
        ffc2dd(cval, dval, status);
    else if (*dtype == 'L')
        ffc2ll(cval, lval, status);
    else 
        ffc2s(cval, sval, status);
        
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2i(char *cval,   /* I - string representation of the value */
          long *ival,   /* O - numerical value of the input string */
          int *status)  /* IO - error status */
/*
  convert formatted string to an integer value, doing implicit
  datatype conversion if necessary.
*/
{
    char dtype, sval[81], msg[81];
    int lval;
    double dval;
    
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* convert the keyword to its native datatype */
    ffc2x(cval, &dtype, ival, &lval, sval, &dval, status);

    if (dtype == 'C')
        *status = BAD_INTKEY;

    if (*status > 0)
    {
            *ival = 0;
            strcpy(msg,"Error in ffc2i evaluating string as an integer: ");
            strncat(msg,cval,30);
            ffpmsg(msg);
            return(*status);
    }

    if (dtype == 'F')
            *ival = (long) dval;
    else if (dtype == 'L')
            *ival = (long) lval;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2l(char *cval,  /* I - string representation of the value */
         int *lval,    /* O - numerical value of the input string */
         int *status)  /* IO - error status */
/*
  convert formatted string to a logical value, doing implicit
  datatype conversion if necessary
*/
{
    char dtype, sval[81], msg[81];
    long ival;
    double dval;
    
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* convert the keyword to its native datatype */
    ffc2x(cval, &dtype, &ival, lval, sval, &dval, status);

    if (dtype == 'C')
        *status = BAD_LOGICALKEY;

    if (*status > 0)
    {
            *lval = 0;
            strcpy(msg,"Error in ffc2l evaluating string as a logical: ");
            strncat(msg,cval,30);
            ffpmsg(msg);
            return(*status);
    }

    if (dtype == 'I')
    {
        if (ival)
            *lval = 1;
        else
            *lval = 0;
    }
    else if (dtype == 'F')
    {
        if (dval)
            *lval = 1;
        else
            *lval = 0;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2r(char *cval,   /* I - string representation of the value */
          float *fval,  /* O - numerical value of the input string */
          int *status)  /* IO - error status */
/*
  convert formatted string to a real float value, doing implicit
  datatype conversion if necessary
*/
{
    char dtype, sval[81], msg[81];
    long ival;
    int lval;
    double dval;
    
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* convert the keyword to its native datatype */
    ffc2x(cval, &dtype, &ival, &lval, sval, &dval, status);

    if (dtype == 'C')
        *status = BAD_FLOATKEY;

    if (*status > 0)
    {
            *fval = 0.;
            strcpy(msg,"Error in ffc2r evaluating string as a float: ");
            strncat(msg,cval,30);
            ffpmsg(msg);
            return(*status);
    }

    if (dtype == 'F')
        *fval = (float) dval;
    else if (dtype == 'I')
        *fval = (float) ival;
    else if (dtype == 'L')
        *fval = (float) lval;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2d(char *cval,   /* I - string representation of the value */
          double *dval, /* O - numerical value of the input string */
          int *status)  /* IO - error status */
/*
  convert formatted string to a double value, doing implicit
  datatype conversion if necessary
*/
{
    char dtype, sval[81], msg[81];
    long ival;
    int lval;
    
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /* convert the keyword to its native datatype */
    ffc2x(cval, &dtype, &ival, &lval, sval, dval, status);

    if (dtype == 'C')
        *status = BAD_DOUBLEKEY;

    if (*status > 0)
    {
            *dval = 0.;
            strcpy(msg,"Error in ffc2d evaluating string as a double: ");
            strncat(msg,cval,30);
            ffpmsg(msg);
            return(*status);
    }

    if (dtype == 'I')
        *dval = (double) ival;
    else if (dtype == 'L')
        *dval = (double) lval;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2ii(char *cval,  /* I - string representation of the value */
          long *ival,   /* O - numerical value of the input string */
          int *status)  /* IO - error status */
/*
  convert null-terminated formatted string to an integer value
*/
{
    char *loc;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    *ival = 0;
    *ival = strtol(cval, &loc, 10);  /* read the string as an integer */

    /* check for read error, or junk following the integer */
    if (*loc != '\0' && *loc != ' ' ) 
        *status = BAD_C2I;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2ll(char *cval,  /* I - string representation of the value: T or F */
           int *lval,   /* O - numerical value of the input string: 1 or 0 */
           int *status) /* IO - error status */
/*
  convert null-terminated formatted string to a logical value
*/
{
    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (cval[0] == 'T')
        *lval = 1;
    else                
        *lval = 0;        /* any character besides T is considered false */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2s(char *instr,  /* I - null terminated quoted input string */
          char *outstr, /* O - null terminated output string without quotes */
          int *status)  /* IO - error status */
/*
    convert an input quoted string to an unquoted string by removing
    the leading and trailing quote character.  Also, replace any
    pairs of single quote characters with just a single quote 
    character (FITS used a pair of single quotes to represent
    a literal quote character within the string).
*/
{
    int jj;
    size_t len, ii;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (instr[0] != '\'')
    {
        strcpy(outstr, instr);  /* no leading quote, so return input string */
        return(*status);
    }

    len = strlen(instr);

    for (ii=1, jj=0; ii < len; ii++, jj++)
    {
        if (instr[ii] == '\'')  /*  is this the closing quote?  */
        {
            if (instr[ii+1] == '\'')  /* 2 successive quotes? */
                ii++;  /* copy only one of the quotes */
            else
                break;   /*  found the closing quote, so exit this loop  */
        }
        outstr[jj] = instr[ii];   /* copy the next character to the output */
    }

    outstr[jj] = '\0';             /*  terminate the output string  */

    if (ii == len)
    {
        ffpmsg("This string value has no closing quote (ffc2s):");
        ffpmsg(instr);
        return(*status = 205);
    }

    for (jj--; jj >= 0; jj--)  /* replace trailing blanks with nulls */
    {
        if (outstr[jj] == ' ')
            outstr[jj] = 0;
        else
            break;
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2rr(char *cval,   /* I - string representation of the value */
           float *fval,  /* O - numerical value of the input string */
           int *status)  /* IO - error status */
/*
  convert null-terminated formatted string to a float value
*/
{
    char *loc, msg[81];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    *fval = 0.;
    *fval = (float) strtod(cval, &loc);  /* read the string as an float */

    /* check for read error, or junk following the value */
    if (*loc != '\0' && *loc != ' ' )
    {
        strcpy(msg,"Error in ffc2rr converting string to float: ");
        strncat(msg,cval,30);
        ffpmsg(msg);

        *status = BAD_C2F;   
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffc2dd(char *cval,   /* I - string representation of the value */
           double *dval, /* O - numerical value of the input string */
           int *status)  /* IO - error status */
/*
  convert null-terminated formatted string to a double value
*/
{
    char msg[81], tval[73], *loc;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);


    strcpy(tval, cval);
    loc = strchr(tval, 'D');

    if (loc)            /*  The C language does not support a 'D' */
       *loc = 'E';      /*  exponent so replace any D's with E's. */               
    *dval = 0.;
    *dval = strtod(tval, &loc);  /* read the string as an double */

    /* check for read error, or junk following the value */
    if (*loc != '\0' && *loc != ' ' )
    {
        strcpy(msg,"Error in ffc2dd converting string to double: ");
        strncat(msg,cval,30);
        ffpmsg(msg);

        *status = BAD_C2D;   
    }

    return(*status);
}

