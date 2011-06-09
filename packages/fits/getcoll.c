/*  This file, getcoll.c, contains routines that read data elements from   */
/*  a FITS image or table, with logical datatype.                          */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.  Users shall not, without prior written   */
/*  permission of the U.S. Government,  establish a claim to statutory     */
/*  copyright.  The Government and others acting on its behalf, shall have */
/*  a royalty-free, non-exclusive, irrevocable,  worldwide license for     */
/*  Government purposes to publish, distribute, translate, copy, exhibit,  */
/*  and perform such material.                                             */

#include <stdlib.h>
#include <string.h>
#include "fitsio2.h"
/*--------------------------------------------------------------------------*/
int ffgcl(  fitsfile *fptr,   /* I - FITS file pointer                       */
            int  colnum,      /* I - number of column to read (1 = 1st col)  */
            long  firstrow,   /* I - first row to read (1 = 1st row)         */
            long  firstelem,  /* I - first vector element to read (1 = 1st)  */
            long  nelem,      /* I - number of values to read                */
            char *array,      /* O - array of values                         */
            int  *status)     /* IO - error status                           */
/*
  Read an array of logical values from a column in the current FITS HDU.
  No checking for null values will be performed.
*/
{
    char cdummy[2];
    int idummy;

    cdummy[0] = 127;   /* special flag value to signal no null detection */
    idummy = 127;

    ffgcfl( fptr, colnum, firstrow, firstelem, nelem, array,
            cdummy, &idummy, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcfl( fitsfile *fptr,   /* I - FITS file pointer                       */
            int  colnum,      /* I - number of column to read (1 = 1st col)  */
            long  firstrow,   /* I - first row to read (1 = 1st row)         */
            long  firstelem,  /* I - first vector element to read (1 = 1st)  */
            long  nelem,      /* I - number of values to read                */
            char *array,      /* O - array of values                         */
            char *nularray,   /* O - array of flags = 1 if nultyp = 2        */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of logical values from a column in the current FITS HDU.
*/
{
    int tcode, maxelem, hdutype, ii, nulcheck;
    long twidth, offset, incre, repeat, rowlen, rownum, elemnum;
    long tnull, startpos, readptr, remain, next, ntodo;
    double scale, zero;
    char tform[20];
    char message[FLEN_ERRMSG];
    char snull[20];   /*  the FITS null value  */
    unsigned char buffer[DBUFFSIZE];

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (*anynul == 127 && *nularray == 127)
        nulcheck = 0;      /* do not test for null values */

    else
    {
        nulcheck = 1;
        *anynul = 0;
        memset(nularray, 0, nelem);   /* initialize nullarray */
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if (ffgcpr( fptr, colnum, firstrow, firstelem, nelem, 0, &scale, &zero,
        tform, &twidth, &tcode, &maxelem, &startpos,  &elemnum, &incre,
        &repeat, &rowlen, &hdutype, &tnull, snull, status) > 0)
        return(*status);

    if (tcode != TLOGICAL)   
        return(*status = NOT_LOGICAL_COL);
 
    /*---------------------------------------------------------------------*/
    /*  Now read the logical values from the FITS column.                  */
    /*---------------------------------------------------------------------*/

    remain = nelem;           /* remaining number of values to read */
    next = 0;                 /* next element in array to be read   */
    rownum = 0;               /* row number, relative to firstrow   */
    ntodo = remain;           /* max number of elements to read at one time */

    while (ntodo)
    {
      /*
         limit the number of pixels to read at one time to the number that
         remain in the current vector.    
      */
      ntodo = minvalue(ntodo, maxelem);      
      ntodo = minvalue(ntodo, (repeat - elemnum));

      readptr = startpos + (rownum * rowlen) + (elemnum * incre);

      ffgi1b(fptr, readptr, ntodo, incre, buffer, status);

      for (ii = 0; ii < ntodo; ii++, next++) /* convert from T or F to 1 or 0 */
      {
        if (buffer[ii] == 'T')
          array[next] = 1;
        else if (buffer[ii] =='F') 
          array[next] = 0;
        else if (nulcheck)
        {
          nularray[next] = 1;    /* this is an undefined pixel value */
          *anynul = 1;
        }
      }

      if (*status > 0)  /* test for error during previous read operation */
      {
        sprintf(message,
          "Error reading elements %ld thruough %ld of logical array (ffgcl).",
           next+1, next + ntodo);
        ffpmsg(message);
        return(*status);
      }

      /*--------------------------------------------*/
      /*  increment the counters for the next loop  */
      /*--------------------------------------------*/
      remain -= ntodo;
      if (remain)
      {
        next += ntodo;
        elemnum += ntodo;

        if (elemnum == repeat)  /* completed a row; start on later row */
          {
            elemnum = 0;
            rownum++;
          }
      }
      ntodo = remain;  /* this is the maximum number to do in next loop */

    }  /*  End of main while Loop  */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffgcx(  fitsfile *fptr,  /* I - FITS file pointer                       */
            int   colnum,    /* I - number of column to write (1 = 1st col) */
            long  frow,      /* I - first row to write (1 = 1st row)        */
            long  fbit,      /* I - first bit to write (1 = 1st)            */
            long  nbit,      /* I - number of bits to write                 */
            char *larray,    /* O - array of logicals corresponding to bits */
            int  *status)    /* IO - error status                           */
/*
  read an array of logical values from a specified bit or byte
  column of the binary table.   If larray is TRUE, then the corresponding
  bit is set to 1, otherwise it is set to 0.
  The binary table column being read to must have datatype 'B' or 'X'. 
*/
{
    long bstart, offset, fbyte, bitloc, ndone;
    long ii, repeat, rstart, estart;
    int tcode, descrp;
    unsigned char cbuff;
    static unsigned char onbit[8] = {128,  64,  32,  16,   8,   4,   2,   1};
    tcolumn *colptr;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /*  check input parameters */
    if (nbit < 1)
        return(*status);
    else if (frow < 1)
        return(*status = BAD_ROW_NUM);
    else if (fbit < 1)
        return(*status = BAD_ELEM_NUM);

    fbyte = (fbit + 7) / 8;
    bitloc = fbit - 1 - ((fbit - 1) / 8 * 8);
    ndone = 0;
    rstart = frow - 1;
    estart = fbyte - 1;

    colptr  = fptr->tableptr;   /* point to first column */
    colptr += (colnum - 1);     /* offset to correct column structure */

    tcode = colptr->tdatatype;

    if (abs(tcode) > TBYTE)
        return(*status = NOT_LOGICAL_COL); /* not correct datatype column */

    if (tcode > 0)
    {
        descrp = FALSE;  /* not a variable length descriptor column */
        /* N.B: REPEAT is the number of bytes, not number of bits */
        repeat = colptr->trepeat;

        if (tcode == TBIT)
            repeat = (repeat + 7) / 8;  /* convert from bits to bytes */

        if (fbyte > repeat)
            return(*status = BAD_ELEM_NUM);

        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr->datastart + (rstart * fptr->rowlength) +
               colptr->tbcol + estart;
    }
    else
    {
        descrp = TRUE;  /* a variable length descriptor column */
        /* only bit arrays (tform = 'X') are supported for variable */
        /* length arrays.  REPEAT is the number of BITS in the array. */

        ffgdes(fptr, colnum, frow, &repeat, &offset, status);
        repeat = (repeat + 7) / 8;

        if ((fbit + nbit + 6) / 8 > repeat)
            return(*status = BAD_ELEM_NUM);

        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr->datastart + offset + fptr->heapstart + estart;
    }

    /* move the i/o pointer to the start of the pixel sequence */
    if (ffmbyt(fptr, bstart, REPORT_EOF, status) > 0)
        return(*status);

    /* read the next byte */
    while (1)
    {
      if (ffgbyt(fptr, 1, &cbuff, status) > 0)
        return(*status);

      for (ii = bitloc; (ii < 8) && (ndone < nbit); ii++, ndone++)
      {
        if(cbuff & onbit[ii])       /* test if bit is set */
          larray[ndone] = TRUE;
        else
          larray[ndone] = FALSE;
      }

      if (ndone == nbit)   /* finished all the bits */
        return(*status);

      /* not done, so get the next byte */
      if (!descrp)
      {
        estart++;
        if (estart == repeat) 
        {
          /* move the i/o pointer to the next row of pixels */
          estart = 0;
          rstart = rstart + 1;
          bstart = fptr->datastart + (rstart * fptr->rowlength) +
               colptr->tbcol;

          ffmbyt(fptr, bstart, REPORT_EOF, status);
        }
      }
      bitloc = 0;
    }
}

