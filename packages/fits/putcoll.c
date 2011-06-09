/*  This file, putcoll.c, contains routines that write data elements to    */
/*  a FITS image or table, with logical datatype.                          */

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
#include "fitsio2.h"
/*--------------------------------------------------------------------------*/
int ffpcll( fitsfile *fptr,  /* I - FITS file pointer                       */
            int  colnum,     /* I - number of column to write (1 = 1st col) */
            long  firstrow,  /* I - first row to write (1 = 1st row)        */
            long  firstelem, /* I - first vector element to write (1 = 1st) */
            long  nelem,     /* I - number of values to write               */
            char *array,     /* I - array of values to write                */
            int  *status)    /* IO - error status                           */
/*
  Write an array of logical values to a column in the current FITS HDU.
*/
{
    int tcode, maxelem, hdutype;
    long twidth, incre, repeat, rowlen, rownum, elemnum, remain, next;
    long tnull, startpos, wrtptr;
    double scale, zero;
    char tform[20], ctrue = 'T', cfalse = 'F';
    char message[FLEN_ERRMSG];
    char snull[20];   /*  the FITS null value  */

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if (ffgcpr( fptr, colnum, firstrow, firstelem, nelem, 1, &scale, &zero,
        tform, &twidth, &tcode, &maxelem, &startpos,  &elemnum, &incre,
        &repeat, &rowlen, &hdutype, &tnull, snull, status) > 0)
        return(*status);

    if (tcode != TLOGICAL)   
        return(*status = NOT_LOGICAL_COL);

    /*---------------------------------------------------------------------*/
    /*  Now write the logical values one at a time to the FITS column.     */
    /*---------------------------------------------------------------------*/
    remain = nelem;           /* remaining number of values to write  */
    next = 0;                 /* next element in array to be written  */
    rownum = 0;               /* row number, relative to firstrow     */

    while (remain)
    {
      wrtptr = startpos + (rownum * rowlen) + (elemnum * incre);

      ffmbyt(fptr, wrtptr, IGNORE_EOF, status);  /* move to write position */

      if (array[next])
         ffpbyt(fptr, 1, &ctrue, status);
      else
         ffpbyt(fptr, 1, &cfalse, status);

      if (*status > 0)  /* test for error during previous write operation */
      {
        sprintf(message,
           "Error writing element %ld of input array of logicals (ffpcll).",
            next+1);
        ffpmsg(message);
        return(*status);
      }

      /*--------------------------------------------*/
      /*  increment the counters for the next loop  */
      /*--------------------------------------------*/
      remain--;
      if (remain)
      {
        next++;
        elemnum++;
        if (elemnum == repeat)  /* completed a row; start on next row */
        {
           elemnum = 0;
           rownum++;
        }
      }

    }  /*  End of main while Loop  */

    return(*status);
}
/*--------------------------------------------------------------------------*/
int ffpclx( fitsfile *fptr,  /* I - FITS file pointer                       */
            int  colnum,     /* I - number of column to write (1 = 1st col) */
            long  frow,      /* I - first row to write (1 = 1st row)        */
            long  fbit,      /* I - first bit to write (1 = 1st)            */
            long  nbit,      /* I - number of bits to write                 */
            char *larray,    /* I - array of logicals corresponding to bits */
            int  *status)    /* IO - error status                           */
/*
  write an array of logical values to a specified bit or byte
  column of the binary table.   If larray is TRUE, then the corresponding
  bit is set to 1, otherwise it is set to 0.
  The binary table column being written to must have datatype 'B' or 'X'. 
*/
{
    long bstart, offset, fbyte, bitloc, ndone;
    long ii, repeat, rstart, estart;
    int tcode, descrp;
    unsigned char cbuff;
    static unsigned char onbit[8] = {128,  64,  32,  16,   8,   4,   2,   1};
    static unsigned char offbit[8] = {127, 191, 223, 239, 247, 251, 253, 254};
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
        repeat = colptr->trepeat;

        if (tcode == TBIT)
            repeat = (repeat + 7) / 8; /* convert from bits to bytes */

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

        repeat = fbit + nbit -1;
        offset = fptr->heapsize;
        /* write the number of elements and the starting offset */
        ffpdes(fptr, colnum, frow, repeat, offset, status);
        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr->datastart + offset + fptr->heapstart + estart;
        /* increment the empty heap starting address (in bytes) */
        repeat = (repeat + 7) / 8;  /* convert from bits to bytes */
        fptr->heapsize += repeat;
    }

    /* move the i/o pointer to the start of the pixel sequence */
    ffmbyt(fptr, bstart, IGNORE_EOF, status);

    /* read the next byte (we may only be modifying some of the bits) */
    while (1)
    {
      if (ffgbyt(fptr, 1, &cbuff, status) == READ_ERROR)
      {
        /* hit end of file trying to read the byte, so just set byte = 0 */
        *status = 0;
        cbuff = 0;
      }

      /* move back, to be able to overwrite the byte */
      ffmbyt(fptr, bstart, REPORT_EOF, status);
 
      for (ii = bitloc; (ii < 8) && (ndone < nbit); ii++, ndone++)
      {
        if(larray[ndone])
          cbuff = cbuff | onbit[ii];
        else
          cbuff = cbuff & offbit[ii];
      }

      ffpbyt(fptr, 1, &cbuff, status); /* write the modified byte */

      if (ndone == nbit)  /* finished all the bits */
        return(*status);

      /* not done, so get the next byte */
      bstart++;
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

          ffmbyt(fptr, bstart, IGNORE_EOF, status);
        }
      }
      bitloc = 0;
    }
}

