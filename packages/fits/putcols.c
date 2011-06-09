/*  This file, putcols.c, contains routines that write data elements to    */
/*  a FITS image or table, of type character string.                       */

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
int ffpcls( fitsfile *fptr,  /* I - FITS file pointer                       */
            int  colnum,     /* I - number of column to write (1 = 1st col) */
            long  firstrow,  /* I - first row to write (1 = 1st row)        */
            long  firstelem, /* I - first vector element to write (1 = 1st) */
            long  nelem,     /* I - number of strings to write              */
            char  **array,   /* I - array of pointers to strings            */
            int  *status)    /* IO - error status                           */
/*
  Write an array of string values to a column in the current FITS HDU.
*/
{
    int tcode, maxelem, hdutype, nchar;
    long twidth, incre, repeat, rowlen, rownum, elemnum, remain, next;
    long ii, jj, ntodo, tnull, startpos, wrtptr;
    double scale, zero;
    char tform[20], *blanks;
    char message[FLEN_ERRMSG];
    char snull[20];   /*  the FITS null value  */
    tcolumn *colptr;

    double cbuff[DBUFFSIZE / sizeof(double)]; /* align cbuff on word boundary */
    char *buffer, *arrayptr;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (fptr->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               /* rescan header */
            return(*status);

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if (colnum < 1 || colnum > fptr->tfield)
    {
        sprintf(message, "Specified column number is out of range: %d",
                colnum);
        ffpmsg(message);
        return(*status = BAD_COL_NUM);
    }

    colptr  = fptr->tableptr;   /* point to first column */
    colptr += (colnum - 1);     /* offset to correct column structure */
    tcode = colptr->tdatatype;

    if (tcode == -TSTRING) /* variable length column in a binary table? */
    {
      /* only write a single string; ignore value of firstelem */
      nchar = strlen(array[0]); 

      if (ffgcpr( fptr, colnum, firstrow, 1, nchar, 1, &scale, &zero,
        tform, &twidth, &tcode, &maxelem, &startpos,  &elemnum, &incre,
        &repeat, &rowlen, &hdutype, &tnull, snull, status) > 0)
        return(*status);

      remain = 1;
      twidth = nchar;
      blanks = 0;          /* initialize null pointer */
    }
    else if (tcode == TSTRING)
    {
      if (ffgcpr( fptr, colnum, firstrow, firstelem, nelem, 1, &scale, &zero,
        tform, &twidth, &tcode, &maxelem, &startpos,  &elemnum, &incre,
        &repeat, &rowlen, &hdutype, &tnull, snull, status) > 0)
        return(*status);

      blanks = (char *) malloc(twidth); /* string for blank fill values */
      if (!blanks)
      {
        ffpmsg("Could not allocate memory for string (ffpcls)");
        return(*status = ARRAY_TOO_BIG);
      }

      for (ii = 0; ii < twidth; ii++)
          blanks[ii] = ' ';          /* fill string with blanks */

      remain = nelem;           /* remaining number of values to write  */
    }
    else 
      return(*status = NOT_ASCII_COL);
 
    /*-------------------------------------------------------*/
    /*  Now write the strings to the FITS column.            */
    /*-------------------------------------------------------*/

    next = 0;                 /* next element in array to be written  */
    rownum = 0;               /* row number, relative to firstrow     */

    while (remain)
    {
      /* limit the number of pixels to process a one time to the number that
         will fit in the buffer space or to the number of pixels that remain
         in the current vector, which ever is smaller.
      */
      ntodo = minvalue(remain, maxelem);      
      ntodo = minvalue(ntodo, (repeat - elemnum));

      wrtptr = startpos + (rownum * rowlen) + (elemnum * incre);
      ffmbyt(fptr, wrtptr, IGNORE_EOF, status);  /* move to write position */

      buffer = (char *) cbuff;

      /* copy the user's strings into the buffer */
      for (ii = 0; ii < ntodo; ii++)
      {
         arrayptr = array[next];

         for (jj = 0; jj < twidth; jj++)  /*  copy the string, char by char */
         {
            if (*arrayptr)
            {
              *buffer = *arrayptr;
              buffer++;
              arrayptr++;
            }
            else
              break;
         }

         for (;jj < twidth; jj++)    /* fill field with blanks, if needed */
         {
           *buffer = ' ';
           buffer++;
         }

         next++;
      }

      /* write the buffer full of strings to the FITS file */
      if (incre == twidth)
         ffpbyt(fptr, ntodo * twidth, cbuff, status);
      else
         ffpbytoff(fptr, twidth, ntodo, incre - twidth, cbuff, status);

      if (*status > 0)  /* test for error during previous write operation */
      {
         sprintf(message,
          "Error writing elements %ld thru %ld of input data array (ffpcls).",
             next+1, next+ntodo);
         ffpmsg(message);

         if (blanks)
           free(blanks);

         return(*status);
      }

      /*--------------------------------------------*/
      /*  increment the counters for the next loop  */
      /*--------------------------------------------*/
      remain -= ntodo;
      if (remain)
      {
          elemnum += ntodo;
          if (elemnum == repeat)  /* completed a row; start on next row */
          {
              elemnum = 0;
              rownum++;
          }
       }
    }  /*  End of main while Loop  */

    if (blanks)
      free(blanks);

    return(*status);
}

