#include <stdlib.h>
#include <tiffio.h>
#include <image.h>

int WriteTIFFRGB(byte *data, int w, int h, const char *fname, int comp)
{
  TIFF *tif;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) return(-1);

  tif = TIFFFdOpen(dup(fileno(fp)), fname, "w");
  if (!tif) return(-1);
  
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, comp);
  if (comp == COMPRESSION_CCITTFAX3)
      TIFFSetField(tif, TIFFTAG_GROUP3OPTIONS,
		   GROUP3OPT_2DENCODING+GROUP3OPT_FILLBITS);

  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);


  /* write the image data */
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFWriteEncodedStrip(tif, 0, data, 3*w*h);

  TIFFClose(tif);

  if (fclose(fp) == EOF) {
      printf("tiffwriteRGB: somethings wrong unexpected end of file\n");
  }

  return 0;
}

int WriteTIFFint(int *data, int w, int h, int channels, const char *fname, int comp)
{
  TIFF *tif;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) return(-1);

    tif = TIFFFdOpen(dup(fileno(fp)), fname, "w");
    if (!tif) return(-1);
  
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, comp);
    if (comp == COMPRESSION_CCITTFAX3)
	TIFFSetField(tif, TIFFTAG_GROUP3OPTIONS,
		     GROUP3OPT_2DENCODING+GROUP3OPT_FILLBITS);

    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);

    /* write the image data */
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, sizeof(int)*8);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, channels);
    TIFFWriteEncodedStrip(tif, 0, (unsigned char *)data,  sizeof(int)*channels*w*h);

    TIFFClose(tif);

    if (fclose(fp) == EOF) {
	printf("tiffwriteRGB: somethings wrong unexpected end of file\n");
    }

    return 0;
}

int WriteTIFFbyte(byte *data, int w, int h, int channels, const char *fname, int comp)
{
  TIFF *tif;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) return(-1);

    tif = TIFFFdOpen(dup(fileno(fp)), fname, "w");
    if (!tif) return(-1);
  
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, comp);
    if (comp == COMPRESSION_CCITTFAX3)
	TIFFSetField(tif, TIFFTAG_GROUP3OPTIONS,
		     GROUP3OPT_2DENCODING+GROUP3OPT_FILLBITS);

    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);

    /* write the image data */
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, sizeof(byte)*8);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, channels);
    TIFFWriteEncodedStrip(tif, 0, (unsigned char *)data,  sizeof(byte)*channels*w*h);

    TIFFClose(tif);

    if (fclose(fp) == EOF) {
	printf("tiffwriteRGB: somethings wrong unexpected end of file\n");
    }

    return 0;
}





