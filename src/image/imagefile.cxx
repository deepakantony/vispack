#include <fstream>
#include <iostream>
#include "util/defs.h"
#include "image/imagefile.h"
#include "image/image.h"
#include "image/imageRGBA.h"
#ifdef USE_TIFF
//#include "tiffio.h"
#include "image/vispacktiffioCC.h"
#endif
#ifdef USE_FITS
// this is the fitsio C code
extern "C" 
{
#include "fitsio.h"
}
#endif
#ifdef USE_JPEG
extern "C"
{
#include "jerror.h"  
#include "jpeglib.h"
#include "setjmp.h"
struct my_error_mgr {
	    struct jpeg_error_mgr pub;
	    jmp_buf setjmp_buffer;
	};
typedef struct my_error_mgr * my_error_ptr;
METHODDEF(void)
    my_error_exit (j_common_ptr cinfo){
	my_error_ptr myerr = (my_error_ptr) cinfo->err;
	(*cinfo->err->output_message) (cinfo);
	longjmp(myerr->setjmp_buffer, 1);
    }
}
#endif


// #define DEBUG_PRINT(a) {printf("%s \n", a);}
#define DEBUG_PRINT(a) {}

// ---------------------------------------------------------------------------

int VISImageFile::write(const VISIm& image, const char* fname)
{
    int	err;
    if (!(image.isValid()))
	{
	    ERROR("VISImageFile::write got an invalid image\n");
	    err = -1;
	}
    else
    switch(image.type()) {
      case VISIm::OTHER:
	WARN("VISImageFile write - file type is OTHER");
	err = 0;
	break;
      case VISIm::BYTE:
#ifdef USE_FITS
	err = write_fits(VISImage<byte>(image), fname);
#else
	err = write_tiff(VISImage<byte>(image), fname );
#endif
	break;
      case VISIm::FLOAT:
#ifdef USE_FITS
	err = write_fits(VISImage<float>(image), fname);
#endif
	break;
      case VISIm::INT:
#ifdef USE_FITS
	err = write_fits(VISImage<int>(image), fname);
#else
	err = write_tiff(VISImage<int>(image), fname );
#endif
	break;
      case VISIm::RGBA:
#ifdef USE_TIFF
	err = write_tiff(VISImageRGBA(image), fname);
#endif
	break;
      case VISIm::NONE:
	WARN("VISImageFile write - file type is NONE");
	err = 0;
	break;
      default:
	WARN("VISImageFile write - unrecognized file type");
	err = 0;
	break;
    }
    return err;
}

// ---------------------------------------------------------------------------



int VISImageFile::write_tiff(const VISImageRGBA& image, const char* fname)
{
    unsigned i,j;
    unsigned pixels = image.width()*image.height();
    const rgba *pixRGBA = (image.rep())->buffer();
    byte* buf = new byte[3*pixels];
    int err;
  // must copy the data because tiff does not take const...
    for (i=0,j=0; i<pixels; i++) {
	buf[j++] = pixRGBA[i].r();
	buf[j++] = pixRGBA[i].g();
	buf[j++] = pixRGBA[i].b();
    }
#ifdef USE_TIFF
    // COMPRESSION_LZW
    err = WriteTIFFRGB(buf,image.width(),image.height(),
		       fname,COMPRESSION_LZW);
#endif
    delete buf;
    return err;
}

// ---------------------------------------------------------------------------

int VISImageFile::write_tiff(const VISImage<float>& image, const char* fname)
{
    const float *buf = (image.rep())->buffer();
    //    return(WriteTIFFfloat(buf, image.width(), image.height(), 1, fname,
    //			  COMPRESSION_LZW));
}

// ---------------------------------------------------------------------------

int VISImageFile::write_tiff(const VISImage<int>& image, const char* fname)
{
  int w = image.width(), h = image.height();
  int buf[w*h];
  memcpy((void*)buf, (const void*)((image.rep())->buffer()), w*h*sizeof(int));
  return(WriteTIFFint(buf, image.width(), image.height(), 1, fname,
		      COMPRESSION_LZW));
}

// ---------------------------------------------------------------------------

int VISImageFile::write_tiff(const VISImage<byte>& image, const char* fname)
{
  int w = image.width(), h = image.height();
  byte buf[w*h];
  memcpy((void*)buf, (const void*)((image.rep())->buffer()), w*h*sizeof(byte));
  return(WriteTIFFbyte(buf, image.width(), image.height(), 1, fname,
		       COMPRESSION_LZW));
}

#ifdef USE_FITS
// ---------------------------------------------------------------------------

int VISImageFile::write_fits(const VISImage<byte>& image, const char* fname)
{
    fitsfile *fptr;  // pointer to a FITS file 
    int status;
//    int ii, jj;
    long group, fpixel, nelements;
    byte *buff_byte, *buff_byte_ptr;
//    int index, k;
    // initialize FITS image parameters 
    int simple   = TRUE;
    int bitpix   =  8;   //       
    long naxis    =   2;  // 2-dimensional image                                
    long naxes[2];

    long pcount   =   0;  // no group parameters 
    long gcount   =   1;  // only a single image/group 
    int extend   = TRUE;  // there may be additional extension in the file  
    
    status = 0;         // initialize status before calling fitsio routines 

    if ( ffinit(&fptr, (char*)fname, &status) )    // create the new FITS file 
	{
	    printFitsError( status );   // call printerror if error occurs 
	    return(-1);	    
	}

// create the buffer
    naxes[0] = (long)image.width();
    naxes[1] = (long)image.height();
    nelements = naxes[0] * naxes[1];     // number of pixels to write 

    buff_byte_ptr = buff_byte = new byte[nelements];
    
    if ( ffphpr(fptr, simple, bitpix, naxis, // write the required keywords 
         naxes, pcount, gcount, extend, &status) )
	printFitsError( status );  // call printerror if error occurs 

    // load the data into the buffer
    for (int i = 0; i < image.height(); i++)
      for (int j = 0; j < image.width(); j++)
	{
	  *buff_byte_ptr++ = image.itemAt
	    (j, image.height() - (i + 1));
	}
    
    // write the array to the FITS file 
    group  = 0;                      // group to write            
    fpixel = 1;                          // first pixel to write      

    if ( ffpprb(fptr, group, fpixel, nelements, buff_byte, &status) )
        printFitsError( status );        // call printerror if error occurs 

    if ( ffclos(fptr, &status) )         // close the file 
	printFitsError( status );        // call printerror if error occurs 

    delete buff_byte;
    return(0);
}

// ---------------------------------------------------------------------------

int VISImageFile::write_fits(const VISImage<int>& image, const char* fname)
{
    fitsfile *fptr;  // pointer to a FITS file 
    int status;
//    int ii, jj;
    long group, fpixel, nelements;
    long *buff_long, *buff_long_ptr;
//    int index, k;
    // initialize FITS image parameters 
    int simple   = TRUE;
    int bitpix   =  32;   //       
    long naxis    =   2;  // 2-dimensional image                                
    long naxes[2];

    long pcount   =   0;  // no group parameters 
    long gcount   =   1;  // only a single image/group 
    int extend   = TRUE;  // there may be additional extension in the file  
    
    status = 0;         // initialize status before calling fitsio routines 

    if ( ffinit(&fptr, (char*)fname, &status) )    // create the new FITS file 
	{
	    printFitsError( status );   // call printerror if error occurs 
	    return(-1);	    
	}

// create the buffer
    naxes[0] = (long)image.width();
    naxes[1] = (long)image.height();
    nelements = naxes[0] * naxes[1];     // number of pixels to write 

    buff_long_ptr = buff_long = new long[nelements];
    

    if ( ffphpr(fptr, simple, bitpix, naxis, // write the required keywords 
         naxes, pcount, gcount, extend, &status) )
	printFitsError( status );  // call printerror if error occurs 


    // load the data into the buffer
    for (int i = 0; i < image.height(); i++)
	for (int j = 0; j < image.width(); j++)
	    {
		*buff_long_ptr++ = image.itemAt
		     (j, image.height() - (i + 1));
	    }
    
    // write the array to the FITS file 
    group  = 0;                      // group to write            
    fpixel = 1;                          // first pixel to write      

    if ( ffpprj(fptr, group, fpixel, nelements, buff_long, &status) )
        printFitsError( status );        // call printerror if error occurs 

    if ( ffclos(fptr, &status) )         // close the file 
	printFitsError( status );        // call printerror if error occurs 

    delete buff_long;
    return(0);
}

// ---------------------------------------------------------------------------

int VISImageFile::write_fits(const VISImage<float>& image, const char* fname)
{
  fitsfile *fptr;  // pointer to a FITS file 
  int status;
  //    int ii, jj;
  long group, fpixel, nelements;
  float *buff_float, *buff_float_ptr;
  //    int index, k;
  // initialize FITS image parameters 
  int simple   = TRUE;
  int bitpix   =  -32;   //       
  long naxis    =   2;  // 2-dimensional image                                
  long naxes[2];

  long pcount   =   0;  // no group parameters 
  long gcount   =   1;  // only a single image/group 
  int extend   = TRUE;  // there may be additional extension in the file  
    
  status = 0;         // initialize status before calling fitsio routines 

  if ( ffinit(&fptr, (char*)fname, &status) )    // create the new FITS file 
    {
      printFitsError( status ); // call printerror if error occurs 
      return(-1);	    
    }

  // create the buffer
  naxes[0] = image.width();
  naxes[1] = image.height();
  nelements = naxes[0] * naxes[1];     // number of pixels to write 

    
  if ( ffphpr(fptr, simple, bitpix, naxis, // write the required keywords 
	      naxes, pcount, gcount, extend, &status) )
    printFitsError( status );  // call printerror if error occurs 


  buff_float_ptr = buff_float = new float[nelements];

  // load the data into the buffer
  for (int i = 0; i < image.height(); i++)
    for (int j = 0; j < image.width(); j++)
      {
	*buff_float_ptr++ = image.itemAt
	  (j, image.height() - (i + 1));
      }
    
  // write the array to the FITS file 
  group  = 0;                      // group to write            
  fpixel = 1;                          // first pixel to write      

  if ( ffppre(fptr, group, fpixel, nelements, buff_float, &status) )
    printFitsError( status );        // call printerror if error occurs 

  if ( ffclos(fptr, &status) )         // close the file 
    printFitsError( status );        // call printerror if error occurs 

  delete buff_float;
  return(0);
}

// ---------------------------------------------------------------------------

int VISImageFile::write_fits(const VISImage<short>& image, const char* fname)
{
    fitsfile *fptr;  // pointer to a FITS file 
    int status;
//    int  ii, jj;
    long group, fpixel, nelements;
    short *buff_short, *buff_short_ptr;
//    int index, k;
    // initialize FITS image parameters 
    int simple   = TRUE;
    int bitpix   =  -32;   //       
    long naxis    =   2;  // 2-dimensional image                                
    long naxes[2];

    long pcount   =   0;  // no group parameters 
    long gcount   =   1;  // only a single image/group 
    int extend   = TRUE;  // there may be additional extension in the file  
    
    status = 0;         // initialize status before calling fitsio routines 

    if ( ffinit(&fptr, (char*)fname, &status) )    // create the new FITS file 
	{
	    printFitsError( status ); // call printerror if error occurs 
	    return(-1);	    
	}

// create the buffer
    naxes[0] = (long)image.width();
    naxes[1] = (long)image.height();
    nelements = naxes[0] * naxes[1];     // number of pixels to write 

    buff_short_ptr = buff_short = new short[nelements];

    if ( ffphpr(fptr, simple, bitpix, naxis, // write the required keywords 
         naxes, pcount, gcount, extend, &status) )
	printFitsError( status );  // call printerror if error occurs 

    // load the data into the buffer
    for (int i = 0; i < image.height(); i++)
	for (int j = 0; j < image.width(); j++)
	    {
		*buff_short_ptr++ = image.itemAt
		     (j, image.height() - (i + 1));
	    }
    
    // write the array to the FITS file 
    group  = 0;                      // group to write            
    fpixel = 1;                          // first pixel to write      

    if ( ffppri(fptr, group, fpixel, nelements, buff_short, &status) )
        printFitsError( status );        // call printerror if error occurs 

    if ( ffclos(fptr, &status) )         // close the file 
	printFitsError( status );        // call printerror if error occurs 

    delete buff_short;
    return(0);
}

#endif



// ---------------------------------------------------------------------------

int VISImageFile::write_tiff(const VISImage<short>& image, const char* fname)
{
    const short *buf = (image.rep())->buffer();
    //    return(WriteTIFFshort(buf, image.width(), image.height(), 1, fname,
    //			  COMPRESSION_LZW));
}

// ---------------------------------------------------------------------------

VISIm VISImageFile::read(const char* fname)
{
    int int_tmp;

#ifdef USE_FITS
    fitsfile *fptr;        // pointer to the FITS file
    int fits_status = 0;
#endif

#ifdef USE_TIFF
    //    PICINFO pinfo;
//    don't want to use no stinkin file system - Ross 3-22-96
//    if (!VISFileSystem::readable(fname)) 
//    {
//     WARN("file is not readable or does not exist!");
//	return nullImage;
//    }
//    else 

    FILE *file_tmp;

    VISIm nullIm = VISIm();

    if((file_tmp = fopen(fname,"r"))==NULL)
	{
	    printf("Imagefile::read - Unable to open %s\n",fname);
	    return(nullIm);
	}
    else
	fclose(file_tmp);

    // quiet the error handler for the open, in case it's just a bad magic number
    TIFFErrorHandler _my_tiff_error = NULL, _old_error_handler;
    _old_error_handler = TIFFSetErrorHandler(_my_tiff_error);
    TIFF* tif = TIFFOpen(fname, "r");
    // now turn back on the error handler
    TIFFSetErrorHandler(_old_error_handler);    
    if (tif) 
      {
      uint32 w, h;
      int i, j;
      size_t npixels;
      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
      //      short pixel_type;
      //      TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &pixel_type);
      VISImageRGBA im_color;

      //      if (pixel_type == PHOTOMETRIC_RGB)
	{
	  rgba buf_tmp;
	  im_color = VISImageRGBA(w, h);
	  uint32* raster;
	  
	  npixels = w * h;
	  //	  raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
	  raster = new uint32[npixels];
	  if (raster != NULL) 
	    {
	      if (TIFFReadRGBAImage(tif, w, h, raster, 0)) 
		{
		  for (j = (h-1); j >=0; j--)
		    for (i = 0; i < w; i++)
		      {
			//tiffs are in reverse order...
			buf_tmp.r(((byte*)raster)[0]);
			buf_tmp.g(((byte*)raster)[1]);
			buf_tmp.b(((byte*)raster)[2]);
			buf_tmp.a((byte)255);
			im_color.poke(i, j) = buf_tmp;
			raster++;
		      }
		  //		  _TIFFfree(raster);
		}
	      else 
		{
		  cout << "bad tiff read: " << fname << endl;
		  im_color = VISImageRGBA();
		}
	    }
	  else 
	    {
	      cout << "bad tiff read: " << fname << endl;
	      im_color = VISImageRGBA();
	    }
	  //	  delete raster;
	}
	//      else if (pixel_type == PHOTOMETRIC_RGB)
	//	{
	//	  
	//	}
      TIFFClose(tif);
      return im_color;
      }


// OLD TIFF STUFF - DOESN'T WORK WITH NEW TIFF LIBRARY.    
//
//     if (LoadTIFF(fname, &pinfo))
// 	{
// 	    DEBUG_PRINT("done tiff load ");
// 	    if (pinfo.colType == F_FULLCOLOR)
// 		{
// 		    DEBUG_PRINT("got full color ");
// 		    VISImageRGBA image_temp = VISImageRGBA(pinfo.w,pinfo.h);
// 		    rgba *buf_tmp = (image_temp.repRef())->bufferRef();
// 		    byte* byte_ptr = pinfo.pic;
// 		    for (int i = 0; i < (pinfo.w*pinfo.h); i++)
// 			{
// 			    buf_tmp->r(*(byte_ptr++));
// 			    buf_tmp->g(*(byte_ptr++));
// 			    (buf_tmp)->b(*(byte_ptr++));
// 			    (buf_tmp++)->a((byte)255);
// 			}
// 		    return image_temp;
// 		}
// 	    else if (pinfo.colType == F_COLORINDEXED)
// 		{
// 		    DEBUG_PRINT("got indexed color ");
// 		    VISImageRGBA image_temp = VISImageRGBA(pinfo.w,pinfo.h);
// 		    rgba *buf_tmp=(image_temp.repRef())->bufferRef();
// 		    byte* byte_ptr = pinfo.pic;
// 		    for (int i = 0; i < (pinfo.w*pinfo.h); i++)
// 			{
// 			    buf_tmp->r(pinfo.r[*(byte_ptr)]);
// 			    buf_tmp->g(pinfo.g[*(byte_ptr)]);
// 			    (buf_tmp)->b(pinfo.b[*(byte_ptr++)]);
// 			    (buf_tmp++)->a((byte)255);
// 			}
// 		    return image_temp;		
// 		}
// 	    else if (pinfo.colType == F_GREYSCALE)
// 		{
// 		    VISImage<byte> image_tmp;
// 		    DEBUG_PRINT("entered greyscale tiff load");
// 		    switch(pinfo.type)
// 			{
// 			  case PIC8:
// 			    DEBUG_PRINT("entered 8 bit tiff load");
// 			    image_tmp = VISImage<byte>
// 				(pinfo.w, pinfo.h, 1);
// 			    DEBUG_PRINT("created 8 bit image");
// //			    printf("height %d and width %d in read\n", 
// //				   image_tmp->width(), 
// //				   image_tmp->height());
// /* This is done for the time begin in order to prevent mixing malloc and
//  * new.  In the future, file io should proceed in two steps.  First, 
//  * get the header, and then read the data.  This gives the calling routine
//  * a chance to allocate the buffer */
// 			    (image_tmp.repRef())
// 				->copyBuffer((byte*)pinfo.pic);
// 			    free(pinfo.pic);
// 			    DEBUG_PRINT("copied 8 bit image");
// 			    return(image_tmp);
// 			  case PIC16:
// 			    DEBUG_PRINT("entered 16 bit tiff load");
// 			    short *buf_short;
// 			    buf_short = (short*)(pinfo.pic);
// 			    return(VISImage<short>(pinfo.w, pinfo.h, 1, 
// 							 &buf_short));
// 			  default:
// 			    WARN("can't read above 16 bit greyscale tiff images");
// 			    return nullImage;
// 			}
// 		}
// 	    else
// 		WARN("unrecognized file type");
// 	    return nullImage;
// 	}

	 
#endif // use_tiff
#ifdef USE_FITS
	    #define MAXDIM (10)
// try to read fits image
	    if (!(int_tmp 
		  = ffopen(&fptr, (char*)fname, READONLY, &fits_status)))
		{
//		    int nfound;
		    int anynull;
		    long naxes[MAXDIM], group, fpixel, nbuffer, npixels;
//		    long  ii;
//		    float datamin, datamax
		    float nullval;
		    float *buff_float;
		    byte *buff_byte;
		    long *buff_int;
		    int bitspix, simple, extend, naxis;
		    long pcount, gcount; 
	    
		    // read the NAXIS1 and NAXIS2 keyword to get image size 
		    if (ffghpr(fptr, MAXDIM, &simple, &bitspix, 
			       &naxis, naxes, &pcount, &gcount, 
			       &extend, &fits_status))
			ERROR("Imagefile: error reading fits info");
	    
// number of pixels in the image 
		    if (naxis > 3)
			{
			    ERROR("can read only 2d or 3d fits images");
			    return(nullIm);
			}

		    int buffsize;
		    if (naxis == 2)
			naxes[2] = 1;
		    buffsize = npixels  = naxes[0] * naxes[1] *
			naxes[2];
		    group    = 1;
		    fpixel   = 1;
		    nullval  = 0;         
		    //  don't check for null values in the image

//		    printf("bit per pix %d \n", bitspix);

		    switch(bitspix)
			{
			    //float
			  case -32:
			  case -64:
			    {
			      buff_float = new float[npixels];
				// read as many pixels as will fit in buffer
			      while (npixels > 0)
				{
				  nbuffer = npixels;
				  if (npixels > buffsize)
				    nbuffer = buffsize;

				  if ( ffgpve(fptr, group, fpixel, nbuffer, nullval,
					      buff_float, &anynull, &fits_status) )
				    ERROR("ERROR reading fits buffer");
				  // increment remaining number of pixels
				  npixels -= nbuffer;    
				  // next pixel to be read in image 
				  fpixel  += nbuffer;    
				}
			
			      float *buff_float_ptr = buff_float;
			      VISImage<float> r_float(naxes[0], naxes[1], 
						      naxes[2]);
			      for (int k = 0; k < r_float.channels(); k++)
				for (int i = 0; i < r_float.height(); i++)
				  for (int j = 0; j < r_float.width(); 
				       j++)
				    {
				      r_float.at(j, 
						 r_float.height() 
						 - (i + 1), k) 
					= *buff_float_ptr++;
				    }
			      delete buff_float;

				// close the file 
			      fits_status = 0;
			      if ( ffclos(fptr, &fits_status) )   
				// call printerror if error occurs 
				printFitsError( fits_status );

			      //			      VISImageFile im_file;
			      //			      im_file.write(r_float, "fits_read_float.fit");
			      return(r_float);
			    };
// not needed (unreachable)
// 			    break;
			  case 8:
			    {
				buff_byte = new byte[npixels];
// read as many pixels as will fit in buffer
				while (npixels > 0)
				    {
					nbuffer = npixels;
					if (npixels > buffsize)
					    nbuffer = buffsize;

					if ( ffgpvb(fptr, group, fpixel, 
						    nbuffer, (char)nullval,
						    buff_byte, &anynull, 
						    &fits_status) )
					    ERROR("ERROR reading fits buffer");
		    
// increment remaining number of pixels
					npixels -= nbuffer;    
// next pixel to be read in image 
					fpixel  += nbuffer;    
				    }
			
				byte *buff_byte_ptr = buff_byte;
				VISImage<byte> r_byte(naxes[0], naxes[1], 
						       naxes[2]);
				for (int k = 0; 
				     k < r_byte.channels(); k++)
				    for (int i = 0; 
					 i < r_byte.height(); i++)
					for (int j = 0; 
					     j < r_byte.width(); j++)
					    {
						r_byte.at(j, 
							  r_byte.height() 
							  - (i + 1), k) 
						    = *buff_byte_ptr++;
					    }
				delete buff_byte;
				// close the file 
				fits_status = 0;
				if ( ffclos(fptr, &fits_status) )   
				    // call printerror if error occurs 
				    printFitsError( fits_status );
				return(r_byte);
			}
// not needed (unreachable)
//			    break;
			  case 32:
			    {
				buff_int = new long[npixels];
// read as many pixels as will fit in buffer
				while (npixels > 0)
				    {
					nbuffer = npixels;
					if (npixels > buffsize)
					    nbuffer = buffsize;

					if ( ffgpvj(fptr, group, fpixel, 
						    nbuffer, (char)nullval,
						    buff_int, &anynull, 
						    &fits_status) )
					    ERROR("ERROR reading fits buffer");
		    
					npixels -= nbuffer;    
// next pixel to be read in image 
					fpixel  += nbuffer;    
				    }
			
				long *buff_int_ptr = buff_int;
				VISImage<int> r_int(naxes[0], naxes[1], 
						     naxes[2]);
				for (int k = 0; k < r_int.channels(); k++)
				    for (int i = 0; i < r_int.height(); i++)
					for (int j = 0; j < r_int.width(); j++)
					    {
						r_int.at(j, 
							 r_int.height() 
							 - (i + 1), k) 
						    = *buff_int_ptr++;
					    }
				delete buff_int;
				
				// close the file 
				fits_status = 0;
				if ( ffclos(fptr, &fits_status) )   
				    // call printerror if error occurs 
				    printFitsError( fits_status );
				return(r_int);
			    }
// not needed (unreachable)
//			    break;
			}
		    // close the FITS file
		    if ( ffclos(fptr, &fits_status) ) 
			ERROR("Could not close fits file");

		    return(nullIm);
	    
    
		}

	    
	       
#endif //use_fits
//--------------------------------------------------------------------------------
#ifdef USE_JPEG


    struct jpeg_decompress_struct cinfo;
    byte* image_buffer;
    struct my_error_mgr jerr;
    FILE * infile;		
    JSAMPARRAY buffer;	
    int row_stride,j,i;	
    if ((infile = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "Imagefile: read - Unable to open input file %s\n", fname);
	return (nullIm);
    }
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
  
    if (setjmp(jerr.setjmp_buffer)) {
	WARN("Image Read: File type not supported");
	jpeg_destroy_decompress(&cinfo);
	fclose(infile);
	return (nullIm);
    }

    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    (void) jpeg_read_header(&cinfo, TRUE);


    (void) jpeg_start_decompress(&cinfo);

    row_stride = cinfo.output_width * cinfo.output_components;

    buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    image_buffer=new byte[cinfo.output_width*cinfo.output_height*cinfo.output_components];

    i=0;
    while (cinfo.output_scanline < cinfo.output_height) {
       (void) jpeg_read_scanlines(&cinfo, buffer, 1);
       for (j=0;j<row_stride;j++) {
		image_buffer[i++]=buffer[0][j];
	    }
     }


if (row_stride==3*cinfo.output_width) {
    VISImageRGBA image_temp;
    image_temp=VISImageRGBA(cinfo.output_width, cinfo.output_height);
    rgba *buf_tmp = (image_temp.repRef())->bufferRef();
    for (int i = 0; i < (cinfo.output_width*cinfo.output_height); i++){
                            buf_tmp->r(*(image_buffer++));
                            buf_tmp->g(*(image_buffer++));                                                                   (buf_tmp)->b(*(image_buffer++));
                            (buf_tmp++)->a((byte)255);
                         }
 return image_temp;

}


else if (row_stride==cinfo.output_width) {
    
    VISImage<byte> image_tmp;
    image_tmp = VISImage<byte>(cinfo.output_width, cinfo.output_height, 1);
    (image_tmp.repRef())->copyBuffer(image_buffer);
    return(image_tmp);
}


    (void) jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);

 

#endif //use_jpeg
}

//----------------------------------------------------------------------
#ifdef USE_JPEG 


int VISImageFile::write_jpeg(const VISImageRGBA& image, const char* fname, int quality)
{
//Creating a buffer with the image data.

    unsigned i,j;
    unsigned col = image.width();
    unsigned row = image.height();

    const rgba *pixRGBA = (image.rep())->buffer();
    unsigned pixels=row*col;

    byte* outbuf = new byte[3*pixels];

    for (i=0,j=0; i<pixels; i++) {
	outbuf[j++] = pixRGBA[i].r();
	outbuf[j++] = pixRGBA[i].g();
	outbuf[j++] = pixRGBA[i].b();
    }

   
//Initialisation
    
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE * outfile;		
    JSAMPARRAY buffer;
    int row_stride,count;		
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if ((outfile = fopen(fname, "wb")) == NULL) {
    fprintf(stderr, "Imagefile:write - Unable to open output file %s\n",fname);
    return -1;
    }
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_width = col; 
    cinfo.image_height = row;
    cinfo.input_components = 3;		
    cinfo.in_color_space = JCS_RGB; 	
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);
    

    row_stride = cinfo.image_width * cinfo.input_components;
    
    buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    

//Writing on to file.
    j=0;
    while (cinfo.next_scanline < cinfo.image_height) {
	i=0;
	count=0;
	while(count<3*col){
	    buffer[0][i++]=outbuf[j++];count++;
	}
	(void) jpeg_write_scanlines(&cinfo,&buffer[0], 1);
	
    }

    delete outbuf;
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
    return 0;
}
    
#endif //USE-JPEG

//-----------------------------------------------------------------------------------------
//The 'quality' is a variable that controls the compression of the image 
//and is expressed as a percentage. 100% quality implies 0 compression.
#ifdef USE_JPEG

int VISImageFile::write_jpeg(const VISImage<float> image, const char* fname, int quality)
{
//Creating a buffer with the image data.

    unsigned i,j;
    unsigned col = image.width();
    unsigned row = image.height();
    const float *outbuf=(image.rep())->buffer();

//Initialisation
    
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE * outfile;		
    JSAMPARRAY buffer;
    int row_stride,count;		
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if ((outfile = fopen(fname, "wb")) == NULL) {
    fprintf(stderr, "Imagefile:write - Unable to open output file %s\n",fname);
    return -1;
    }
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_width = col; 
    cinfo.image_height = row;
    cinfo.input_components = 1;		
    cinfo.in_color_space = JCS_GRAYSCALE; 	
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);
    

    row_stride = cinfo.image_width * cinfo.input_components;
    
    buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    

//Writing on to file.
    j=0;
    while (cinfo.next_scanline < cinfo.image_height) {
	i=0;
	count=0;
	while(count<col){
	    buffer[0][i++]=(char)outbuf[j++];count++;
	}
	(void) jpeg_write_scanlines(&cinfo,&buffer[0], 1);
	
    }

//    delete [] outbuf;
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
    return 0;
}
    
#endif //USE-JPEG



 
//-----------------------------------------------------------------------
#ifdef USE_FITS

int VISImageFile::printFitsError(int fits_status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (fits_status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(fits_status, status_str);
        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", fits_status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }
    return(fits_status);
}

#endif


int VISImageFile::write_iv(const VISImage<float>& im, const char* filename)
{

    float x, y, z;
    int w = im.width();
    int h = im.height();
    int i, j;

    ofstream out_file;
    out_file.open(filename);
    if (!out_file) {
        cerr << filename << ": cannot open output file for write_iv.\n";
	return(-1);
    }
    
    out_file << "#Inventor V2.0 ascii\n\n"
        << "Separator {\n"
        << "     ShapeHints {\n"
        << "         vertexOrdering COUNTERCLOCKWISE }\n"
        << "     Coordinate3 {\n"
        << "          point [\n";

    // now actual data, omit border values 1 pixel wide, seems to screw things
    // up
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    	    x = (float)i;
			    y = (float)j;
			    z = im.itemAt(i, j);
			    out_file << "                "<<x<<" "
				     <<y<<" "<<-z<<",\n";
		}
	}

    out_file <<"              ]\n"
        <<"     }\n"
        <<"     QuadMesh {\n"
	     <<"          verticesPerColumn    " <<w - 2 << "\n"
	     <<"          verticesPerRow       " <<h - 2 << "\n"
	     <<"     }\n"
	     <<"}\n";
    
    // Close files and remove temporary files created
    out_file.close();
    return(0);
}


int VISImageFile::write_iv(const VISImage<float>& x, 
			    const VISImage<float>& y, 
			    const VISImage<float>& z, 
			    const char* filename)
{
    int w = x.width();
    int h = x.height();
    int i, j;
    float xx, yy, zz;

    ofstream out_file;
    out_file.open(filename);
    if (!out_file) {
        cerr << filename << ": cannot open output file for write_iv.\n";
	return(-1);
    }
    if (!(x.compareSize(y)&&x.compareSize(z)))
	{
        cerr << filename << "size mismatch in write_iv";
	return(-1);
	}
    
    out_file << "#Inventor V2.0 ascii\n\n"
        << "Separator {\n"
        << "     ShapeHints {\n"
        << "         vertexOrdering COUNTERCLOCKWISE }\n"
        << "     Coordinate3 {\n"
        << "          point [\n";

    // now actual data, omit border values 1 pixel wide, seems to screw things
    // up
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    	    xx = x.itemAt(i, j);
			    yy = y.itemAt(i, j);
			    zz = z.itemAt(i, j);
			    out_file << "                "<<xx<<" "
				     <<yy<<" "<<-zz<<",\n";
		}
	}

    out_file <<"              ]\n"
        <<"     }\n"
        <<"     QuadMesh {\n"
	     <<"          verticesPerColumn    " <<w - 2 << "\n"
	     <<"          verticesPerRow       " <<h - 2 << "\n"
	     <<"     }\n"
	     <<"}\n";
    
    // Close files and remove temporary files created
    out_file.close();
    return(0);
}


int VISImageFile::write_iv(const VISImage<float>& x, 
			    const VISImage<float>& y, 
			    const VISImage<float>& z, 
			   const VISImageRGBA& im_color, 
			    const char* filename)
{
    int w = x.width();
    int h = x.height();
    int i, j;

    ofstream out_file;
    out_file.open(filename);
    if (!out_file) {
        cerr << filename << ": cannot open output file for write_iv.\n";
	return(-1);
    }

    if (!(im_color.compareSize(x)&&im_color.compareSize(y)&&
	  im_color.compareSize(z)))
	{
	    cerr << "Size mismatch: write_iv.\n";
	    return(-1);
	}
    
    out_file << "#Inventor V2.0 ascii\n\n"
        << "Separator {\n"
        << "     ShapeHints {\n"
        << "         vertexOrdering COUNTERCLOCKWISE }\n"
        << "     Coordinate3 {\n"
        << "          point [\n";

    // now actual data, omit border values 1 pixel wide, seems to screw things
    // up
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		  out_file << "                "<< x.peek(i, j) <<" "
				     << y.peek(i, j) <<" "<< z.peek(i, j) << ",\n";
		}
	}

    out_file <<"              ]\n"
        <<"     }\n"

	<< "     Material {\n"
        << "          diffuseColor [\n";

    rgba color;
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    color = im_color.itemAt(i, j);
		    out_file << "                "<< (float)color.r()/255.0f
			     <<" " << (float)color.g()/255.0f 
			     <<" "<<(float)color.b()/255.0f<<",\n";
		}
	}
		    
	    out_file <<"              ]\n"
		     <<"     }\n";

		out_file << "MaterialBinding { \n"
		     <<	"         value PER_VERTEX_INDEXED \n" <<
		"        }" << endl;

		out_file <<"     QuadMesh {\n"
		     <<"          verticesPerColumn    " <<w - 2 << "\n"
	     <<"          verticesPerRow       " <<h - 2 << "\n"
	     <<"     }\n"
	     <<"}\n";
    
    // Close files and remove temporary files created
    out_file.close();
    return(0);
}

int VISImageFile::write_iv_trans(const VISImage<float>& x, 
			    const VISImage<float>& y, 
			    const VISImage<float>& z, 
			   const VISImage<float>& im_trans, 
			    const char* filename)
{
    int w = x.width();
    int h = x.height();
    int i, j;

    ofstream out_file;
    out_file.open(filename);
    if (!out_file) {
        cerr << filename << ": cannot open output file for write_iv.\n";
	return(-1);
    }

    if (!(im_trans.compareSize(x)&&im_trans.compareSize(y)&&
	  im_trans.compareSize(z)))
	{
	  cerr << "Size mismatch: write_iv.\n"; 
	  im_trans.print(); x.print(); y.print(); z.print(); 
	  return(-1);
	}
    
    out_file << "#Inventor V2.0 ascii\n\n"
        << "Separator {\n"
        << "     ShapeHints {\n"
        << "         vertexOrdering COUNTERCLOCKWISE }\n"
        << "     Coordinate3 {\n"
        << "          point [\n";

    // now actual data, omit border values 1 pixel wide, seems to screw things
    // up
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		  out_file << "                "<< x.peek(i, j) <<" "
				     << y.peek(i, j) <<" "<< z.peek(i, j) << ",\n";
		}
	}

    out_file <<"              ]\n"
        <<"     }\n"

	<< "     Material {\n"
        << "          transparency [\n";

    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    out_file << "                "<< im_trans.itemAt(i, j) <<",\n";
		}
	}
		    
	    out_file <<"              ]\n"
		     <<"     }\n";

		out_file << "MaterialBinding { \n"
		     <<	"         value PER_VERTEX_INDEXED \n" <<
		"        }" << endl;

		out_file <<"     QuadMesh {\n"
		     <<"          verticesPerColumn    " <<w - 2 << "\n"
	     <<"          verticesPerRow       " <<h - 2 << "\n"
	     <<"     }\n"
	     <<"}\n";
    
    // Close files and remove temporary files created
    out_file.close();
    return(0);
}



int VISImageFile::write_iv(const VISImage<float>& im, 
			    const VISImageRGBA& im_color, 
			    const char* filename)
{
    float x, y, z;
    int w = im.width();
    int h = im.height();
    int i, j;

    ofstream out_file;
    out_file.open(filename);
    if (!out_file) {
        cerr << filename << ": cannot open output file for write_iv.\n";
	return(-1);
    }

    if (!im_color.compareSize(im))
	{
	    cerr << "Size mismatch: write_iv.\n";
	    return(-1);
	}
    
    out_file << "#Inventor V2.0 ascii\n\n"
        << "Separator {\n"
        << "     ShapeHints {\n"
        << "         vertexOrdering COUNTERCLOCKWISE }\n"
        << "     Coordinate3 {\n"
        << "          point [\n";

    // now actual data, omit border values 1 pixel wide, seems to screw things
    // up
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    	    x = (float)i;
			    y = (float)j;
			    z = im.itemAt(i, j);
			    out_file << "                "<<x<<" "
				     <<y<<" "<<-z<<",\n";
		}
	}

    out_file <<"              ]\n"
        <<"     }\n"

	<< "     Material {\n"
        << "          diffuseColor [\n";

    rgba color;
    for (i = 1; i < w-1; i++) 
	{
	    for (j = 1; j < h-1; j++) 
		{
		    color = im_color.itemAt(i, j);
		    out_file << "                "<< (float)color.r()/255.0f
			     <<" " << (float)color.g()/255.0f 
			     <<" "<<(float)color.b()/255.0f<<",\n";
		}
	}
		    
	    out_file <<"              ]\n"
		     <<"     }\n";

		out_file << "MaterialBinding { \n"
		     <<	"         value PER_VERTEX_INDEXED \n" <<
		"        }" << endl;

		out_file <<"     QuadMesh {\n"
		     <<"          verticesPerColumn    " <<w - 2 << "\n"
	     <<"          verticesPerRow       " <<h - 2 << "\n"
	     <<"     }\n"
	     <<"}\n";
    
    // Close files and remove temporary files created
    out_file.close();
    return(0);
}







// ---------------------------------------------------------------------------


