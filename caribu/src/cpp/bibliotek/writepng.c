#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef WIN32
#include <alloc.h>
#endif

#include <png.h>
#include <cantools.h>

////////// / / / / / / / / Write_Png / / / / / / / / / ///////////
void Write_Png_File(FILE *out, 
	       int iViewHeight,
	       int iViewWidth ) {

  char *ImageBuffer ;
  int  iblig, ibcol, ibpix ;

  png_structp psWritePtr ;
  png_infop piWriteInfoPtr
    /*, piEndInfoPtr*/ ;
  png_bytep row_buf = (png_bytep)NULL;

#define SIZEOF_IEND 14
  /* piqué dans une image Gimp */
  char My_IEND[SIZEOF_IEND]={0x00,0x00,0x00,0x00,0x49,0x45,0x4e,
                   0x44,0xae, 0x42,0x60,0x82,0x0a,0x00 };
  

#define BYTE_PER_PIXEL 4
#define BIT_DEPTH 8
#define COLOR_TYPE 6
#define COMPRESS_TYPE 0
#define FILTER_TYPE 0

#ifdef DEBUG_ON
  png_debug(0, "Ecriving info struct\n");
  fprintf (stderr, "Entering Output (never mind...)\n");
#endif

  /* warnig: pbs potentiels de compatiblité GLint/int
   * compilation conditionelle ? */
  ImageBuffer = (char*) calloc(iViewWidth * iViewHeight * BYTE_PER_PIXEL ,
			       sizeof(char)) ; /* ou sizeof(char) */
  row_buf = (char*)calloc(iViewWidth * BYTE_PER_PIXEL, sizeof (png_byte));
  if ( ( ImageBuffer == NULL ) || ( row_buf == NULL )) {
    fprintf(stderr, "Plus de mémoire dans Output()\n");
    exit (Quitter( -1 )) ;
  }

  psWritePtr= png_create_write_struct( PNG_LIBPNG_VER_STRING,
    (png_voidp)NULL, (png_error_ptr)NULL, (png_error_ptr) NULL ) ;
  piWriteInfoPtr= png_create_info_struct (psWritePtr) ;
  png_init_io (psWritePtr, out ) ; 
  png_set_IHDR(psWritePtr, piWriteInfoPtr, iViewWidth, iViewHeight,
    BIT_DEPTH,  COLOR_TYPE, PNG_INTERLACE_NONE, COMPRESS_TYPE, FILTER_TYPE ) ;
  png_write_info(psWritePtr, piWriteInfoPtr ) ;

  // transfère le buffer en appliquant des modifs (coloration du fond...)
  Get_Image( ImageBuffer, iViewHeight, iViewWidth  );

  /* sortie Raw */
  for (iblig = iViewHeight-1 ; iblig >=0 ; iblig -- )
  {
    for( ibcol = 0; ibcol < iViewWidth; ibcol ++ ) {
      for (ibpix = 0; ibpix < BYTE_PER_PIXEL ; ibpix++ ) {
        row_buf[ibcol*BYTE_PER_PIXEL + ibpix] =
        ImageBuffer[iblig*iViewWidth*BYTE_PER_PIXEL+ibcol*BYTE_PER_PIXEL + ibpix ]*2 ;
      }
    }
    png_write_rows(psWritePtr, (png_bytepp) &row_buf, 1 ) ;
  }
  /* for*/
  /* ecrire l'info de fin IEND 
   * Version 0: piquer l'info ds une image et l'écrire brute */
  fwrite(My_IEND, SIZEOF_IEND, 1, out ) ;

  free (ImageBuffer) ;
}
/* Write_Png() */
