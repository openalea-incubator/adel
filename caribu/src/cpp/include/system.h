#ifndef _SYSTEM_H
#define _SYSTEM_H

#define CANOPY_SHMID 68 // must be below 100 (-1<id<100)

#define BVIS_EXT ".bvis"
#define PNG_EXT  ".png"

#ifdef WIN32
#define PATH_SEP '\\'
#define MAX_PATH_LEN 512
#define DEVNUL "C:\\nul"
#define COPY   "copy"
#define KILL	"kill"
#define RM "del"
#define MKDIR "md"

#define S2V    "s2v10.exe"
#define MCSAIL "mcsail19.exe"
#define CHRONO ""
#define CRASH_TARGET	"cpfg.exe"
#ifdef IRRADIANCE
#define RADEXEC "canestrad44.exe"
#else
#define RADEXEC "canestrad.exe"
#endif
#define AWK	"my_awk.exe"
#define AWK_QUOTE	'"'

#else

#define PATH_SEP '/'
#define MAX_PATH_LEN 1024
#define DEVNUL "/dev/null"
#define COPY   "cp"
#define KILL	"kill.exe"
#define RM "rm -f"
#define TEMP "/tmp/"
#define MKDIR "mkdir"

#define S2V    "s2v"
#define MCSAIL "mcsail-1.9"
#define CHRONO "/usr/bin/time"
#define RADEXEC "canestrad"
#define AWK	"awk"
#define AWK_QUOTE	'\''


#endif

// echanges caribu <==> canestra
#define DG_NAME  "diag_"
#define NZ_NAME  "nz_"
#define BF_NAME  "Bfar_"


// 1 standard compatible CD iso-9660
#define NAME_LENGTH 130

// renvoie TMP/nom sous DOS
// TMP est recherché à l'exécution
char *GetAllFileName(char *nom) ;

#endif
