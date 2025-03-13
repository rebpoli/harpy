#ifndef __TERM_COLORS_HXX
#define __TERM_COLORS_HXX

/* FOREGROUND */
#define RST  "\x1B[0m"
#define KFGBLK  "\x1B[30m"
#define KFGRED  "\x1B[31m"
#define KFGGRN  "\x1B[32m"
#define KFGYEL  "\x1B[33m"
#define KFGBLU  "\x1B[34m"
#define KFGMAG  "\x1B[35m"
#define KFGCYN  "\x1B[36m"
#define KFGWHT  "\x1B[37m"

/* BACKGROUND */
#define KBGBLK  "\x1B[40m"
#define KBGRED  "\x1B[41m"
#define KBGGRN  "\x1B[42m"
#define KBGYEL  "\x1B[43m"
#define KBGBLU  "\x1B[44m"
#define KBGMAG  "\x1B[45m"
#define KBGCYN  "\x1B[46m"
#define KBGWHT  "\x1B[47m"

// FOREGROUND
#define FFBLK(x) KFGBLK x RST
#define FFRED(x) KFGRED x RST
#define FFGRN(x) KFGGRN x RST
#define FFYEL(x) KFGYEL x RST
#define FFBLU(x) KFGBLU x RST
#define FFMAG(x) KFGMAG x RST
#define FFCYN(x) KFGCYN x RST
#define FFWHT(x) KFGWHT x RST
// BACKGROUND
#define FBBLK(x) KBGBLK x RST
#define FBRED(x) KBGRED x RST
#define FBGRN(x) KBGGRN x RST
#define FBYEL(x) KBGYEL x RST
#define FBBLU(x) KBGBLU x RST
#define FBMAG(x) KBGMAG x RST
#define FBCYN(x) KBGCYN x RST
#define FBWHT(x) KBGWHT x RST


#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

#endif 
