#include "harpy/PetscLogger.h"

#include <petscsys.h>
#include <petsc/private/petscimpl.h>

#include <iostream>

/** 
 *
 * Encaminha a saida do Petsc para o logger
 *
 **/
PetscErrorCode mypetscvfprintf(FILE *fd, const char format[], va_list Argp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (fd != stdout && fd != stderr) { /* handle regular files */
    ierr = PetscVFPrintfDefault(fd, format, Argp);CHKERRQ(ierr);
  } else {
    char buff[1024]; /* Make sure to assign a large enough buffer */
    size_t  length;

    ierr = PetscVSNPrintf(buff, 1024, format, &length, Argp);CHKERRQ(ierr);

    /* now send buff to whatever stream or whatever you want */
    if ( fd == stdout ) std::cout << buff;
    if ( fd == stderr ) std::cerr << buff;

    // Log file
//    *(Log::_dfile) << buff;
  }
  PetscFunctionReturn(0);
}


/**
 *
 * Isso tem que rodar antes do init do LibMesh e do Petsc.
 *
 */
PetscLogger::PetscLogger() 
{
  PetscVFPrintf = mypetscvfprintf; 
}
