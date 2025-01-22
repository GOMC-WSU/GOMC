/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   dcdlib contains C routines for reading and writing binary DCD
   files.  The output format of these files is based on binary FORTRAN
   output, so its pretty ugly.  If you are squeamish, don't look!
*/

#include "DCDlib.h"

#ifndef OUTPUT_SINGLE_FILE
#define OUTPUT_SINGLE_FILE 1
// including Output.h causes compile error on BlueGeneQ
#endif

#ifdef WIN32
#define access(PATH, MODE) _access(PATH, 00)
#endif

#ifndef O_LARGEFILE
#define O_LARGEFILE 0x0
#endif

// same as write, only does error checking internally
void NAMD_write(int fd, const char *buf, size_t count, const char *errmsg) {
  int bad_cntr = 0;
  while (count) {
#if defined(WIN32) && !defined(__CYGWIN__)
    long retval = _write(fd, buf, count);
#else
    ssize_t retval = write(fd, buf, count);
#endif
    if (retval < 0 && errno == EINTR)
      retval = 0;
    if (retval < 0) {
      // include the error description in the message
      std::string s(errmsg);
      s += ": ";
      s += strerror(retval);
      // handle transient errors by retrying a few times
      if (bad_cntr < 5) {
        NAMD_warn(s.c_str());
        bad_cntr++;
      } else
        NAMD_err(s.c_str());
    }
    if (retval > (long)count)
      NAMD_bug("extra bytes written in NAMD_write()");
    buf += retval;
    count -= retval;
  }
}

// same as open, only does error checking internally
int NAMD_open(const char *fname) {
  int fd;

  //  open the file and die if the open fails
#ifdef WIN32
  while (
      (fd = _open(fname, O_WRONLY | O_CREAT | O_EXCL | O_BINARY | O_LARGEFILE,
                  _S_IREAD | _S_IWRITE)) < 0) {
#else
#ifdef NAMD_NO_O_EXCL
  while ((fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
#else
  while ((fd = open(fname, O_WRONLY | O_CREAT | O_EXCL | O_LARGEFILE,
#endif
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) < 0) {
#endif
    if (errno != EINTR) {
      char errmsg[1024];
      sprintf(errmsg, "Unable to open binary file %s", fname);
      NAMD_err(errmsg);
    }
  }

  return fd;
}

// same as close, only does error checking internally
void NAMD_close(int fd, const char *fname) {
#ifdef WIN32
  while (_close(fd)) {
#else
  while (close(fd)) {
#endif
    if (errno != EINTR) {
      char errmsg[1024];
      sprintf(errmsg, "Error on closing file %s", fname);
      NAMD_err(errmsg);
    }
  }
}

// move filename to filename.BAK
void NAMD_backup_file(const char *filename, const char *extension) {
  if (NAMD_file_exists(filename)) {
    if (!extension)
      extension = ".BAK";
    char *backup = new char[strlen(filename) + strlen(extension) + 1];
    strcpy(backup, filename);
    strcat(backup, extension);
#if defined(WIN32) && !defined(__CYGWIN__)
    if (remove(backup))
      if (errno != ENOENT) {
        const char *sys_err_msg = strerror(errno);
        if (!sys_err_msg)
          sys_err_msg = "(unknown error)";
        std::cout << "Error: "
                  << "Error on removing file " << backup << ": " << sys_err_msg
                  << std::endl;
        fflush(stdout);
      }
#endif
    while (rename(filename, backup)) {
      if (errno == EINTR)
        continue;
      const char *sys_err_msg = strerror(errno);
      if (!sys_err_msg)
        sys_err_msg = "(unknown error)";
      std::cout << "Error: "
                << "Error on renaming file " << filename << " to " << backup
                << ": " << sys_err_msg << std::endl;
      fflush(stdout);
      if (errno == EXDEV)
        continue;
      break;
    }
    delete[] backup;
  }
}

int NAMD_file_exists(const char *filename) {
  int rval;
  do {
    rval = access(filename, F_OK);
  } while (rval != 0 && errno == EINTR);
  if (rval != 0 && errno != ENOENT) {
    const char *sys_err_msg = strerror(errno);
    if (!sys_err_msg)
      sys_err_msg = "(unknown error)";
    std::cout << "Error: "
              << "Error on checking file " << filename << ": " << sys_err_msg
              << std::endl;
    fflush(stdout);
  }
  return !rval;
}

void NAMD_die(const char *err_msg) {
  if (!err_msg)
    err_msg = "(unknown error)";
  printf("ERROR: %s\n", err_msg);
  fflush(stdout);
  exit(EXIT_FAILURE);
}

void NAMD_err(const char *err_msg) {
  if (!err_msg)
    err_msg = "(unknown error)";
  printf("ERROR: %s\n", err_msg);
  fflush(stdout);
  exit(EXIT_FAILURE);
}

void NAMD_warn(const char *err_msg) {
  if (!err_msg)
    err_msg = "(unknown warning)";
  printf("Warning: %s\n", err_msg);
  fflush(stdout);
}

void NAMD_bug(const char *err_msg) {
  if (!err_msg)
    err_msg = "(unknown error)";
  const char *bug_msg =
      "ERROR: Report to https://github.com/GOMC-WSU/GOMC/issues";
  printf("FATAL ERROR: %s\n%s\n", err_msg, bug_msg);
  fflush(stdout);
  exit(EXIT_FAILURE);
}

#ifdef WIN32
#define LSEEK _lseeki64
#define READ _read
#else
#define LSEEK lseek
#define READ read
#endif

OFF_T NAMD_seek(int file, OFF_T offset, int whence) {
  OFF_T retval = LSEEK(file, offset, whence);
  if (retval < 0)
    NAMD_err("seek failed while writing DCD file");
  if (whence == SEEK_SET && retval != offset) {
    char buf[256];
    sprintf(buf,
            "seek failed while writing DCD file: SEEK_SET %ld returned %ld\n",
            offset, retval);
    NAMD_die(buf);
  }
  return retval;
}

#undef LSEEK
#define LSEEK NAMD_seek

/************************************************************************/
/*                                    */
/*                MACRO CHECK_FREAD            */
/*                                    */
/*    CHECK_FREAD is used to check the status of a file read.  It     */
/*   is passed the return code from read and a string to print out if   */
/*   an error is detected.  If an error is found, an error message is   */
/*   printed out and the program terminates.  This was made into a      */
/*   macro because it had to be done over and over and over . . .    */
/*                                    */
/************************************************************************/

#define CHECK_FREAD(X, msg)                                                    \
  if (X == -1) {                                                               \
    return (DCD_BADREAD);                                                      \
  }

/************************************************************************/
/*                                    */
/*            MACRO CHECK_FEOF                */
/*                                    */
/*    CHECK_FEOF checks for an abnormal EOF on a file read.  It is    */
/*   passed the return status from read and a message to print on error.*/
/*   if an EOF is detected, the error message is printed and the program*/
/*   terminates.  This is used for reads that should find something in  */
/*   the file.                                */
/*                                    */
/************************************************************************/

#define CHECK_FEOF(X, msg)                                                     \
  if (X == 0) {                                                                \
    return (DCD_BADEOF);                                                       \
  }

/*            GLOBAL VARIABLES                */

void pad(char *s, int len) {
  int curlen;
  int i;

  curlen = strlen(s);

  if (curlen > len) {
    s[len] = '\0';
    return;
  }

  for (i = curlen; i < len; i++) {
    s[i] = ' ';
  }

  s[i] = '\0';
}

/************************************************************************/
/*                                    */
/*            FUNCTION open_dcd_read                */
/*                                    */
/*   INPUTS:                                */
/*    filename - filename to read in                    */
/*                                    */
/*   OUTPUTS:                                */
/*    an open filedesriptor (integer) is returned if the call is      */
/*   successful, otherwise a negative number is returned        */
/*                                    */
/*    open_dcd_read opens a dcd file for input.  It first does a check*/
/*   to see if the file really exists.  If the file does not exist,     */
/*   a DCD_DNE is returned, if the file exists but can' be opened for   */
/*   some reason, a DCD_OPENFAILED is returned.  If the open is         */
/*   successful, the file descriptor is returned.            */
/*                                    */
/************************************************************************/

int open_dcd_read(char *filename)

{
  int dcdfd; /*  file descriptor for dcd file    */

  /*  Do a stat just to see if the file really exists    */
  if (access(filename, F_OK) != 0) {
    if (errno == ENOENT) {
      return (DCD_DNE);
    }
  }

  /*  Try and open the file                */
#ifdef WIN32
  dcdfd = _open(filename, O_RDONLY | O_BINARY | O_LARGEFILE);
#else
  dcdfd = open(filename, O_RDONLY | O_LARGEFILE);
#endif

  if (dcdfd == -1) {
    return (DCD_OPENFAILED);
  }

  return (dcdfd);
}

/****************************************************************/
/*                                */
/*            FUNCTION read_dcdheader            */
/*                                */
/*   INPUTS:                            */
/*    fd - file descriptor for the dcd file            */
/*    N - Number of atoms                    */
/*    NSET - Number of sets of coordinates            */
/*    ISTART - Starting timestep of DCD file            */
/*    NSAVC - Timesteps between DCD saves            */
/*    DELTA - length of a timestep                */
/*                                */
/*   OUTPUTS:                            */
/*    N, NSET, ISTART, NSAVC, and DELTA are returned populated*/
/*   a 0 is returned if the call is successful, or a negative   */
/*   number if errors are detected                */
/*                                */
/*    read_header reads in the header from a DCD file and     */
/*   returns the timestep size and the number of atoms in the   */
/*   system.  A 0 is returned if the header is successfully     */
/*   read, or a negative error code is returned otherwise.      */
/*                                */
/****************************************************************/

int read_dcdheader(int fd, int *N, int *NSET, int *ISTART, int *NSAVC,
                   double *DELTA, int *NAMNF, int **FREEINDEXES)

{
  int32 input_integer; /*  Integer buffer space    */
  float input_float;
  char bigbuf[256]; /*  A large string buffer    */
  int ret_val;      /*  Return value from read    */
  int i;            /*  Loop counter        */
  char HDR[5];      /*  Header = "CORD"        */
  int32 I;
  int32 NTITLE;

  /*  First thing in the file should be an 84        */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading first int from dcd file");
  CHECK_FEOF(ret_val, "reading first int from dcd file");

  if (input_integer != 84) {
    return (DCD_BADFORMAT);
  }

  /*  Read in the string "CORD"            */
  ret_val = read(fd, HDR, 4);

  CHECK_FREAD(ret_val, "reading CORD from dcd file");
  CHECK_FEOF(ret_val, "reading CORD from dcd file");

  HDR[ret_val] = '\0';

  if (strcmp(HDR, "CORD") != 0) {
    return (DCD_BADFORMAT);
  }

  /*  Read in the number of Sets of coordinates, NSET  */
  ret_val = read(fd, &input_integer, sizeof(int32));
  *NSET = input_integer;

  CHECK_FREAD(ret_val, "reading NSET");
  CHECK_FEOF(ret_val, "reading NSET");

  /*  Read in ISTART, the starting timestep         */
  ret_val = read(fd, &input_integer, sizeof(int32));
  *ISTART = input_integer;

  CHECK_FREAD(ret_val, "reading ISTART");
  CHECK_FEOF(ret_val, "reading ISTART");

  /*  Read in NSAVC, the number of timesteps between   */
  /*  dcd saves                         */
  ret_val = read(fd, &input_integer, sizeof(int32));
  *NSAVC = input_integer;

  CHECK_FREAD(ret_val, "reading NSAVC");
  CHECK_FEOF(ret_val, "reading NSAVC");

  /*  Skip blank integers                     */
  for (i = 0; i < 5; i++) {
    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading I");
    CHECK_FEOF(ret_val, "reading I");
  }

  /*  Read NAMNF, the number of free atoms         */
  ret_val = read(fd, &input_integer, sizeof(int32));
  *NAMNF = input_integer;

  CHECK_FREAD(ret_val, "reading NAMNF");
  CHECK_FEOF(ret_val, "reading NAMNF");

  /*  Read in the timestep, DELTA                */
  ret_val = read(fd, &input_float, sizeof(float));

  CHECK_FREAD(ret_val, "reading DELTA");
  CHECK_FEOF(ret_val, "reading DELTA");
  *DELTA = input_float;

  /*  Skip blank integers                    */
  for (i = 0; i < 10; i++) {
    ret_val = read(fd, &I, sizeof(int32));

    CHECK_FREAD(ret_val, "reading I");
    CHECK_FEOF(ret_val, "reading I");
  }

  /*  Get the end size of the first block            */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading second 84 from dcd file");
  CHECK_FEOF(ret_val, "reading second 84 from dcd file");

  if (input_integer != 84) {
    return (DCD_BADFORMAT);
  }

  /*  Read in the size of the next block            */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading size of title block");
  CHECK_FEOF(ret_val, "reading size of title block");

  if (((input_integer - 4) % 80) == 0) {
    /*  Read NTITLE, the number of 80 characeter    */
    /*  title strings there are            */
    ret_val = read(fd, &NTITLE, sizeof(int32));

    CHECK_FREAD(ret_val, "reading NTITLE");
    CHECK_FEOF(ret_val, "reading NTITLE");

    for (i = 0; i < NTITLE; i++) {
      ret_val = read(fd, bigbuf, 80);

      CHECK_FREAD(ret_val, "reading TITLE");
      CHECK_FEOF(ret_val, "reading TITLE");
    }

    /*  Get the ending size for this block        */
    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading size of title block");
    CHECK_FEOF(ret_val, "reading size of title block");
  } else {
    return (DCD_BADFORMAT);
  }

  /*  Read in an 4                */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading an 4");
  CHECK_FEOF(ret_val, "reading an 4");

  if (input_integer != 4) {
    return (DCD_BADFORMAT);
  }

  /*  Read in the number of atoms            */
  ret_val = read(fd, &input_integer, sizeof(int32));
  *N = input_integer;

  CHECK_FREAD(ret_val, "reading number of atoms");
  CHECK_FEOF(ret_val, "reading number of atoms");

  /*  Read in an 4                */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading an 4");
  CHECK_FEOF(ret_val, "reading an 4");

  if (input_integer != 4) {
    return (DCD_BADFORMAT);
  }

  if (*NAMNF != 0) {
    (*FREEINDEXES) = new int[(*N) - (*NAMNF)];
    int32 *freeindexes32 = new int32[(*N) - (*NAMNF)];

    if (*FREEINDEXES == NULL || freeindexes32 == NULL)
      return (DCD_BADMALLOC);

    /*  Read in an size                */
    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");

    if (input_integer != ((*N) - (*NAMNF)) * 4) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, freeindexes32, ((*N) - (*NAMNF)) * sizeof(int32));

    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");

    for (i = 0; i < ((*N) - (*NAMNF)); ++i)
      (*FREEINDEXES)[i] = freeindexes32[i];

    delete[] freeindexes32;

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading size of index array");
    CHECK_FEOF(ret_val, "reading size of index array");

    if (input_integer != ((*N) - (*NAMNF)) * 4) {
      return (DCD_BADFORMAT);
    }
  }

  return (0);
}

/************************************************************************/
/*                                    */
/*            FUNCTION read_dcdstep                */
/*                                    */
/*   INPUTS:                                */
/*    fd - file descriptor to use                    */
/*    N - Number of atoms                        */
/*    X - array of X coordinates                    */
/*    Y - array of Y coordinates                    */
/*    Z - array of Z coordinates                    */
/*    num_fixed - Number of fixed atoms                */
/*    first - flag 1->first time called                */
/*    indexes - indexes for free atoms                */
/*                                    */
/*   OUTPUTS:                                */
/*    read step populates the arrays X, Y, Z and returns a 0 if the   */
/*   read is successful.  If the read fails, a negative error code is   */
/*   returned.                                */
/*                                    */
/*    read_step reads in the coordinates from one step.  It places    */
/*   these coordinates into the arrays X, Y, and Z.            */
/*                                    */
/************************************************************************/

int read_dcdstep(int fd, int N, float *X, float *Y, float *Z, int num_fixed,
                 int first, int *indexes)

{
  int ret_val;         /*  Return value from read        */
  int32 input_integer; /*  Integer buffer space        */
  int i;               /*  Loop counter            */
  static float *tmpX;

  if (first && num_fixed) {
    tmpX = new float[N - num_fixed];

    if (tmpX == NULL) {
      return (DCD_BADMALLOC);
    }
  }

  /*  Get the first size from the file                */
  ret_val = read(fd, &input_integer, sizeof(int32));

  CHECK_FREAD(ret_val, "reading number of atoms at begining of step");

  /*  See if we've reached the end of the file            */
  if (ret_val == 0) {
    delete[] tmpX;

    return (-1);
  }

  if ((num_fixed == 0) || first) {
    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, X, 4 * N);

    CHECK_FREAD(ret_val, "reading X array");
    CHECK_FEOF(ret_val, "reading X array");

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");

    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");

    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, Y, 4 * N);

    CHECK_FREAD(ret_val, "reading Y array");
    CHECK_FEOF(ret_val, "reading Y array");

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");

    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");

    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, Z, 4 * N);

    CHECK_FREAD(ret_val, "reading Z array");
    CHECK_FEOF(ret_val, "reading Z array");

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after Z array");
    CHECK_FEOF(ret_val, "reading number of atoms after Z array");

    if (input_integer != 4 * N) {
      return (DCD_BADFORMAT);
    }
  } else {
    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, tmpX, 4 * (N - num_fixed));

    CHECK_FREAD(ret_val, "reading Xtmp array");
    CHECK_FEOF(ret_val, "reading Xtmp array");

    for (i = 0; i < N - num_fixed; i++) {
      X[indexes[i] - 1] = tmpX[i];
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after X array");
    CHECK_FEOF(ret_val, "reading number of atoms after X array");

    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, tmpX, 4 * (N - num_fixed));

    CHECK_FREAD(ret_val, "reading Xtmp array");
    CHECK_FEOF(ret_val, "reading Xtmp array");

    for (i = 0; i < N - num_fixed; i++) {
      Y[indexes[i] - 1] = tmpX[i];
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after Y array");
    CHECK_FEOF(ret_val, "reading number of atoms after Y array");

    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }

    ret_val = read(fd, tmpX, 4 * (N - num_fixed));

    CHECK_FREAD(ret_val, "reading Xtmp array");
    CHECK_FEOF(ret_val, "reading Xtmp array");

    for (i = 0; i < N - num_fixed; i++) {
      Z[indexes[i] - 1] = tmpX[i];
    }

    ret_val = read(fd, &input_integer, sizeof(int32));

    CHECK_FREAD(ret_val, "reading number of atoms after Z array");
    CHECK_FEOF(ret_val, "reading number of atoms after Z array");

    if (input_integer != 4 * (N - num_fixed)) {
      return (DCD_BADFORMAT);
    }
  }

  return (0);
}

#ifndef OUTPUT_SINGLE_FILE
#error OUTPUT_SINGLE_FILE not defined!
#endif

#if OUTPUT_SINGLE_FILE
#define NFILE_POS ((OFF_T)8)
#define NPRIV_POS ((OFF_T)12)
#define NSAVC_POS ((OFF_T)16)
#define NSTEP_POS ((OFF_T)20)
#else
// Need to consider extra fields:
// 1. magic number (int32)
// 2. file version (float)
// 3. number of files that contain portion of data (int32)
// So the total offset is 12 bytes.
#define NFILE_POS ((OFF_T)20)
#define NPRIV_POS ((OFF_T)24)
#define NSAVC_POS ((OFF_T)28)
#define NSTEP_POS ((OFF_T)32)
#endif

/*********************************************************************/
/*                                     */
/*            FUNCTION open_dcd_write                 */
/*                                     */
/*   INPUTS:                                 */
/*    dcdfile - Name of the dcd file                     */
/*                                     */
/*   OUTPUTS:                                 */
/*    returns an open file descriptor for writing             */
/*                                     */
/*    This function will open a dcd file for writing.  It takes    */
/*   the filename to open as its only argument.     It will return a    */
/*   valid file descriptor if successful or DCD_OPENFAILED if the    */
/*   open fails for some reason.  If the file specified already       */
/*   exists, it is renamed by appending .BAK to it.             */
/*                                     */
/*********************************************************************/

int open_dcd_write(const char *dcdname)

{
  int dcdfd;
  NAMD_backup_file(dcdname, ".BAK");

#ifdef WIN32
  while ((dcdfd =
              _open(dcdname, O_RDWR | O_CREAT | O_EXCL | O_BINARY | O_LARGEFILE,
                    _S_IREAD | _S_IWRITE)) < 0)
#else
#ifdef NAMD_NO_O_EXCL
  while ((dcdfd = open(dcdname, O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE,
#else
  while ((dcdfd = open(dcdname, O_RDWR | O_CREAT | O_EXCL | O_LARGEFILE,
#endif
                       S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) < 0)
#endif
  {
    if (errno != EINTR)
      return (DCD_OPENFAILED);
  }

  return (dcdfd);
}

// non master nodes open an existing file
// master will have taken care of backup and creation before the
// broadcast which triggers the slaves
int open_dcd_write_par_slave(char *dcdname)

{
#if OUTPUT_SINGLE_FILE
  // In the case of single file, the dcd output by slaves has been created
  // by the master, so the dcd file doesn't need to be created again and
  // backed up. --Chao Mei
  int dcdfd;
#ifdef WIN32
  while ((dcdfd = _open(dcdname, O_WRONLY | O_BINARY | O_LARGEFILE,
                        _S_IREAD | _S_IWRITE)) < 0)
#else
  while ((dcdfd = open(dcdname, O_WRONLY | O_LARGEFILE,
                       S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) < 0)
#endif
  {
    if (errno != EINTR)
      return (DCD_OPENFAILED);
  }

  return (dcdfd);
#else
  // In the case of multiple output files, each slave has its own output file
  // and it needs to be created. --Chao Mei
  return open_dcd_write(dcdname);
#endif
}

/************************************************************************/
/*                                    */
/*                FUNCTION write_dcdstep            */
/*                                    */
/*   INPUTS:                                */
/*    fd - file descriptor for the DCD file to write to        */
/*    N - Number of atoms                        */
/*    X - X coordinates                        */
/*    Y - Y coordinates                        */
/*    Z - Z coordinates                        */
/*  unitcell - a, b, c, alpha, beta, gamma of unit cell */
/*                                    */
/*   OUTPUTS:                                */
/*    none                                */
/*                                    */
/*    write_dcdstep writes the coordinates out for a given timestep   */
/*   to the specified DCD file.                        */
/*                                                                      */
/************************************************************************/

int write_dcdstep(int fd, int N, float *X, float *Y, float *Z, double *cell)

{
  int32 NSAVC, NSTEP, NFILE;
  int32 out_integer;

  /* Unit cell */
  if (cell) {
    out_integer = 48;
    NAMD_write(fd, (char *)&out_integer, sizeof(int32));
    NAMD_write(fd, (char *)cell, out_integer);
    NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  }

  /* Coordinates */
  // Note: the value of out_integer wraps for N >= 2^30.
  out_integer = N * 4;
  // Use a separate byte count stored without overflow.
  size_t nbytes = ((size_t)N) * 4;

  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)X, nbytes);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)Y, nbytes);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)Z, nbytes);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  /* don't update header until after write succeeds */
  OFF_T end = LSEEK(fd, 0, SEEK_CUR);
  LSEEK(fd, NSAVC_POS, SEEK_SET);
  size_t res = READ(fd, (void *)&NSAVC, sizeof(int32));
  LSEEK(fd, NSTEP_POS, SEEK_SET);
  res = READ(fd, (void *)&NSTEP, sizeof(int32));
  LSEEK(fd, NFILE_POS, SEEK_SET);
  res = READ(fd, (void *)&NFILE, sizeof(int32));
  NSTEP += NSAVC;
  NFILE += 1;
  LSEEK(fd, NSTEP_POS, SEEK_SET);
  NAMD_write(fd, (char *)&NSTEP, sizeof(int32));
  LSEEK(fd, NFILE_POS, SEEK_SET);
  NAMD_write(fd, (char *)&NFILE, sizeof(int32));
  LSEEK(fd, end, SEEK_SET);

  return (0);
}

int write_dcdstep_par_cell(int fd, double *cell) {
  if (cell) {
    int32 out_integer = 48;
    NAMD_write(fd, (char *)&out_integer, sizeof(int32));
    NAMD_write(fd, (char *)cell, out_integer);
    NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  }
  return 0;
}

/*   No useful file format description of the DCD format can be found.  */
/*   The closest approximation purporting to be a format description is */
/*   a block of fortran formated statements.  Ultra lame.               */

/*   Therefore I simply reverse engineered the sequential output.       */
/*    -EJB */

/*   Notionally, one creates the file on the master with header and cell*/
/*   Master then broadcasts a command for the slaves to write.  They    */
/*   contribute to a reduction which roots at the master.  Which then   */
/*   updates the header.                                                */
/*   Update the X/Y/Z headers for each X,Y and Z fields starting the current */
/*   position of the file                                                    */
int write_dcdstep_par_XYZUnits(int fd, int N) {
  int32 out_integer;

  // number of elements
  // Note: the value of out_integer wraps for N >= 2^30.
  out_integer = N * sizeof(float);

  // For byte seeking, use a properly sized variable.
  OFF_T nbytes = ((OFF_T)N) * sizeof(float);

  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  // seek to the end of each x y z block and write out the count
  LSEEK(fd, nbytes, SEEK_CUR);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  LSEEK(fd, nbytes, SEEK_CUR);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  LSEEK(fd, nbytes, SEEK_CUR);
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  return (0);
}

/*  A function to update the header.                                      */
/*  Once all slaves have written the master updates the header.           */
int update_dcdstep_par_header(int fd) {
  int32 NSAVC, NSTEP, NFILE;
  /* don't update header until after write succeeds */
  OFF_T end = LSEEK(fd, 0, SEEK_CUR);
  LSEEK(fd, NSAVC_POS, SEEK_SET);
  size_t res = READ(fd, (void *)&NSAVC, sizeof(int32));
  LSEEK(fd, NSTEP_POS, SEEK_SET);
  res = READ(fd, (void *)&NSTEP, sizeof(int32));
  LSEEK(fd, NFILE_POS, SEEK_SET);
  res = READ(fd, (void *)&NFILE, sizeof(int32));
  NSTEP += NSAVC;
  NFILE += 1;
  LSEEK(fd, NSTEP_POS, SEEK_SET);
  NAMD_write(fd, (char *)&NSTEP, sizeof(int32));
  LSEEK(fd, NFILE_POS, SEEK_SET);
  NAMD_write(fd, (char *)&NFILE, sizeof(int32));
  LSEEK(fd, end, SEEK_SET);

  return (0);
}

/* Each slave writer has a sequential block of atoms to output */
/* The stream file position needs to be put at the beginning of each */
/* timestep X/Y/Z output                                             */
int write_dcdstep_par_slave(int fd, int parL, int parU, int N, float *X,
                            float *Y, float *Z) {
  /* Coordinates for the N elements handled by this writer */
  int parN = parU - parL + 1;
  size_t nbytes = ((size_t)parN) * 4;

  /* x's 1st number of Xs */
  /* skip field for the bytes in X, and the first parL atoms in X*/
  OFF_T xoffset = sizeof(int) + sizeof(float) * ((OFF_T)parL);
  LSEEK(fd, xoffset, SEEK_CUR);
  NAMD_write(fd, (char *)X, nbytes);

  /* skip field for the bytes in X and Y; */
  /* skip the remaining atoms in X at number of (N-1)-(parU+1)+1
   * where N-1 is the last atom id, parU+1 is the next atom id. */
  /* skip the first parL atoms in Y; */
  OFF_T yoffset =
      2 * sizeof(int) + sizeof(float) * ((OFF_T)(N - parU + parL - 1));
  LSEEK(fd, yoffset, SEEK_CUR);
  NAMD_write(fd, (char *)Y, nbytes);

  OFF_T zoffset = yoffset;
  LSEEK(fd, zoffset, SEEK_CUR);
  NAMD_write(fd, (char *)Z, nbytes);
  return (0);
}

/*****************************************************************************/
/*                                         */
/*                FUNCTION write_dcdheader             */
/*                                         */
/*   INPUTS:                                     */
/*    fd - file descriptor for the dcd file                     */
/*    filename - filename for output                         */
/*    N - Number of atoms                             */
/*    NFILE - Number of sets of coordinates                     */
/*    NPRIV - Starting timestep of DCD file - NOT ZERO             */
/*    NSAVC - Timesteps between DCD saves                     */
/*    NSTEP - Number of timesteps                         */
/*    DELTA - length of a timestep                         */
/*                                         */
/*   OUTPUTS:                                     */
/*    none                                     */
/*                                         */
/*    This function prints the "header" information to the DCD file.  Since*/
/*   this is duplicating an unformatted binary output from FORTRAN, its ugly.*/
/*   So if you're squeamish, don't look.                         */
/* Whenever you are updating this function by removing or adding fields, please
 */
/* also update get_dcdheader_size function to reflect the changes! */
/*                                         */
/*****************************************************************************/

int write_dcdheader(int fd, const char *filename, int N, int NFILE, int NPRIV,
                    int NSAVC, int NSTEP, double DELTA, int with_unitcell) {
  int32 out_integer;
  float out_float;
  char title_string[200];
  int user_id;
#ifndef WIN32
  struct passwd *pwbuf;
#endif
  time_t cur_time;
  struct tm *tmbuf;
  char time_str[11];

  out_integer = 84;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  strcpy(title_string, "CORD");
  NAMD_write(fd, title_string, 4);
  out_integer = NFILE; /* located at fpos 8 */
  out_integer = 0;     /* ignore the lies */
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = NPRIV; /* located at fpos 12 */
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = NSAVC; /* located at fpos 16 */
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = NSTEP;         /* located at fpos 20 */
  out_integer = NPRIV - NSAVC; /* ignore the lies */
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 0;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_float = DELTA;
  NAMD_write(fd, (char *)&out_float, sizeof(float));
  out_integer = with_unitcell ? 1 : 0;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 0;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 24; // PRETEND TO BE CHARMM24 -JCP
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 84;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  out_integer = 164;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 2;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  sprintf(title_string, "REMARKS FILENAME=%s CREATED BY NAMD", filename);
  pad(title_string, 80);
  NAMD_write(fd, title_string, 80);

  char username[100];
#if defined(WIN32) || defined(NO_GETPWUID)
  sprintf(username, "Win32");
#else
  user_id = (int)getuid();
  pwbuf = getpwuid(user_id);
  if (pwbuf)
    sprintf(username, "%s", pwbuf->pw_name);
  else
    sprintf(username, "%d", user_id);
#endif
  cur_time = time(NULL);
  tmbuf = localtime(&cur_time);
  strftime(time_str, 10, "%m/%d/%y", tmbuf);

  sprintf(title_string, "REMARKS DATE: %s CREATED BY USER: %s", time_str,
          username);
  pad(title_string, 80);
  NAMD_write(fd, title_string, 80);
  out_integer = 164;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 4;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = N;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));
  out_integer = 4;
  NAMD_write(fd, (char *)&out_integer, sizeof(int32));

  return (0);
}

int get_dcdheader_size() {
  int headersize = 0;
  int totalInt32s = 27; /* 27 writes from out_integer */
  int totalChars = 164; /* 3 writes from title_string */
  int totalFloats = 1;  /* 1 write from out_float */
  headersize =
      sizeof(int32) * totalInt32s + totalChars + sizeof(float) * totalFloats;
  return headersize;
}

/****************************************************************/
/*                                */
/*            FUNCTION close_dcd_read            */
/*                                */
/*   INPUTS:                            */
/*    fd - file descriptor to close                */
/*    num_fixed - the number of fixed atoms            */
/*    indexes - Array of indexes to be deallocated        */
/*                                */
/*   OUTPUTS:                            */
/*    the file pointed to by fd is closed and the memory      */
/*   pointed to by indexes is freed if any was allocated    */
/*                                */
/*    close_dcd_read closes a dcd file that was opened for    */
/*   reading.  It also deallocated memory used for the indexes  */
/*   used for the free atom list, if there were fixed atoms.    */
/*                                */
/****************************************************************/

void close_dcd_read(int fd, int num_fixed, int *indexes)

{
#ifdef WIN32
  _close(fd);
#else
  close(fd);
#endif

  if (num_fixed) {
    delete[] indexes;
  }
}

/****************************************************************/
/*                                */
/*            FUNCTION close_dcd_write        */
/*                                */
/*   INPUTS:                            */
/*    fd - file descriptor to close                */
/*                                */
/*   OUTPUTS:                            */
/*    the file pointed to by fd                */
/*                                */
/*    close_dcd_write close a dcd file that was opened for    */
/*   writing                            */
/*                                */
/****************************************************************/

void close_dcd_write(int fd)

{
#ifdef WIN32
  if (_close(fd))
#else
  if (fsync(fd) || close(fd))
#endif
  {
    NAMD_err("Error closing DCD file");
  }
}

void read_binary_file(const char *fname, XYZ *data, int n) {
  int32 filen; //  Number of atoms read from file
  FILE *fp;    //  File descriptor
  int needToFlip = 0;

  std::cout << "Info: Reading from binary file " << fname << "\n" << std::endl;

  //  Open the file and die if the open fails
  if ((fp = Fopen(fname, "rb")) == NULL) {
    char errmsg[256];
    sprintf(errmsg, "Unable to open binary file %s", fname);
    NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  if (fread(&filen, sizeof(int32), 1, fp) != (size_t)1) {
    char errmsg[256];
    sprintf(errmsg, "Error reading binary file %s", fname);
    NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  //  check for palindromic number of atoms
  char lenbuf[4];
  memcpy(lenbuf, (const char *)&filen, 4);
  char tmpc;
  tmpc = lenbuf[0];
  lenbuf[0] = lenbuf[3];
  lenbuf[3] = tmpc;
  tmpc = lenbuf[1];
  lenbuf[1] = lenbuf[2];
  lenbuf[2] = tmpc;
  if (!memcmp((const char *)&filen, lenbuf, 4)) {
    std::cout << "Warning: Number of atoms in binary file " << fname
              << " is palindromic, assuming same endian.\n"
              << std::endl;
  }

  //  Die if this doesn't match the number in our system
  if (filen != n) {
    needToFlip = 1;
    memcpy((char *)&filen, lenbuf, 4);
  }
  if (filen != n) {
    char errmsg[256];
    sprintf(errmsg, "Incorrect atom count in binary file %s", fname);
    NAMD_die(errmsg);
  }

  if (fread(data, sizeof(XYZ), n, fp) != (size_t)n) {
    char errmsg[256];
    sprintf(errmsg, "Error reading binary file %s", fname);
    NAMD_die(errmsg);
  }

  Fclose(fp);

  if (needToFlip) {
    std::cout << "Converting binary file " << fname << "\n" << std::endl;
    int i;
    char *cdata = (char *)data;
    for (i = 0; i < 3 * n; ++i, cdata += 8) {
      char tmp0, tmp1, tmp2, tmp3;
      tmp0 = cdata[0];
      tmp1 = cdata[1];
      tmp2 = cdata[2];
      tmp3 = cdata[3];
      cdata[0] = cdata[7];
      cdata[1] = cdata[6];
      cdata[2] = cdata[5];
      cdata[3] = cdata[4];
      cdata[7] = tmp0;
      cdata[6] = tmp1;
      cdata[5] = tmp2;
      cdata[4] = tmp3;
    }
  }

  std::cout << "Info: Finished reading from binary file " << fname << "\n"
            << std::endl;
}

/***************************************************************************
 Fopen(char *Filename, char *mode):  similar to fopen(filename,mode) except
 it checks for compressed file names too.
 For example:  Fopen("config");
   This will first look for the filename "config" (and return "r" file handle
   if it is found).
   Then it will look for "config.Z" (and run "zcat config.Z", returning
   a file handle to the uncompressed data if found).
   Then it will look for "config.gz" (and run "gzip -d -c config.gz", returning
   a file handle to the uncompressed data if found).
 ***************************************************************************/
FILE *Fopen(const char *filename, const char *mode) {
  struct stat buf;
  // check if basic filename exists (and not a directory)

#if defined(NOCOMPRESSED)
  if (!stat(filename, &buf)) {
    FILE *rval;
    while (!(rval = fopen(filename, mode))) {
      if (errno != EINTR)
        break;
    }
    return (rval);
  }
#else
  if (!stat(filename, &buf)) {
    if (!S_ISDIR(buf.st_mode)) {
      FILE *rval;
      while (!(rval = fopen(filename, mode))) {
        if (errno != EINTR)
          break;
      }
      return (rval);
    }
  }
  // check for a compressed file
  /*
  char *realfilename;
  char *command;
  FILE *fout;
  command = (char *)malloc(strlen(filename)+25);
  // check for .Z (unix compress)
  sprintf(command,"zcat %s.Z",filename);
  realfilename = command+5;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.Z = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
    {
    if (!S_ISDIR(buf.st_mode))
        {
        fout = popen(command,mode);
        // on HP-UX, the first character(s) out of pipe may be
        // garbage!  (Argh!)
        int C;
        do
          {
          C = fgetc(fout);
          // iout << "C is " << C << "\n" << endi;
          if (isalnum(C) || isspace(C))
            {
            ungetc(C,fout);
            C = -1;    // outta loop
            }
          } while(C != -1);
        free(command);
        return(fout);
        }
    }
  // check for .gz (gzip)
  sprintf(command,"gzip -d -c %s.gz",filename);
  realfilename = command+11;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.gz = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
    {
    if (!S_ISDIR(buf.st_mode))
        {
        fout = popen(command,mode);
        // on HP-UX, the first character(s) out of pipe may be
        // garbage!  (Argh!)
        int C;
        do
          {
          C = fgetc(fout);
          // iout << "C is " << C << "\n" << endi;
          if (isalnum(C) || isspace(C))
            {
            ungetc(C,fout);
            C = -1;    // outta loop
            }
          } while(C != -1);
        free(command);
        return(fout);
        }
    }
  free(command);
  */
#endif /* !defined(NOCOMPRESSED) */

  return (NULL);
} /* Fopen() */

/***************************************************************************
 Fclose(FILE *fout):  similar to fclose(fout) except it first checks if the
 file handle fout is a named pipe.
 ***************************************************************************/
int Fclose(FILE *fout) {
  int rc = -1;
#if !defined(NOCOMPRESSED)
  rc = pclose(fout);
#endif
  if (rc == -1) // stream not associated with a popen()
  {
    rc = fclose(fout);
  }
  return rc;
} /* Fclose() */