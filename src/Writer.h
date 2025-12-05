/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef WRITER_H
#define WRITER_H

#include <stdlib.h> //for exit

#include <fstream>
#include <iostream>
#include <ostream>
#include <string> //for file names etc.

class Writer {
public:
  // Init later...
  Writer() { isOpen = false; }
  // Init now.
  Writer(std::string const &name, std::string const &alias,
         const bool crit = true, const bool note = true) {
    Init(name, alias, crit, note);
  }

  ~Writer(void) {
    if (isOpen)
      close();
  }

  // Set main class vars.
  void Init(std::string const &name, std::string const &alias, const bool crit,
            const bool note) {
    fileName = name;
    fileAlias = alias;
    critical = crit;
    notify = note;
    isOpen = false;
    firstWrite = true;
    nameWAlias = fileAlias + " file: ./" + fileName;
  }

  // Open or close a file, with basic protections
  void open(void) {
    if (isOpen)
      return;
    file.open(fileName.c_str(), (firstWrite ? std::ios::out : std::ios::app));
    CheckFileState(true, "...could not be opened.", "Writing to ");
  }

  // Open or close a file, with basic protections
  void openOverwrite(void) {
    if (isOpen)
      return;
    file.open(fileName.c_str(), std::ios::out);
    CheckFileState(true, "...could not be opened.", "Writing to ");
  }

  void close(void) {
    if (!isOpen)
      return;
    file.close();
    CheckFileState(false, "...could not be closed.", "Finished writing ");
    // If first write is complete, unset firstWrite flag.
    if (firstWrite)
      firstWrite = false;
  }

  // Make public to expose object.
  std::ofstream file;

protected:
  bool GoodFileWData(void) { return file.is_open() && file.good(); }

  void HandleError(std::string const &msg) {
    std::cerr << ((critical) ? "Error " : "Warning ") << nameWAlias << std::endl
              << msg << std::endl;
    if (critical)
      exit(1);
  }

  void HandleNote(std::string const &msg) {
    std::cout << msg << nameWAlias << std::endl;
  }

  void CheckFileState(const bool expected, std::string const &errMessage,
                      std::string const &note) {
    isOpen = GoodFileWData();
    if (isOpen == expected && notify)
      HandleNote(note);
    else if (isOpen != expected)
      HandleError(errMessage);
  }

  std::string fileName, fileAlias, nameWAlias;
  bool critical, notify, isOpen, firstWrite;
};

#endif /*WRITER_H*/
