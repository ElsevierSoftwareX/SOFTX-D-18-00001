/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/** @file FileWritable.cc Implementation of the abstract class FileWritable 
 *  for simple writing and reading of objects to/from file.
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2006
 *
 */
#include <stdlib.h> /* system */
#include <unistd.h>
#include <stdio.h>  /* For FILE sprintf*/
#include <string.h> /* For strcpy */
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ios>
#include <stdexcept>
#include "FileWritable.h"
#include "Failure.h"
#include "matInclude.h"
#include <errno.h>

namespace mat {
  unsigned int FileWritable::nObjects = 0;
  char* FileWritable::path = NULL;
  bool FileWritable::active = false;

  void FileWritable::setPath(char const * const newPath) {
    if (nObjects != 0)
      throw Failure("FileWritable::set_path(char*) : It is not allowed "
		    "to set the path after instantiation of objects.");
    if (newPath == NULL)
      throw Failure("FileWritable::set_path(char*) : newPath == NULL.");
    if (path != NULL)
      delete[] path;
    path = new char[strlen(newPath)+2];
    strcpy(path, newPath);
    strcat(path,"/");
  }
  void FileWritable::activate() {
    if (nObjects != 0)
      throw Failure("FileWritable::activate() : It is not allowed to "
		    "activate filewritable after instantiation of objects.");
    active = true;
  }

  void FileWritable::writeToFile() {
    if (objectIsOnFile)
      throw Failure("FileWritable::writeToFile(): "
		    "Object is already on file.");
    Time t; t.tic();
    if (FileWritable::active) {
      std::ofstream file(fileName, std::ios::out|std::ios::binary);
      file.exceptions( std::ofstream::eofbit  | 
		       std::ofstream::failbit | 
		       std::ofstream::badbit );
      if (!file.is_open())
	throw "FileWritable::writeToFile(): unable to open file. ";
      try {
	this->writeToFileProt(file);
	file.close();
      }
      catch (std::ofstream::failure& e) {
	/* FIXME: We should not be doing memory allocation in error handling. */
	std::string errstr("FileWritable::writeToFile(): write to file failed : " );
	errstr += e.what();
	errstr += "\n";
	errstr += "Error (strerror(errno)): ";
	errstr += strerror(errno);
	throw strdup(errstr.c_str());
      }
      // Free memory used by object.
      this->clear();
    }
    this->inMemorySet(false);
    objectIsOnFile = true;
    Stats::instance().wallTimeWrite[obj_type_id()] += t.toc();
    ++Stats::instance().countWrite[obj_type_id()];
  }

  void FileWritable::readFromFile() {
    if (!objectIsOnFile)
      throw Failure("FileWritable::readFromFile(): Object is not on file.");
    Time t; t.tic();
    this->inMemorySet(true);
    if (FileWritable::active) {
      std::ifstream file;
      // ELIAS NOTE 2012-02-24: There was a strange error where 
      // "std::ios_base::failure was thrown with what(): 'basic_ios::clear'".
      // Another try/catch was added here to try to understand what was happening.
      try {
	// open file
	file.exceptions( std::ifstream::eofbit  | 
			 std::ifstream::failbit | 
			 std::ifstream::badbit );
	file.open(fileName, std::ios::in|std::ios::binary);
      }
      catch (std::ifstream::failure& e) {
	std::stringstream ss;
	ss << "Exception std::ifstream::failure caught in FileWritable::readFromFile() when trying to open file. fileName = '" << fileName << "'. what(): '" << e.what() << "'.";
	std::string errstr = ss.str();
	errstr += "\n";
        errstr += "Error (strerror(errno)): ";
        errstr += strerror(errno);
	std::cerr << errstr << std::endl;
	throw std::runtime_error(errstr);
      }
      if (!file.is_open())
	throw "FileWritable::readFromFile(): unable to open file. ";
      try {
	this->readFromFileProt(file);
	file.close();
      }
      catch (std::ifstream::failure& e) {
	/* FIXME: We should not be doing memory allocation in error handling. */
	std::string errstr("FileWritable::readFromFile(): read from file failed : " );
	errstr += e.what();
	errstr += "\n";
        errstr += "Error (strerror(errno)): ";
        errstr += strerror(errno);
	throw strdup(errstr.c_str());
      }
      // delete file
      unlink(fileName);  /* <unistd.h> */
    }
    objectIsOnFile = false;
    Stats::instance().wallTimeRead[obj_type_id()] += t.toc();
    ++Stats::instance().countRead[obj_type_id()];
  }


  static long int get_file_size(const char* fileName);

 void FileWritable::copyFromFile(const char* sourceFileName)
 {
   clear();

   if (FileWritable::active) 
     {
       // Open source file for reading.
       std::ifstream inFile;
       inFile.exceptions( std::ifstream::eofbit  | 
			  std::ifstream::failbit | 
			  std::ifstream::badbit );
       inFile.open(sourceFileName, std::ios::in|std::ios::binary);
       if (!inFile.is_open())
	 throw std::runtime_error("FileWritable copyFromFile error: unable to open file for reading.");
       // Open destination file for writing.
       std::ofstream outFile;
         outFile.exceptions( std::ifstream::eofbit  | 
  			   std::ifstream::failbit | 
  			   std::ifstream::badbit );
          try {
            outFile.open(fileName, std::ios::out|std::ios::binary);
          }
          catch(std::ofstream::failure& e){
              throw std::runtime_error("FileWritable copyFromFile error: an error "
                                 "has occurred on the associated stream. "
                                 "Did you set a wrong path in setPath()?");}
         if (!outFile.is_open())
  	         throw std::runtime_error("FileWritable copyFromFile error: unable to open file for writing.");
         // Copy data.
         try {
	         outFile << inFile.rdbuf();
        }
       catch (std::ofstream::failure& e) {
	         throw std::runtime_error("FileWritable copyFromFile error: failed to copy data. Out of disk space?");
       }
       inFile.close();
       outFile.close();
       /* Verify that file was really copied. */
       long int sz1 = get_file_size(fileName);
       long int sz2 = get_file_size(sourceFileName);
       if(sz1 < 0 || sz2 < 0 || sz1 != sz2)
	 throw std::runtime_error("FileWritable copyFromFile error: file sizes do not match after copying.");

     }
   else
     {
       std::ifstream file;
       // open file
       file.exceptions( std::ifstream::eofbit  | 
			std::ifstream::failbit | 
			std::ifstream::badbit );
       file.open(sourceFileName, std::ios::in|std::ios::binary);

       if (!file.is_open())
	 throw "FileWritable copyFromFile error: unable to open file.";

       this->readFromFileProt(file);
       file.close();
     }

     this->inMemorySet(false);
     objectIsOnFile = true;

 }

  void FileWritable::copyToFile(const char* destFileName)
  {
    if (!objectIsOnFile)
      throw Failure("FileWritable::copyToFile(): Object is not on file.");

    if (FileWritable::active) // if ergo flag is set, it means matrix is in the temporary file, 
                              // thus just copy a file
      {
	// Open source file for reading.
	std::ifstream inFile;
	inFile.exceptions( std::ifstream::eofbit  | 
			   std::ifstream::failbit | 
			   std::ifstream::badbit );
	inFile.open(fileName, std::ios::in|std::ios::binary);
	if (!inFile.is_open())
	  throw std::runtime_error("FileWritable copyToFile: unable to open file for reading.");
	// Open destination file for writing.
	std::ofstream outFile;
	outFile.exceptions( std::ifstream::eofbit  | 
			    std::ifstream::failbit | 
			    std::ifstream::badbit );
	outFile.open(destFileName, std::ios::out|std::ios::binary);
	if (!outFile.is_open())
	  throw std::runtime_error("FileWritable copyToFile error: unable to open file for writing.");
	// Copy data.
	try {
	  outFile << inFile.rdbuf();
	}
	catch (std::ofstream::failure& e) {
	  throw std::runtime_error("FileWritable copyToFile error: failed to copy data. Out of disk space?");
	}
	inFile.close();
	outFile.close();
	/* Verify that file was really copied. */
	long int sz1 = get_file_size(fileName);
	long int sz2 = get_file_size(destFileName);
	if(sz1 < 0 || sz2 < 0 || sz1 != sz2)
	  throw std::runtime_error("FileWritable copyToFile error: file sizes do not match after copying.");
      }
    else // flag is not set, create a new file
      {
	readFromFile(); // since active is false, it will just change a flag objectIsOnFile to false
	std::ofstream file(destFileName, std::ios::out|std::ios::binary);
	file.exceptions( std::ofstream::eofbit  | 
			 std::ofstream::failbit | 
			 std::ofstream::badbit );
	if (!file.is_open())
	  throw "FileWritable copyToFile: unable to open file. ";
	try {
	  this->writeToFileProt(file);
	  file.close();
	}
	catch (std::ofstream::failure& e) {
	  throw std::runtime_error("FileWritable copyToFile: failed to copy data.");
	}
	writeToFile(); // since active is false, it will just change a flag objectIsOnFile to true
      }
  }




  long int FileWritable::fileSize() {
    if ( !isOnFile() )
      throw std::runtime_error("Attempt to get file size for object "
			       "not on file");
    if ( !FileWritable::active ) 
      return 0;
    FILE * fptr;
    fptr = fopen ( fileName , "rb" );
    if ( fptr == NULL )
      throw std::runtime_error("FileWritable::fileSize(): fptr == NULL");
    fseek( fptr, (long int)0, SEEK_END );
    long int sz = ftell( fptr );
    fclose ( fptr );
    return sz;
  }
  
  void FileWritable::resetStats() {
    Stats::instance().wallTimeWrite.clear();
    Stats::instance().wallTimeRead.clear();
    Stats::instance().wallTimeCopyAndAssign.clear();
    Stats::instance().countWrite.clear();
    Stats::instance().countRead.clear();
    Stats::instance().countCopyAndAssign.clear();
  }
  std::string FileWritable::getStatsTime(TypeTimeMap & theMap) {
    double totalTime = 0;
    std::stringstream ss;
    TypeTimeMap::iterator it;
    bool firstIter = true;
    for ( it=theMap.begin() ; it != theMap.end(); it++ ) {
      if (firstIter)
	firstIter = false;
      else
	ss << " + ";
      ss << std::setprecision(2) << (*it).second << " s (" << (*it).first << ")";
      totalTime += (*it).second;
    }    
    if (!firstIter)
      ss << " = ";
    ss << std::setprecision(2) << totalTime << " s";
    return ss.str();
  } 
  std::string FileWritable::getStatsCount(TypeCountMap & theMap) {
    int totalCount = 0;
    std::stringstream ss;
    TypeCountMap::iterator it;
    bool firstIter = true;
    for ( it=theMap.begin() ; it != theMap.end(); it++ ) {
      if (firstIter)
	firstIter = false;
      else
	ss << " + ";
      ss << (*it).second << " (" << (*it).first << ")";
      totalCount += (*it).second;
    }    
    if (!firstIter)
      ss << " = ";
    ss << totalCount << " ";
    return ss.str();
  } 
  std::string FileWritable::getStatsTimeWrite() {
    return getStatsTime( Stats::instance().wallTimeWrite );
  }
  std::string FileWritable::getStatsTimeRead() {
    return getStatsTime( Stats::instance().wallTimeRead );
  }
  std::string FileWritable::getStatsTimeCopyAndAssign() {
    return getStatsTime( Stats::instance().wallTimeCopyAndAssign );
  }
  std::string FileWritable::getStatsCountWrite() {
    return getStatsCount( Stats::instance().countWrite );
  }
  std::string FileWritable::getStatsCountRead() {
    return getStatsCount( Stats::instance().countRead );
  }
  std::string FileWritable::getStatsCountCopyAndAssign() {
    return getStatsCount( Stats::instance().countCopyAndAssign );
  }

  FileWritable::FileWritable()
    : IDNumber(nObjects), objectIsOnFile(false) {
    nObjects++;
    if (nObjects >= 2147483647) /* 2^31 - 1 */
      throw Failure("FileWritable::FileWritable(): To many (2^31 - 1) "
		    "allocated objects");
    /* Assign filename: */
    char buffer [30]; 
    sprintf(buffer, "tmp%010i.obj", IDNumber);
    //    itoa(IDNumber,buffer,10); /* stdlib.h */
    if (path != NULL) {
      fileName = new char[strlen(path) + strlen(buffer) + 10];
      strcpy(fileName, path);
      strcat(fileName, buffer);
    }
    else {
      fileName = new char[strlen(buffer) + 10];
      strcpy(fileName, buffer);
    }
    Manager::registerObj(this);
  }
  
  FileWritable::~FileWritable() {
    Manager::unRegisterObj(this);
    // Remove file if it exists
    if (FileWritable::active && objectIsOnFile)
      unlink(fileName);
    // free memory used for file name string.
    delete []fileName;
  }

  static long int get_file_size(const char* fileName) {
    FILE * fptr = fopen ( fileName , "rb" );
    if ( fptr == NULL )
      return -1;
    fseek( fptr, (long int)0, SEEK_END );
    long int sz = ftell( fptr );
    fclose ( fptr );
    return sz;
  }

  /* ELIAS NOTE 2012-02-29: Earlier, the file copying for FileWritable
     was done using a "cp" call to the system() function, but that
     turned out to be unreliable. Copying that way failed for some
     large cases like the Umeda molecule. Doing it using "outFile <<
     inFile.rdbuf()" seems to work better. */
  static void copy_file(const char* sourceFileName, const char* destFileName) {
    // Open source file for reading.
    std::ifstream inFile;
    inFile.exceptions( std::ifstream::eofbit  | 
		       std::ifstream::failbit | 
		       std::ifstream::badbit );
    inFile.open(sourceFileName, std::ios::in|std::ios::binary);
    if (!inFile.is_open())
      throw std::runtime_error("FileWritable copy_file error: unable to open file for reading.");
    // Open destination file for writing.
    std::ofstream outFile;
    outFile.exceptions( std::ifstream::eofbit  | 
			std::ifstream::failbit | 
			std::ifstream::badbit );
    outFile.open(destFileName, std::ios::out|std::ios::binary);
    if (!outFile.is_open())
      throw std::runtime_error("FileWritable copy_file error: unable to open file for writing.");
    // Copy data.
    try {
      outFile << inFile.rdbuf();
    }
    catch (std::ofstream::failure& e) {
      throw std::runtime_error("FileWritable copy_file error: failed to copy data. Out of disk space?");
    }
    inFile.close();
    outFile.close();
    /* Verify that file was really copied. */
    long int sz1 = get_file_size(sourceFileName);
    long int sz2 = get_file_size(destFileName);
    if(sz1 < 0 || sz2 < 0 || sz1 != sz2)
      throw std::runtime_error("FileWritable copy_file error: file sizes do not match after copying.");
  }

  FileWritable::FileWritable(FileWritable const & other) 
    : IDNumber(nObjects), objectIsOnFile(other.objectIsOnFile) {
    nObjects++;
    if (nObjects >= 2147483647)
      throw Failure("FileWritable::FileWritable(): To many (2^31 - 1) "
		    "allocated objects");
    Time t; t.tic();
    /* Assign filename: */
    char buffer [30]; 
    sprintf(buffer, "tmp%010i.obj", IDNumber);
    //    itoa(IDNumber,buffer,10); /* stdlib.h */
    if (path != NULL) {
      fileName = new char[strlen(path) + strlen(buffer) + 10];
      strcpy(fileName, path);
      strcat(fileName, buffer);
    }
    else {
      fileName = new char[strlen(buffer) + 10];
      strcpy(fileName, buffer);
    }
    if (FileWritable::active && objectIsOnFile) {
      copy_file(other.fileName, fileName);
    }
    Stats::instance().wallTimeCopyAndAssign[other.obj_type_id()] += t.toc();
    ++Stats::instance().countCopyAndAssign[other.obj_type_id()];
    Manager::registerObj(this);
  }
  
  FileWritable& FileWritable::operator=(FileWritable const & other) {
    Time t; t.tic();
    /* validSet(bool valid) should be called */
    if (FileWritable::active && objectIsOnFile) {
      unlink(fileName);
    }
    objectIsOnFile = other.objectIsOnFile;
    if (FileWritable::active && objectIsOnFile) {
      copy_file(other.fileName, fileName);
    }
    Stats::instance().wallTimeCopyAndAssign[other.obj_type_id()] += t.toc();
    ++Stats::instance().countCopyAndAssign[other.obj_type_id()];
    return *this;
  }


  void FileWritable::Manager::registerObj(FileWritable* objPtr) {
    ObjPtrSet & obj_ptr_set = instance_prot().obj_ptr_set;
    std::pair<ObjPtrSet::iterator,bool> ret = obj_ptr_set.insert( objPtr );
    if ( ret.second == false )
      throw std::runtime_error("Attempt to register object already in "
			       "object set.");
  }
  
  void FileWritable::Manager::unRegisterObj(FileWritable* objPtr) {
    ObjPtrSet & obj_ptr_set = instance_prot().obj_ptr_set;
    if ( obj_ptr_set.erase( objPtr ) != 1 )
      throw std::runtime_error("Attempt to unregister object not "
			       "in object set.");
  }
  
  std::string FileWritable::writeAndReadAll() {
    ObjPtrSet const & obj_ptr_set = Manager::instance().obj_ptr_set;
    ObjPtrSet objs_to_read;
    ObjPtrSet::const_iterator it;
    for ( it=obj_ptr_set.begin() ; it != obj_ptr_set.end() ; it++ ) {
      if ( !(*it)->isOnFile() ) {
	objs_to_read.insert( *it );
	(*it)->writeToFile();
      }	    
    } // end for
    // Now all objects should be on file.
    std::string str = getStatsFileSizes ( objs_to_read );
    // Now we read in the ones that were previously in memory.
    for ( it=objs_to_read.begin() ; it != objs_to_read.end() ; it++ ) 
      (*it)->readFromFile();	
    return str;
  } // end writeAndReadAll()
  

  std::string FileWritable::getStatsFileSizes() {
    return getStatsFileSizes( Manager::instance().obj_ptr_set );
  }
  
  std::string FileWritable::getStatsFileSizes( ObjPtrSet const & oset ) {
    typedef std::map<std::string, long int> TypeSizeMap;
    TypeSizeMap sizes;
    typedef std::map<std::string, int> TypeCountMap;
    TypeCountMap counts;
    TypeCountMap counts_large;
    // Collect information
    {
      ObjPtrSet::const_iterator it;
      for ( it=oset.begin() ; it != oset.end() ; it++ ) {
	if ( (*it)->isOnFile() ) {
	  std::string obj_type = (*it)->obj_type_id();
	  long int fs = (*it)->fileSize();
	  sizes[obj_type] += fs;
	  ++counts[obj_type];
	  if (fs > 1000000)
	    ++counts_large[obj_type];
	}	    
      } // end for
    }
    std::stringstream ss;
    // Print info to string
    {
      long int totalSize = 0;
      int totalCount = 0;
      int totalLargeCount = 0;
      TypeSizeMap::iterator it;
      bool firstIter = true;
      for ( it=sizes.begin() ; it != sizes.end(); it++ ) {
	if (firstIter)
	  firstIter = false;
	else
	  ss << " + ";
	ss << std::setprecision(2) << (*it).second / (double)1000000 << " MB (" << counts[(*it).first] << " "
	   << (*it).first << ", " << counts_large[(*it).first] << " > 1 MB)";
	totalSize += (*it).second;
	totalCount += counts[(*it).first];
	totalLargeCount += counts_large[(*it).first];
      }    
      if (!firstIter)
	ss << " = ";
      ss << std::setprecision(2) << totalSize / (double)1000000 << " MB (" << totalCount << " total, " << totalLargeCount << " > 1 MB)";
    }
    return ss.str();
  }
  

  

} /* end namespace mat */
