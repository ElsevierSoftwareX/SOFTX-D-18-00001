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

/** @file FileWritable.h Abstract class for simple writing and reading of 
 *  objects to/from file.
 *
 * @see mat::FileWritable
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2006
 *
 */
#ifndef MAT_FILEWRITABLE
#define MAT_FILEWRITABLE
#include <map>
#include <set>
namespace mat {
  /** Write and read objects to/from file. 
   *
   * This is an abstract class.
   * Classes that are derived from this class must define the 
   * following pure virtual functions to be able to instantiate objects:
   * - clear()
   * - write_to_buffer_count(int&) const
   * - write_to_buffer(void*, int const) const
   * - read_from_buffer(void*, int const)
   */  
  class FileWritable {
  public:
    /** Set the path to which the objects will be written. 
     *  This function can only be called before instantiation of objects. 
     */
    static void setPath(char const * const newPath);
   
    /** Activate the filewriting.
     *  Without calling this function no filewriting will occur. 
     *  This function can only be called before instantiation of objects. 
     * 
     */
    static void activate();
    /* FIXME: Make it possible to call activate() and deactivate() at any 
     *        time. These functions will then go through the list of objects
     *        and check the objectIsOnFile flag for each of them. Some 
     *        objects will be put on file when activate() is called and some
     *        be taken from file when deactivate() is called.
     *        A static list of objects is needed for this and for the 
     *        defragmentation function.
     */

    /** Write object to file if filewrite is active.
     *  Object is "cleared" in this call.
     */
    void writeToFile();

    /** Read object from file if filewrite is active.
     */
    void readFromFile();

    void copyToFile(const char* destFileName);

    void copyFromFile(const char* sourceFileName);

    /** Check if object is on file.
     */
    bool isOnFile() { return objectIsOnFile; }
    
    /** Return file size. Call only if obj is on file. */
    long int fileSize();

    static std::string getStatsFileSizes();
    static std::string writeAndReadAll();

    static void resetStats();
    static std::string getStatsTimeWrite();
    static std::string getStatsTimeRead();
    static std::string getStatsTimeCopyAndAssign();
    static std::string getStatsCountWrite();
    static std::string getStatsCountRead();
    static std::string getStatsCountCopyAndAssign();

    
  protected:
    /** Release memory for the information written to file.
     */
    virtual void clear() = 0;
    /** Make object invalid (false) via this function when object is 
     *  written to file and valid (true) when object is read from file.
     */
    virtual void inMemorySet(bool) = 0;

    /** Write object to file. Defined in derived class. */
    virtual void writeToFileProt(std::ofstream &) const = 0;
    /** Read object from file. Defined in derived class. */
    virtual void readFromFileProt(std::ifstream &) = 0;

    FileWritable(); /**< Gives each object a unique ID-number and filename. */
    virtual ~FileWritable(); /**< Removes file, if any. */
    
    FileWritable(FileWritable const &);     
    /* Remember to call me (operator=) explicitly in derived class! */
    FileWritable& operator=(FileWritable const &);
    
    virtual std::string obj_type_id() const = 0;
    typedef std::map<std::string, double> TypeTimeMap;
    typedef std::map<std::string, int> TypeCountMap;
    static std::string getStatsTime( TypeTimeMap & theMap );
    static std::string getStatsCount( TypeCountMap & theMap );
    struct Stats {
      // This should be a singleton
      static Stats& instance() {
	static Stats stats;
	return stats;
      }
      TypeTimeMap wallTimeWrite;
      TypeTimeMap wallTimeRead;
      TypeTimeMap wallTimeCopyAndAssign;
      TypeCountMap countWrite;
      TypeCountMap countRead;
      TypeCountMap countCopyAndAssign;
    protected:
      Stats() {}
    private:
      Stats(Stats const &);
    };

    typedef std::set<FileWritable*> ObjPtrSet;
    static std::string getStatsFileSizes( ObjPtrSet const & set );
    struct Manager {
      static Manager const & instance() {
	return instance_prot();
      }
      static void registerObj(FileWritable* objPtr);
      static void unRegisterObj(FileWritable* objPtr);
      ObjPtrSet obj_ptr_set;
    protected:
      // Only members can reach a non-const set
      static Manager& instance_prot() {
	static Manager manager;
	return manager;
      }
      Manager() {}
      Manager(Manager const &);
      //      std::map<FileWritable*, bool> obj_onFile_map;
    };
    
  private:
    static unsigned int nObjects; /**<  The number of instantiated objects. 
				   * Note that the objects may be of different 
				   * types derived from this base class. */
    static char* path;   /**< The path to which files will be written. */
    static bool active;  /**< States whether the filewriting is active. */
    unsigned int const IDNumber; /**< Each object has its unique ID-number. */
    char * fileName;     /**< Each object has its unique filename. */
    bool objectIsOnFile; /**< States whether the object is on file or not. */
    
  };
  
} /* end namespace mat */

#endif
