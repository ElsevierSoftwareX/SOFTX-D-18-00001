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

/** @file ValidPtr.h Smart pointer class to control access to object.
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date November 2006
 *
 */
#ifndef MAT_VALIDPTR
#define MAT_VALIDPTR
namespace mat {


  /** Smart pointer class to control access to object 
   *
   * Primary use:
   * Control access to objects that may be written to file.
   *
   */
  template <typename Tobj>
    class ValidPtr {
  public:
    /** Copy ordinary pointer constructor */
    explicit ValidPtr(Tobj * p) 
      : ptr(p), inMemory(true), haveDataStructure(false){}
    ~ValidPtr() {
      delete ptr;
    }
    
      /* Pointer can not be changed only object pointed to. 
       * Therefore this is a const operation.
       * Note that if Tobj is const it can not be changed of course.
       */
      Tobj & operator*() const {
	if (!inMemory)
	  throw Failure("ValidPtr::operator*() const: "
			"Attempt to access invalid object. "
			"Object is on file.");
	if (!haveDataStructure)
	  throw Failure("ValidPtr::operator*() const: "
			"Attempt to access invalid object. "
			"Do not have data structure.");
	return *ptr;
      }

      Tobj * operator->() const {
	if (!inMemory)
	  throw Failure("ValidPtr::operator->() const: "
			"Attempt to access invalid pointer."
			"Object is on file.");
	if (!haveDataStructure)
	  throw Failure("ValidPtr::operator->() const: "
			"Attempt to access invalid pointer. "
			"Do not have data structure.");
	return ptr;
      }

      /** getConstRefForCopying() is provided to make it possible to
	  copy the object also when it is written to file. */
      const Tobj & getConstRefForCopying() const {
	return *ptr;
      }

      inline void inMemorySet(bool val) {
	inMemory = val;
      }
      inline bool inMemoryGet() const {
	return inMemory;
      }
      inline void haveDataStructureSet(bool val) {
	haveDataStructure = val;
      }
      inline bool haveDataStructureGet() const {
	return haveDataStructure;
      }

      static void swap( ValidPtr<Tobj> & ptrA, ValidPtr<Tobj> & ptrB ) {
	// For the moment, we do not allow swapping ptrs with objs
	// written to file. This could be a feature to add but would
	// require swapping filenames.
	if ( !ptrA.inMemoryGet() ||  !ptrB.inMemoryGet() ) 
	  throw "Swap called for objects not in memory";
	if ( !ptrA.haveDataStructureGet() ||  !ptrB.haveDataStructureGet() )
	  throw "Swap called for objects without data structure";
	Tobj * tmpPtr = ptrA.ptr;
	ptrA.ptr = ptrB.ptr;
	ptrB.ptr = tmpPtr;
      }
  protected:
      Tobj * ptr;
      /** Access to ptr forbidden if inMemory is false */ 
      bool inMemory; 
      /** Access to ptr forbidden if haveDataStructure is false */
      bool haveDataStructure; 
  private:
      ValidPtr<Tobj>(ValidPtr<Tobj> const &) {}
      ValidPtr<Tobj>& operator=(ValidPtr<Tobj> const &) {}
  };
  
}  /* end namespace mat */
#endif
