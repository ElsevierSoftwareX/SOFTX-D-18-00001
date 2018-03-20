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

/** @file matrix_proxy.h Proxy structs used by the matrix API
 *
 * This file contains proxy structs that are used by the matrix
 * API classes to enable operator syntax when using the API.
 *
 * Copyright(c) Emanuel Rubensson 2005
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date February 2005
 *
 */
#ifndef MAT_MATRIX_PROXY
#define MAT_MATRIX_PROXY

namespace mat {
  /*********** New code */
  /** This proxy expresses the result of multiplication of two objects, 
   *  of possibly different types, TX and TY.
   *  Primary application is for scalars, matrices, and transposed matrices.
   */
  template<typename TX, typename TY>
    struct XY {
      TX const & A;
      TY const & B;
      bool const tA;
      bool const tB;
      XY(TX const & AA, TY const & BB, 
	 bool const tAA = false, bool const tBB = false)
	:A(AA), B(BB), tA(tAA), tB(tBB)
      {}
    };

  /** This proxy expresses the result of multiplication of three objects, 
   *  of possibly different types, TX, TY, and TZ.
   *  Primary application is for scalars, matrices, and transposed matrices.
   */
  template<typename TX, typename TY, typename TZ>
    struct XYZ {
      TX const & A;
      TY const & B;
      TZ const & C;
      bool const tA;
      bool const tB;
      bool const tC;
      XYZ(TX const & AA, TY const & BB, TZ const & CC, 
	  bool const tAA = false, 
	  bool const tBB = false, 
	  bool const tCC = false)
	:A(AA), B(BB), C(CC), tA(tAA), tB(tBB), tC(tCC)
      {}
    };

  /** This proxy expresses the result of multiplication of three objects
   *  added to two other multiplied objects. 
   *  All objects may have different types, TX, TY, TZ, TU, and TV.
   *  Primary application is for scalars, matrices, and transposed matrices.
   */
 template<typename TX, typename TY, typename TZ, typename TU, typename TV>
   struct XYZpUV {
     TX const & A;
     TY const & B;
     TZ const & C;
     TU const & D;
     TV const & E;
     bool const tA;
     bool const tB;
     bool const tC;
     bool const tD;
     bool const tE;
     XYZpUV(TX const & AA, TY const & BB, TZ const & CC, 
	    TU const & DD, TV const & EE, 
	    bool const tAA = false, 
	    bool const tBB = false, 
	    bool const tCC = false,
	    bool const tDD = false, 
	    bool const tEE = false)
       :A(AA), B(BB), C(CC), D(DD), E(EE), 
	tA(tAA), tB(tBB), tC(tCC), tD(tDD), tE(tEE)
     {}
   };


  /** This proxy expresses the result of transposition of an object
   *  of type TX.
   *  Primary application is for matrices and transposed matrices.
   *  @see transpose(TX const &)
   */
  template<typename TX>
    struct Xtrans {
      TX const & A;
      bool const tA;
      explicit Xtrans(TX const & AA, bool const tAA = false)
	:A(AA), tA(tAA)
      {}
    };
  
  /** Transposition.
   *  Returns a transposition proxy of an object of type TX.
   *  @see Xtrans 
   */
  template<typename TX>
    inline Xtrans<TX> transpose(TX const & A) {
    return Xtrans<TX>(A,true);
  }
  /** Transposition.
   *  Returns a transposition proxy of an object of type Xtrans<TX>.
   *  Only for correct treatment of repeated transposition, 
   *  e.g. transpose(transpose(A))
   *  @see Xtrans 
   *  @see transpose(TX const &)
   */
  template<typename TX>
    inline Xtrans<TX> transpose(const Xtrans<TX>& xtrans) {
    return Xtrans<TX>(xtrans.A, !xtrans.tA);
  }

  /* Some operators   */
  /** Multiplication of two transposition proxys holding objects of
   *  type TX and TY respectively.
   *  Returns multiplication proxy XY.
   *  @see XY
   *  @see Xtrans 
   *  @see operator*(TX const &, Xtrans<TY> const &)
   *  @see operator*(Xtrans<TX> const &, TY const &)
   *  @see operator*(TX const &, TY const &)
   */
  template<typename TX, typename TY>
    inline XY<TX, TY> operator*(Xtrans<TX> const & trAA, 
				Xtrans<TY> const & trBB) {
    return XY<TX, TY>(trAA.A, trBB.A, trAA.tA, trBB.tA);
  }

  /** Multiplication of an object of type TX with a tranposition proxy  
   *  holding an object of type TY.
   *  Returns multiplication proxy XY.
   *  @see XY
   *  @see Xtrans 
   *  @see operator*(Xtrans<TX> const &, Xtrans<TY> const &)
   *  @see operator*(Xtrans<TX> const &, TY const &)
   *  @see operator*(TX const &, TY const &)
   */
  template<typename TX, typename TY>
    inline XY<TX, TY> operator*(TX const & AA, 
				Xtrans<TY> const & trBB) {
    return XY<TX, TY>(AA, trBB.A, false, trBB.tA);
  }

  /** Multiplication of a tranposition proxy holding an object of type TX 
   *  with an object of type TY.
   *  Returns multiplication proxy XY.
   *  @see XY
   *  @see Xtrans 
   *  @see operator*(Xtrans<TX> const &, Xtrans<TY> const &)
   *  @see operator*(TX const &, Xtrans<TY> const &)
   *  @see operator*(TX const &, TY const &)
   */
  template<typename TX, typename TY>
    inline XY<TX, TY> operator*(Xtrans<TX> const & trAA, 
				TY const & BB) {
    return XY<TX, TY>(trAA.A, BB, trAA.tA, false);
  }
    
  /** Multiplication of an object of type TX with an object of type TY.
   *  Returns multiplication proxy XY.
   *  @see XY
   *  @see operator*(Xtrans<TX> const &, Xtrans<TY> const &)
   *  @see operator*(TX const &, Xtrans<TY> const &)
   *  @see operator*(Xtrans<TX> const &, TY const &)
   */
  template<typename TX, typename TY>
    inline XY<TX, TY> operator*(TX const & AA, 
				TY const & BB) {
    return XY<TX, TY>(AA, BB, false, false);
  }
  
  /** Multiplication of a multiplication proxy XY with a transposition 
   *  proxy Xtrans. 
   *  Returns multiplication proxy XYZ.
   *  @see XY
   *  @see XYZ
   *  @see Xtrans
   */
  template<typename TX, typename TY, typename TZ>
    inline XYZ<TX, TY, TZ> 
    operator*(XY<TX, TY> const & AB, Xtrans<TZ> const & trCC) {
    return XYZ<TX, TY, TZ>(AB.A, AB.B, trCC.A, AB.tA, AB.tB, trCC.tA); 
  }

  /** Multiplication of a multiplication proxy XY with an object of type TZ. 
   *  Returns multiplication proxy XYZ.
   *  @see XY
   *  @see XYZ
   */
  template<typename TX, typename TY, typename TZ>
    inline XYZ<TX, TY, TZ> 
    operator*(XY<TX, TY> const & AB, TZ const & CC) {
    return XYZ<TX, TY, TZ>(AB.A, AB.B, CC, AB.tA, AB.tB, false); 
  }
  
  /** Addition of two multiplication proxys XYZ and XY.
   * Returns multiplication and addition proxy XYZpUV.
   *  @see XY
   *  @see XYZ
   *  @see XYZpUV
   */
  template<typename TX, typename TY, typename TZ, typename TU, typename TV>
    inline XYZpUV<TX, TY, TZ, TU, TV> 
    operator+(XYZ<TX, TY, TZ> const & ABC, XY<TU, TV> const & DE) {
    return XYZpUV<TX, TY, TZ, TU, TV>(ABC.A, ABC.B, ABC.C, DE.A, DE.B, ABC.tA, ABC.tB, ABC.tC, DE.tA, DE.tB);
  }

  /** This proxy expresses the result of addition of two objects, 
   *  of possibly different types, TX and TY.
   *  Primary application is for scalars, matrices, and transposed matrices.
   */
  template<typename TX, typename TY>
    struct XpY {
      const TX& A;
      const TY& B;
      XpY(const TX& AA,const TY& BB)
	:A(AA),B(BB)
      {}
    };
  /** Addition of two objects of type TX and TY.
   *  @see XpY
   */
  template<typename TX, typename TY>
    inline XpY<TX, TY> operator+(TX const & AA, TY const & BB) {
    return XpY<TX, TY>(AA, BB);
  }

  /** This proxy expresses the result of substraction of two objects, 
   *  of possibly different types, TX and TY.
   *  Primary application is for scalars, matrices, and transposed matrices.
   */
  template<typename TX, typename TY>
    struct XmY {
      const TX& A;
      const TY& B;
      XmY(const TX& AA,const TY& BB)
	:A(AA),B(BB)
      {}
    };
  /** Substraction of two objects of type TX and TY.
   *  @see XmY
   */
  template<typename TX, typename TY>
    inline XmY<TX, TY> operator-(TX const & AA, TY const & BB) {
    return XmY<TX, TY>(AA, BB);
  }


  /************* New code ends */


#if 0
  template<class MAT>
    struct M2 {
      const MAT& A;
      M2(const MAT& AA)
	:A(AA)
      {}
    };
  
  template<class MAT>
    inline M2<MAT> square(const MAT& A) {
    return M2<MAT>(A);
  }

  template<class SCAL, class MAT>
    struct SM2 {
      const SCAL alpha;
      const MAT& A;
      SM2(const MAT& AA, const SCAL a = 1)
	: A(AA), alpha(a)
      {}
      SM2(const M2<MAT>& m2)
	:A(m2.A), alpha(1)
      {}
    };

  template<class SCAL, class MAT>
    inline SM2<SCAL, MAT> 
    operator*(const SCAL s, const M2<MAT>& m2) {
    return SM2<SCAL, MAT>(m2.A, s); 
  }
  
  

  
  template<class MAT>
    struct MT {
      const MAT& A;
      const bool tA;
      MT(const MAT& AA, const bool tAA = false)
	:A(AA), tA(tAA)
      {}
    };

  template<class MAT>
    inline MT<MAT> transpose(const MAT& A) {
    return MT<MAT>(A,true);
  }
  template<class MAT>
    inline MT<MAT> transpose(const MT<MAT>& mt) {
    return MT<MAT>(mt.A, !mt.tA);
  }

  
  template<class SCAL, class MAT>
    struct SM {
      const SCAL alpha;
      const MAT& A;  
      const bool tA;
    SM(const MAT& AA, const SCAL scalar = 1, const bool tAA = false)
      :A(AA),alpha(scalar), tA(tAA)
      {}
    }; 

  template<class SCAL, class MAT>
    inline SM<SCAL, MAT> operator*(const SCAL scalar, const MT<MAT>& mta) {
    return SM<SCAL, MAT>(mta.A,scalar, mta.tA);
  }
  
  template<class SCAL, class MAT>
    inline SM<SCAL, MAT> operator*(const SCAL scalar, const MAT& AA) {
    return SM<SCAL, MAT>(AA, scalar, false);
  }

  

  template<class MAT, class MATB = MAT>
    struct MM {
      const MAT& A;
      const MATB& B;
      const bool tA;
      const bool tB;
      MM(const MAT& AA,const MATB& BB, const bool tAA, const bool tBB)
	:A(AA),B(BB), tA(tAA), tB(tBB)
      {}
    };

  template<class MAT, class MATB = MAT>
    struct MpM {
      const MAT& A;
      const MATB& B;
      MpM(const MAT& AA,const MATB& BB)
	:A(AA),B(BB)
      {}
    };

  template<class MAT, class MATB>
    inline MpM<MAT, MATB> operator+(const MAT& AA, const MATB& BB) {
    return MpM<MAT, MATB>(AA, BB);
  }
  
  /*
  template<class MAT, class MATB>
    inline MM<MAT, MATB> operator*(const MT<MAT>& mta, const MT<MATB>& mtb) { 
    return MM<MAT, MATB>(mta.A, mtb.A, mta.tA, mtb.tA);
  }
  */
  /*
  template<class MAT, class MATB>
    inline MM<MAT, MATB> operator*(const MAT& AA, const MT<MATB>& mtb) { 
    return MM<MAT, MATB>(AA, mtb.A, false, mtb.tA);
  }
  template<class MAT, class MATB>
    inline MM<MAT, MATB> operator*(const MT<MAT>& mta, const MATB& BB) { 
    return MM<MAT, MATB>(mta.A, BB, mta.tA, false);
  }
  template<class MAT, class MATB>
    inline MM<MAT, MATB> operator*(const MAT& AA, const MATB& BB) { 
    return MM<MAT, MATB>(AA, BB, false, false);
  }
  */

  template<class SCAL, class MAT, class MATB = MAT>
    struct SMM {
      const SCAL alpha;
      const MAT& A;
      const MATB& B;
      const bool tA;
      const bool tB;
      SMM(const MM<MAT, MATB>& mm)
	:A(mm.A),B(mm.B),alpha(1), tA(mm.tA), tB(mm.tB)
      {}
      SMM(const MAT& AA,const MATB& BB, 
	  const bool tAA, const bool tBB,
	  const SCAL a = 1)
	:A(AA), B(BB), tA(tAA), tB(tBB), alpha(a)
      {}
    };
  
  template<class SCAL, class MAT, class MATB>
    inline SMM<SCAL, MAT, MATB> 
    operator*(const SM<SCAL, MAT>& sm,const MT<MATB>& mtb) {
    return SMM<SCAL, MAT, MATB>(sm.A, mtb.A, sm.tA, mtb.tA, sm.alpha); 
  }

  template<class SCAL, class MAT, class MATB>
    inline SMM<SCAL, MAT, MATB> 
    operator*(const SM<SCAL, MAT>& sm,const MATB& BB) {
    return SMM<SCAL, MAT, MATB>(sm.A, BB, sm.tA, false, sm.alpha); 
  }
  
 template<class SCAL, class MATC, class MATA = MATC, class MATB = MATC>
   struct SMMpSM {
     const SCAL alpha;
     const MATA& A;
     const MATB& B;
     const SCAL beta;
     const MATC& C;
     const bool tA;
     const bool tB;
     SMMpSM(const MATA& AA, const MATB& BB, const MATC& CC, 
	    const bool tAA, const bool tBB,
	    const SCAL a=1, const SCAL b=1)
       :A(AA), B(BB), C(CC), alpha(a), beta(b), tA(tAA), tB(tBB)
     {} 
   };

 template<class SCAL, class MATC, class MATA, class MATB>
   inline SMMpSM<SCAL, MATC, MATA, MATB> 
   operator+(const SMM<SCAL, MATA, MATB>& smm, const SM<SCAL, MATC>& sm) {
   return SMMpSM<SCAL, MATC, MATA, MATB>
     (smm.A, smm.B, sm.A, smm.tA, smm.tB, smm.alpha, sm.alpha);
 }
#if 0

  template<class SCAL, class MATC, class MATA, class MATB>
    inline SMMpSM<SCAL, MATC, MATA, MATB> 
    operator+(const SMM<SCAL, MATA, MATB>& smm, MATC& CC) {
    return SMMpSM<SCAL, MATC, MATA, MATB>
      (smm.A, smm.B, CC, smm.tA, smm.tB, smm.alpha, 1);
  }
  template<class SCAL, class MATC, class MATA, class MATB>
    inline SMMpSM<SCAL, MATC, MATA, MATB> 
    operator+(const MM<MATA, MATB>& mm, const SM<SCAL, MATC>& sm) {
    return SMMpSM<SCAL, MATC, MATA, MATB>
      (mm.A, mm.B, sm.A, mm.tA, mm.tB, 1, sm.alpha);
  }
#endif

  
  template<class SCAL, class MAT>
    struct SM2pSM {
      const SCAL alpha;
      const MAT& A;
      const SCAL beta;
      const MAT& C;
      SM2pSM(const MAT& AA, const MAT& CC, const SCAL a = 1, const SCAL b = 0)
	: A(AA), alpha(a), C(CC), beta(b)
      {}
    };
  
  template<class SCAL, class MAT>
    inline SM2pSM<SCAL, MAT> 
    operator+(const SM2<SCAL, MAT>& sm2, const SM<SCAL, MAT> sm) {
    return SM2pSM<SCAL, MAT>(sm2.A, sm.A, sm2.alpha, sm.alpha); 
  }



  /* Done so far with new transpose */

template<class MAT>
    struct MMpM {
    const MAT& A;
    const MAT& B;
    const MAT& C;
    MMpM(const MAT& AA, const MAT& BB, const MAT& CC)
      :A(AA),B(BB),C(CC)
      {} 
  };


  template<class SCAL, class MAT>
    struct SMpSM {
    const SCAL alpha, beta;
    const MAT& A, B;
    SMpSM(const MAT& AA, const MAT& BB,
	  const SCAL scalar_a=1, const SCAL scalar_b=1)
      :A(AA), B(BB), alpha(scalar_a), beta(scalar_b)
      {}
  };
  
  template<class SCAL, class MAT>
    inline SMpSM<SCAL, MAT> 
    operator+(const SM<SCAL, MAT> sm1, const SM<SCAL, MAT> sm2 ) {
    return SMpSM<SCAL, MAT>(sm1.A, sm2.A, sm1.alpha, sm2.alpha);
  }

  /*
  template<class MAT>
    struct MpM {
    const MAT& A;
    const MAT& B;
    MpM(const MAT& AA,const MAT& BB)
      :A(AA),B(BB)
      {}
  };
  
  template<class MAT>
    inline MpM<MAT> operator+(const MAT& A, const MAT& B) { 
    return MpM<MAT>(A,B);
  }
  */  
  template<class MAT>
    struct MmM {
    const MAT& A;
    const MAT& B;
    MmM(const MAT& AA,const MAT& BB)
      :A(AA),B(BB)
      {}
  };
  
  template<class MAT>
    inline MmM<MAT> operator-(const MAT& A, const MAT& B) { 
    return MmM<MAT>(A,B);
  }
  

  



  /*onodig finns redan för SMM
    template<class MAT>
    inline MMpM<MAT> operator+(const MM<MAT>& mm, const MAT& CC)
    {
    return MMpM<MAT>(mm.A,mm.B,CC);
    }*/


  /*Maste ligga i arvda klassen!!*/
  /*
    Matrix::Matrix(const sMMmul& mm)
    :nrofrows(mm.A.nrofrows),nrofcols(mm.B.nrofcols)
    {
    this.multiply(mm.A,mm.B,*this,mm.tA,mm.tB,mm.alpha,0);
    }
    Matrix::Matrix(const sMMmulsMadd& mm)
    :nrofrows(mm.A.nrofrows),nrofcols(mm.B.nrofcols)
    {
    this->multiply(mm.A,mm.B,mm.C,mm.tA,mm.tB,mm.alpha,mm.beta);
    }

  */
#endif
} /* end namespace mat */
#endif
