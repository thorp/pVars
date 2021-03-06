
/*******************************************************************************

pVars - The Parallel Variables Library

File: pvars.h

Copyright 2021, John Thorp

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.

Ref: MIT License - https://opensource.org/licenses/MIT

******************************************************************************/

#include <iostream> 
#include <new>
#include <random>
#include <cmath>
#include <climits>
#include <ctime>
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include <execution>
#include <array>

#define _PVARS_MACRO_STRING(y) #y

namespace pvarLib{
  
  enum class shapeName {none, array_1d, array_2d, array_3d, square, cube};
  enum class dimensionName { all, x, y, z };
  enum class opCode { add, mult };

  
  // Forward declaration of class pvar
  template<typename T> class pvar;

  //----------------------------------------------------------------------
  // Declare relational test classes for use in the
  //   void GenericPvarTest_argPS(pvar<T>& a, P funcTest_argPS )
  // policy template.

#define _COMP_SCALAR_PREAMBLE(FUNC,EQN)				\
  template<typename T>						\
  class internal_argPS_##FUNC { T compVal; public:		\
    internal_argPS_##FUNC(T& sv) :compVal{sv} { }		\
    inline bool operator()(T& x) { return ( EQN ); }  }
  
  _COMP_SCALAR_PREAMBLE(eq, (x == compVal) );
  _COMP_SCALAR_PREAMBLE(ne, (x != compVal) );
  _COMP_SCALAR_PREAMBLE(gt, (x > compVal) );
  _COMP_SCALAR_PREAMBLE(ge, (x >= compVal) );
  _COMP_SCALAR_PREAMBLE(lt, (x < compVal) );
  _COMP_SCALAR_PREAMBLE(le, (x <= compVal) );


  //----------------------------------------------------------------------
  // Declare relational test classes for use in the
  //   void GenericPvarTest_argPP(pvar<T>& a, pvar<T>b, P funcTest_argPP )
  // policy template.
  
#define _COMP_PVAR_PREAMBLE(FUNC,EQN)				\
  template<typename T>						\
  class internal_argPP_##FUNC { public:				\
    internal_argPP_##FUNC() { }					\
    inline bool operator()(T& x, T& y) { return ( EQN ); } }

  _COMP_PVAR_PREAMBLE(eq, (x == y) );
  _COMP_PVAR_PREAMBLE(ne, (x != y) );
  _COMP_PVAR_PREAMBLE(gt, (x > y) );
  _COMP_PVAR_PREAMBLE(ge, (x >= y) );
  _COMP_PVAR_PREAMBLE(lt, (x < y) );
  _COMP_PVAR_PREAMBLE(le, (x <= y) );


  //----------------------------------------------------------------------
  // Declare uniary function classes for use in the
  //   void GenericPvarFunc_Uniary( P funcTest )
  // policy template.
 
#define UNIARY_PREAMBLE(FUNC,EQN)\
template<typename T>\
class internal_uniary_##FUNC {public:\
internal_uniary_##FUNC(){}\
inline T operator()(T& x) { return (T) EQN; }\
inline std::string name(){return #FUNC;} }

  UNIARY_PREAMBLE(sin, std::sin(x));
  UNIARY_PREAMBLE(cos, std::cos(x));
  UNIARY_PREAMBLE(tan, std::tan(x));
  UNIARY_PREAMBLE(asin, std::asin(x));
  UNIARY_PREAMBLE(acos, std::acos(x));
  UNIARY_PREAMBLE(atan, std::atan(x));

  UNIARY_PREAMBLE(sinh, std::sinh(x));
  UNIARY_PREAMBLE(cosh, std::cosh(x));
  UNIARY_PREAMBLE(tanh, std::tanh(x));
  UNIARY_PREAMBLE(asinh, std::asinh(x));
  UNIARY_PREAMBLE(acosh, std::acosh(x));
  UNIARY_PREAMBLE(atanh, std::atanh(x));

  UNIARY_PREAMBLE(sigmoid, (1/(1+std::exp(-x))) );
  UNIARY_PREAMBLE( signum, ((T(0) < x) - (x < T(0))) );
  UNIARY_PREAMBLE(     ln, (std::log(x)) );
  UNIARY_PREAMBLE(   log2, (std::log2(x)) );
  UNIARY_PREAMBLE(  log10, (std::log10(x)) );
  UNIARY_PREAMBLE(    exp, (std::exp(x)) );
  UNIARY_PREAMBLE(   sqrt, (std::sqrt(x)) );
  UNIARY_PREAMBLE(   cbrt, (std::cbrt(x)) );
  UNIARY_PREAMBLE(   ceil, (std::ceil(x)) );
  UNIARY_PREAMBLE(  floor, (std::floor(x)) );
  UNIARY_PREAMBLE(  trunc, (std::trunc(x)) );
  UNIARY_PREAMBLE(  round, (std::round(x)) );
  // UNIARY_PREAMBLE(    abs, (std::abs(x)) ); Requires special case, see impl.
  UNIARY_PREAMBLE(    erf, (std::erf(x)) );
  UNIARY_PREAMBLE(   erfc, (std::erfc(x)) );

		  
  //----------------------------------------------------------------------
  /// Class Shape: defines the shape of the calculation.
  ///
  /// Supported shapes: none, array_1d, array_2d, array_3d, square, cube
  /// All shapes are listed in the enum **shapeName**.
  //----------------------------------------------------------------------
  
  class shape{
  private:
    shapeName _form{shapeName::none};
    int _nDimensions{0};
    int* _dimensionSizes{nullptr};
    int _containingPowerOfTwo{0}; // size of array that is >= to current size.
    int _powerOfTwoDoublings{0}; // Number of doublings to reach powerOfTwoSize.
  
    // Is this element participating in the operation?
    uint* _activeSet{nullptr};
    bool _allElementsActive{true};

    //----------------------------------------------------------------------
    // Related to random handling
    //----------------------------------------------------------------------

    /// The actual max size for int pvar Rands  
    const int _halfMaxRandomInt{RAND_MAX/2};

    /// Set up internal variables related to the power of two 
    /// 
    /// Private. that are; the nearest larger power of two and
    /// and the log2 of that size
    inline void setPowerOfTwoVars(){
      _containingPowerOfTwo=1;
      _powerOfTwoDoublings=1;
      while ( _containingPowerOfTwo < countElements() ){
	_containingPowerOfTwo <<= 1;
	_powerOfTwoDoublings += 1;
      }
    }
    
  public:
    /// A vector of random ints, one per element of the shape.
    // int* _vectorOfRandomInts{nullptr};
    std::vector<int> _vectorOfRandomInts;
    
    /// Constructor- new shape from data length
    shape(int x_size);
    //shape(int x_size, int y_size);
    //shape(int x_size, int y_size, int z_size);

    /// Constructor- new shape derived from old
    shape(shape* s);

    /// Count the number of elements in the pvar.
    inline int countElements() const{
      assert ( this );
      //assert( shapeNameIs( shapeName::array_1d ) );
      return _dimensionSizes[0];
    }

    /// Count the number of dimensions in the pvar.
    inline int countDimensions() const{
      assert(shapeNameIs( shapeName::array_1d ));
      return _nDimensions;
    }

    /// Count the number of elements in a dimension in the pvar.
    inline int countElementsInDimension( dimensionName n ) const{
      assert(shapeNameIs( shapeName::array_1d ));
      return _dimensionSizes[0];
    }

    /// Return the name (type shapeName) of the shape type. 
    inline shapeName getShapeName() const{
      return _form;
    }
  
    inline bool shapeNameIs( shapeName n ) const{
      return (_form == n);
    }

    inline bool isRectilinear(){
      bool b=false;
      b |= ( _form == shapeName::array_1d );
      b |= ( _form == shapeName::array_2d );
      b |= ( _form == shapeName::array_3d );
      b |= ( _form == shapeName::square );
      b |= ( _form == shapeName::cube );
  
      return b;
    }

    // Compare shapes.
    inline bool isTheSameAs( shape* s ){
      bool b=true;
      b &= ( _form == s-> _form );
      b &= ( countDimensions() == s-> countDimensions() ) ;
      b &= ( _dimensionSizes == s-> _dimensionSizes );
  
      return b;
    }


    /// Get pointer to self
    inline shape* getPointerToSelf() { return this;}
    
    // **********************************************************************
    // * Related to Power of two data sets
    // **********************************************************************

    /// Checks if data is a power of two.
    inline bool isPowerOfTwo() const{
      return _containingPowerOfTwo==countElements();
    }

    /// Returns the smallest power of two that is >= the data set size.
    inline int smallestContainingPowerOfTwo() const{
      return _containingPowerOfTwo;
    }

    /// Returns the number of doublings to reach the containing power of two.
    inline int log2Steps() const{
      return _powerOfTwoDoublings;
    }


    // **********************************************************************
    // * Related to Random ints
    // **********************************************************************

    
    /// Sets of the pseudorandom number generator
    ///
    /// Should be called before use of the random() function.
    void initPseudoRandom();

    /// Gets pointer to array of RandomInts 
    //    int* getRandomInts();

    /// Updates the array with ( r(i) + (0.001)Rand_Max ) % (0.5 * rand_max)
    void pseudorandomize();


    // **********************************************************************
    // Related to the Active Set
    // **********************************************************************

    /// Will this element participate in parallel computations?
    inline uint activityStatus(int i){
      return _activeSet[i];
    }
    
    /// Retrieve a pointer to the active set.
    inline uint* activeSet(){
      return &_activeSet[0];
    }
    
    /// Turn on each element of the active set.
    /// Changes the indicator that all elements are participating
    /// in the computation to true then sets all of the elements to true.

    inline void resetActiveSet(){
      _allElementsActive= true;
      
#pragma acc kernels loop independent
      for (int i=0; i!=countElements(); ++i)
	_activeSet[i]= (uint) 1;
    }
    
    /// Changes the indicator that all elements are participating
    /// in the computation to false.
    inline void touchActiveSet(){ _allElementsActive= false;}

    /// Are all elements participating in the computation?
    inline bool areAllElementsActive()
    {
      uint* active= _activeSet;
      uint result{1};

      if (!_allElementsActive){
	// ActiveSet may have been touched but all may still be active.
	
#pragma acc parallel loop reduction(&:result)
	for (int i=0; i!=countElements(); ++i){
	  result &= (active[i]);
	}
      }

      return result;
    }

    /// returns the number of active elements in the active set.
    int countActiveElements()
    {
      uint* active= _activeSet;
      int result{0};

#pragma acc parallel loop reduction(+:result)
      for (int i=0; i!=countElements(); ++i){
	result += (int) (active[i]);
      }
    
      return(result);
    }

    /// returns true if there are no elements in the active set.
    bool areAllElementsInactive()
    {
      uint* active= _activeSet;
      uint result{0};

#pragma acc parallel loop reduction(|:result)
      for (int i=0; i!=countElements(); ++i){
	result |= (active[i]);
      }
    
      return(result==0);
    }


    

    // Policy object for testing the relationship of pvars to scalar
    template<typename T, typename P>
    void GenericPvarTest_argPS(pvar<T>& a, P funcTest_argPS ){
      T* tmpA{a.getLinearizedData()};
      _allElementsActive= false;
      
#pragma acc kernels loop independent
      for (int i=0; i!=countElements(); ++i)
	_activeSet[i] &= funcTest_argPS(tmpA[i]);
    }

    // Policy object for testing the relationship of pvars to pvar
    template<typename T, typename P>
    void GenericPvarTest_argPP(pvar<T>& a, pvar<T>b, P funcTest_argPP ){
      T* tmpA{a.getLinearizedData()};
      T* tmpB{b.getLinearizedData()};
      _allElementsActive= false;

#pragma acc kernels loop independent
      for (int i=0; i!=countElements(); ++i)
	_activeSet[i] &= funcTest_argPP(tmpA[i], tmpB[i]);
    }
    
    /// Test for equality to pvar
    template<typename T>
    void eq(pvar<T>& a, pvar<T>& b);

    /// Test for equality to scalar
    template<typename T>
    void eq(pvar<T>& a, T bScalar );

    /// Test for equality to 0
    template<typename T>
    void eq(pvar<T>& a);

    /// Test for inequality to pvar
    template<typename T>
    void ne(pvar<T>& a, pvar<T>& b);

    /// Test for inequality to scalar
    template<typename T>
    void ne(pvar<T>& a, T bScalar );

    /// Test for inequality to 0
    template<typename T>
    void ne(pvar<T>& a);

    /// Test for pvarA[i] > pvarB[i]
    template<typename T>
    void gt(pvar<T>& a, pvar<T>& b);

    /// Test for x > scalar
    template<typename T>
    void gt(pvar<T>& a, T bScalar );

    /// Test for x > 0
    template<typename T>
    void gt(pvar<T>& a);

    /// Test for pvarA[i] >= pvarB[i]
    template<typename T>
    void ge(pvar<T>& a, pvar<T>& b);

    /// Test for x >= scalar
    template<typename T>
    void ge(pvar<T>& a, T bScalar );

    /// Test for x >= 0
    template<typename T>
    void ge(pvar<T>& a);

    /// Test for pvarA[i] < pvarB[i]
    template<typename T>
    void lt(pvar<T>& a, pvar<T>& b);

    /// Test for x < scalar
    template<typename T>
    void lt(pvar<T>& a, T bScalar );

    /// Test for x < 0
    template<typename T>
    void lt(pvar<T>& a);

    /// Test for pvarA[i] <= pvarB[i]
    template<typename T>
    void le(pvar<T>& a, pvar<T>& b);

    /// Test for x <= scalar
    template<typename T>
    void le(pvar<T>& a, T bScalar );

    /// Test for x <= 0
    template<typename T>
    void le(pvar<T>& a);

    
    // **********************************************************************
    // Reductions using the Active Set
    // **********************************************************************

    /// Sums the results of element-wise multiplication of pvars
    ///
    /// Fused multiply add on active set of indexes.
    /// Performs a pairwise multiply on active set of indexes
    /// then returns the sum.
    /// May use fused multiply add for better performance. 
    /// @param a The second pvar of the multiplication pair.
    /// @return Scalar value same type as pvar data.
    template<typename T>
    T subsetReduce_MultAdd(pvar<T>& a, pvar<T>& b);

    /// Sums the pvar elements in the active set.
    template<typename T>
    T subsetReduce_Sum(pvar<T>& a);

    /// Finds the min() the pvar elements in the active set.
    template<typename T>
    T subsetReduce_Min(pvar<T>& a);

    /// Finds the max() the pvar elements in the active set.
    template<typename T>
    T subsetReduce_Max(pvar<T>& a);

    /// Finds the max() the pvar elements in the active set.
    template<typename T>
    T subsetReduce_Or(pvar<T>& a);

    // End of shape class
  };


  /*
    **********************************************************************
    **********************************************************************
    pvar
    **********************************************************************
    **********************************************************************
  */

  template<typename T>
  class pvar{
  private:
    T* _v;
    int length; // const?
    shape* myShape;

    int inline _snarkleSrcAddr(int idx, int st, int stPrev) const{
      return (st*(idx/stPrev)) + (stPrev-1);
    }

    int inline _snarkleDestAddr(int idx, int st, int stPrev) const{
      return (st*(idx/stPrev)) + stPrev+ (idx%stPrev);
    }

    // Internal use only. not presently part of API.
    pvar<T>(shape *s)
      : _v{new T[s-> countElements()]},
	_pvec{std::vector<T>(s->countElements(),(T) 0)},
	length{s-> countElements()},
	myShape{s}
    {
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= (T) 0;
      }

    }
    

  public:
    std::vector<T> _pvec;

    /*!
     *  Creates a pvar that conforms to the specified shape and type.
     *  
     *  This constructor will make a pvar of the specified shape.
     *  It's operations will work seamlessly with other pvars of the 
     *  same shape. The new pvar's data will be initiallized to (T)0.
     */
    pvar<T>(shape& s)
      : _v{new T[s.countElements()]},
	_pvec{std::vector<T>(s.countElements(),(T) 0)},
	length{s.countElements()},
	myShape{s.getPointerToSelf()}
	//myShape{s}
    {

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= (T) 0;
      }
    }


    /*!
     *  Creates a pvar that conforms to the specified shape and type.
     *  
     *  This constructor will make a pvar of the specified shape.
     *  It's operations will work seamlessly with other pvars of the 
     *  same shape. The new pvar's data will be initiallized to the
     *  scalar value (T)initVal.
     */
    template<typename U>
    pvar<T>(shape& s, U initVal)
      : _v{new T[s.countElements()]},
	_pvec{std::vector<T>(s.countElements(),(T) initVal)},
	length{s.countElements()},
	myShape{s.getPointerToSelf()}
    {
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= initVal;
      }
    }
      

    /*!
     *  Creates a pvar that conforms to the specified shape and type.
     *  
     *  This constructor will make a pvar of the specified shape.
     *  It's operations will work seamlessly with other pvars of the 
     *  same shape. The new pvar's data will be initiallized to the
     *  axis index in dimension n.
     */
    pvar(shape& s, dimensionName n)
      : _v{new T[s.countElements()]},
	_pvec{std::vector<T>(s.countElements(),(T) 0)},
	length  {s.countElements()},
	myShape{s.getPointerToSelf()}
    {

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= (T)i;
      }

      std::vector<T> _tmpvec{std::vector<T>(s.countElements(),(T) 1)};
      std::exclusive_scan(std::execution::par_unseq,
			  _tmpvec.begin(),
			  _tmpvec.end(),
			  _pvec.begin(), (T) 0);
    }


    /*!
     *  Destroys this pvar removing all associated data.
     *  
     *  This action conforms with RAII principles and cleans up
     *  internally heap allocated data and pointers.
     */
    inline ~pvar() {delete[] _v;}


    /*!
     *  Creates a pvar from another pvar.
     *  
     *  This constructor will make a pvar based on the pvar argument.
     *  It's operations will work seamlessly with other pvars of the 
     *  same shape. The new pvar's data will be initiallized to the
     *  the values in the original pvar.
     *
     *  The const designator indicates that the copy constructor will
     *  not change the original.
     */
    template <typename U>
    pvar<T>( pvar<U>& a )
      : _v{new T[a.length]},
	_pvec{std::vector<T>( a.size(), (T) 0 )},
	myShape{a.getShape()},
	length{a.size()}
    {
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= (T) a._v[i];
      }

      std::transform( std::execution::par_unseq,
		      a._pvec.begin(),
		      a._pvec.end(),
		      _pvec.begin(),
		      [](U tmp){return (T)tmp;} );
    }
    

    // Copy assignment - const
    template <typename U>
    pvar<T>& operator=(const pvar<U>& a )
    {
      T* p = new T[length];

#pragma acc kernels loop independent
      for (int i=0; i!=a.length; ++i)
	p[i]= (T) a._v[i];
  
      delete[] _v;
      _v= p;
      length= a.length;
  
      return(*this);
    }


    // Copy assignment - not const
    template <typename U>
    pvar<T>& operator=(pvar<U>& a )
    {
      assert(myShape);
      assert(a.getShape());
      assert( myShape->countElements() == a.getShape()->countElements() );

      T* p = new T[length];

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	p[i]= (T) a[i];
  
      delete[] _v;
      _v= p;

      std::transform( std::execution::par_unseq,
		      a._pvec.begin(),
		      a._pvec.end(),
		      _pvec.begin(),
		      [](U tmp){return (T)tmp;} );

      return(*this);
    }

    
    // Constant assignment
    pvar& operator=(const T a ){

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	_v[i]= a;

      std::for_each(std::execution::par_unseq,
	       _pvec.begin(),
	       _pvec.end(),
	       [a](T &n){ n= a; });
  
      return(*this);
    }

    
    /*!
     *  Creates a pvar by moving the data pointer from the old pvar.
     *  
     *  This move constructor will make a pvar from the original
     *  pvar. It is used in calculations on the RHS of an assignment.
     *  (Questions? Research "C++ Move Constructors")
     */
    pvar(pvar&& a)
      :_v{a._v}, length{a.length}
    {
      std::cout <<"move cons rhs start";

      a._v= nullptr;
      a.length= 0;
      
    }

    
    /*!
     *  Assigns a pvar to the LHS by moving the data pointer from the RHS pvar.
     *  
     *  This move constructor will make a pvar from the original
     *  pvar. It is used in calculations on the RHS of an assignment.
     *  (Questions? Research "C++ Move Assignment")
     */
    pvar& operator=(pvar&& a)
    {
      std::cout <<"move rhs start";
      _v= std::move(a._v);
      length=std::move(a.length);
      
      a._v= nullptr;
      a.length= 0;
      std::cout <<"move rhs start";
	
      return (*this);
    }


    // Returns single element
    T& operator[](int index){return _v[index];}
    T& operator[](int index) const {return _v[index];}

    // Returns the number of elements. (redundant with shape.countElements()).lements
    int size() {return (length);}

    // Returns the shape of this pvar.
    shape* getShape() const{return myShape;}

    // Returns pointer to data_vector.
    inline T* getLinearizedData() const { return &_v[0]; }

    /*!
     *  Gets the value at a specific index in the pvar.
     *  
     *  Included temporarily for testing _pvec.
     */
    inline T getValueAt(int i) const { return _pvec.at(i); }



    /*
      ----------------------------------------------------------------------
      vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      Parallel Operations
      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ----------------------------------------------------------------------
    */
    

    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Index Commands
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */

    /*!
     *  Set each element in the pVar to the element's index.
     *  
     *  Details...
     */
    void setToIndex(dimensionName dim){
      assert ( dim == dimensionName::x);
 
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	_v[i]= (T)i;
      }

      std::vector<T> _tmpvec{std::vector<T>( length, (T) 1 )};
      std::exclusive_scan(std::execution::par_unseq,
			  _tmpvec.begin(),
			  _tmpvec.end(),
			  _pvec.begin(), (T) 0);
    }


    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Overloaded Operator Commands
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */


    /// Addition operator (+): add constant to pvar
    ///
    /// Add a constant number to each pvar element. C++ type coercion
    /// handles type mismatches.
    pvar& operator+(T x)
    {
      // Make a new pvar to hold the RHS result.
      pvar* plv{ new pvar<T>(myShape) };

      T* tmpV{plv -> _v};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] + x;


      std::transform( std::execution::par_unseq,
		      _pvec.begin(),       // Src First
		      _pvec.end(),         // Src Last
		      plv-> _pvec.begin(), //Dest First
		      [x](T tmp){return tmp + x;} );
      
      return(*plv);
    }


    
    /// Subtraction operator (-): subtract constant from pvar
    ///
    /// Subtract a constant number to each pvar element. C++ type coercion
    /// handles type mismatches.
    pvar& operator-(T x)
    {
      //make a new pvar.
      pvar* plv = new pvar(this-> getShape());
      plv-> length= length;
      T* tmpV{plv -> _v};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] - x;
  
      return(*plv);
    }

    
    /// Multiplicion operator (*): Multiply pvar by constant
    ///
    /// Multiply a constant number with each pvar element. C++ type coercion
    /// handles type mismatches.
    pvar& operator*(T x)
    {
      //make a new pvar.
      pvar* plv = new pvar(this-> getShape());
      plv-> length= length;
      T* tmpV{plv -> _v};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] * x;
  
      return(*plv);
    }


    /// Division operator (/): Divide pvar by constant
    ///
    /// Divide each pvar element by a constant number.
    /// C++ type coercion
    /// handles type mismatches.
    template <typename U>
    pvar<T>& operator/(U x)
    {
      //make a new pvar.
      pvar* plv = new pvar(this-> getShape());
      T* tmpV{plv -> getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!= length; ++i)
	tmpV[i]= _v[i] / x;
  
      return(*plv);
    }
 
    
    /// Addition operator (+): add pvar to pvar
    ///
    /// For each element,
    /// add the corresponding pvar data elements. C++ type coercion
    /// handles type mismatches.
    template <typename U>
    pvar<T>& operator+(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(a);
      T* tmpV{plv -> _v};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] + (T) tmpW[i];

      std::transform( std::execution::par_unseq,
		      _pvec.begin(),       // Src1 First
		      _pvec.end(),         // Src1 Last
		      a._pvec.begin(),     // Src2 First
		      plv-> _pvec.begin(), // Dest First
		      [](T tmp1, U tmp2){return tmp1 + ((T)tmp2);} );
      
      return(*plv);
    }


    /// Subtraction operator (-): subtract pvar to pvar
    ///
    /// For each element,
    /// subtract the corresponding pvar data elements. C++ type coercion
    /// handles type mismatches.
    template <typename U>
    pvar<T>& operator-(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(a);
      plv-> length= length;  
      T* tmpV{plv -> _v};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] - (T)tmpW[i];
  
      std::transform( std::execution::par_unseq,
		      _pvec.begin(),       // Src1 First
		      _pvec.end(),         // Src1 Last
		      a._pvec.begin(),     // Src2 First
		      plv-> _pvec.begin(), // Dest First
		      [](T tmp1, U tmp2){return tmp1 - ((T)tmp2);} );

      return(*plv);
    }


    /// Multiplication operator (*): Multiply pvar with pvar
    ///
    /// For each element,
    /// multiply the corresponding pvar data elements. C++ type coercion
    /// handles type mismatches. 
    template <typename U>
    pvar<T>& operator*(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(a);
      plv-> length= length;
      T* tmpV{plv -> _v};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] * (T) tmpW[i];

      std::transform( std::execution::par_unseq,
		      _pvec.begin(),       // Src1 First
		      _pvec.end(),         // Src1 Last
		      a._pvec.begin(),     // Src2 First
		      plv-> _pvec.begin(), // Dest First
		      [](T tmp1, U tmp2){return tmp1 * ((T)tmp2);} );
    
      return(*plv);
    }


    /// Division operator (*): Divide pvar by pvar
    ///
    /// For each element,
    /// divide the corresponding pvar data elements. C++ type coercion
    /// handles type mismatches. 
    template <typename U>
    pvar<T>& operator/(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(a);
      T* tmpV{plv -> _v};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[i] / tmpW[i];

      std::transform( std::execution::par_unseq,
		      _pvec.begin(),       // Src1 First
		      _pvec.end(),         // Src1 Last
		      a._pvec.begin(),     // Src2 First
		      plv-> _pvec.begin(), // Dest First
		      [](T tmp1, U tmp2){return tmp1 / ((T)tmp2);} );

      return(*plv);
    }



    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REDUCE Commands
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */

    /// Sums the data elements in 'this'.
    ///
    /// Returns the sum of the data elements. All indexes participate.
    /// @param a The pvar to be reduced.
    /// @return Scalar value same type as pvar data.
    T reduce_Sum(){

      // Make a new T.
      //T* result = new T{0.0};
      T result{0.0};

#pragma acc parallel loop reduction(+:result)
      for (int i=0; i != length; ++i)
	result += _v[i];
    
      return(result);
    } 
    

    /// Finds the Max of the data elements in 'this'.
    ///
    /// Returns the maximum value of the data elements.
    /// All indexes participate.
    /// @param a The pvar to be reduced.
    /// @return Scalar value same type as pvar data.
    T reduce_MaxElement(){

      // Make a new T.
      //T* result = new T{0.0};
      T result{_v[0]};

#pragma acc parallel loop reduction(max:result)
      for (int i=0; i != length; ++i)
	result = std::max(result, _v[i]);
    
      return(result);
    }
    
    
    /// Sums the results of element-wise multiplication of pvars;
    /// 'this' and 'a'.
    ///
    /// Performs a pairwise multiply of the two pvars then returns the sum.
    /// May use fused multiply add for better performance. 
    /// @param a The pvar to be reduced.
    /// @return Scalar value same type as pvar data.
    template <typename U>
    T reduce_MultAdd(pvar<U>& a){

      T result{0.0};
      U* tmpW{a.getLinearizedData()};

      // fmadd
#pragma acc parallel loop reduction(+:result)
      for (int i=0; i != length; ++i)
	result += _v[i] * (T) tmpW[i];
    
      return(result);
    } 
    

    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Uniary Commands
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */

    /// Policy object for testing the relationship of pvars to scalar
    ///
    /// This template is used interally to unify the actions of
    /// uniary functions. Basically the pvar.uniary() definition
    /// calls this function with the correct uniary, like sin(),
    /// then it calls the apropriate internal serial function.
    ///
    /// @param typename_P is the class of the serial function.
    /// This template calls either a fast All Elements version or
    /// a subset of the elements if not all elements of the shape
    /// are participating.
    template<typename U, typename P>
    pvar<U>& GenericPvarFunc_Uniary( P funcTest_arg0 ){
      shape* currentShape{this->getShape()};
      
      //make a new pvar.
      pvar* plv = new pvar<U>(currentShape);

      if ( currentShape-> areAllElementsActive() ){
	if constexpr (std::is_floating_point<U>::value)
		       {
			 U* tmpV{plv -> _v};

#pragma acc kernels loop independent
			 for (int i=0; i!=length; ++i)
			   tmpV[i]= (U) funcTest_arg0(_v[i]);
		       }
	else
	  std::cerr <<"\n\nERROR: pvar." <<funcTest_arg0.name()
		    <<"() is only implemented for floating point pvars.\n\n";

      }
      else{
	uint* active= currentShape-> activeSet();

	if constexpr (std::is_floating_point<U>::value)
		       {
			 U* tmpV{plv -> _v};

#pragma acc kernels loop independent
			 for (int i=0; i!=length; ++i)
			   tmpV[i]= (U) ( active[i]? funcTest_arg0(_v[i]): 0 );
		       }
	else
	  std::cerr <<"\n\nERROR: pvar.function() is only implemented for floating point pvars.\n\n";
      }

      return(*plv);
    }

    // Standard math lib functions that are implemented with
    // GenericPvarFunc_Uniary()
    pvar& cbrt();
    pvar& erf();
    pvar& erfc();
    pvar& exp();
    pvar& ln();
    pvar& log2();
    pvar& log10();
    pvar& sigmoid();
    pvar& signum();
    pvar& sqrt();

    pvar& ceil();
    pvar& floor();
    pvar& round();
    pvar& trunc();

    /// Compute trig functions for each element in the activeSet. 
    ///
    /// If the active set is all on, the routines can run a different algo.
    /// Computes function(x) for each element, where x is in radians.
    pvar& sin();
    pvar& cos();
    pvar& tan();
    pvar& asin();
    pvar& acos();
    pvar& atan();
    
    //----------------------------------------------------------------------
    // Hyperbolic functions
    //----------------------------------------------------------------------
    pvar& sinh();
    pvar& cosh();
    pvar& tanh();
    pvar& asinh();
    pvar& acosh();
    pvar& atanh();


    /// ABS() has a fcn for fp and ints
    pvar<T>& abs(){
      shape* currentShape{this->getShape()};
      
      //make a new pvar.
      pvar* plv = new pvar<T>(currentShape);
      T* tmpV{plv -> _v};

      if ( currentShape-> areAllElementsActive() ){
	if constexpr (std::is_floating_point<T>::value)
		       {
#pragma acc kernels loop independent
			 for (int i=0; i!=length; ++i)
			   tmpV[i]= (T) std::fabs(_v[i]);
		       }
	else
		       {
#pragma acc kernels loop independent
			 for (int i=0; i!=length; ++i)
			   tmpV[i]= (T) std::abs(_v[i]);
		       }
      }
      else{
	uint* active= currentShape-> activeSet();

	if constexpr (std::is_floating_point<T>::value)
		       {
#pragma acc kernels loop independent
			 for (int i=0; i!=length; ++i)
			   tmpV[i]= (T) ( active[i]? std::fabs(_v[i]): 0 );
		       }
	else
	  {
#pragma acc kernels loop independent
	    for (int i=0; i!=length; ++i)
	      tmpV[i]= (T) ( active[i]? std::abs(_v[i]): 0 );
	  }
      }

      return(*plv);
    }

    
    /// Generate a pvar filled with random numbers. 
    ///
    /// For floating point, the resulting vector is between 0.0 and 1.0.
    /// For integral values, the resulting vector is between 0 and RANDOM_MAX/2.
    void setToRandom()
    {
      setToIndex( dimensionName::x );

      int r{std::rand()};
      double _pvarSize{ (double) length };
  
      if constexpr (std::is_floating_point<T>::value){
	  double _shift{ (double) r };

	  // Multiply each element id by a large random number then take
	  // remainder of division by number of elements which gives
	  // a number 0<=n<=number of elements. Thne divide by the number
	  // of elements get number in range 0<=n<=1.
	  std::for_each(std::execution::par_unseq,
			_pvec.begin(),
			_pvec.end(),
			//[r](int& n){ n= ((n^r)&0xEFFF)+1; });
			[_pvarSize, _shift](T& n){
			  n= (T) std::fmod((n*_shift), _pvarSize)/_pvarSize;
			});
	}

      if constexpr (std::is_integral<T>::value){
	  std::for_each(std::execution::par_unseq,
			_pvec.begin(),
			_pvec.end(),
			[_pvarSize,r](T& n){
			  n= (((int)(std::fmod((double)n*(double)r, (double) INT_MAX)))^r);
			});
	}

      /*
      vector<int>* random_v;
      random_v= this-> getShape()-> getRandomInts();
      
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i){
	if constexpr (std::is_integral<T>::value)
		       _v[i]= (T) random_v[i];
	
	if constexpr (std::is_floating_point<T>::value)
		       _v[i]= ((T) random_v[i])/((T) (RAND_MAX/2));
      }
      */
      /*
      if constexpr (std::is_integral<T>::value){
	  std::transform( std::execution::par_unseq,
			  myShape->_vectorOfRandomInts.begin(),
			  myShape->_vectorOfRandomInts.end(),
			  _pvec.begin(),
			  [](int tmp1){return (T) tmp1;} );
	}


      if constexpr (std::is_floating_point<T>::value){
	  std::transform( std::execution::par_unseq,
			  myShape->_vectorOfRandomInts.begin(),
			  myShape->_vectorOfRandomInts.end(),
			  _pvec.begin(),
			  [](int tmp1){
			    return (T)tmp1; //((T) tmp1)/((T) (RAND_MAX/2));
			  } );
	}	  
      

      // Pseudorandomize for the next call
      this-> getShape()-> pseudorandomize();
      */

      
      
    }


    /// Compute gompertz function for each element. 
    ///
    /// The formula for the Gompertz function is:
    /// a*e^( -b*e^(-c*t))
    /// @param a asymptote
    /// @param b x displacement
    /// @param c y scale (growth rate)
    /// The source pvar contains the time variable.
    pvar& gompertz(double a, double b, double c)
    {
      //make a new pvar.
      pvar* plv = new pvar(this->getShape());
      plv-> length= length;
      T* tmpV{plv -> _v};

      if constexpr (std::is_floating_point<T>::value)
		     {

#pragma acc kernels loop independent
		       for (int i=0; i!=length; ++i)
			 tmpV[i]= (T) (a* std::exp( -b* std::exp(-c*v[i])));
		     }
      else
	std::cerr <<"\n\nERROR: pvar.gompertz() is only implemented for floating point pvars.\n\n";
      
      return(*plv);
    }


    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Binary Commands: pvar op constant
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */

    /// Compute power(double x) for each element. 
    ///
    /// Compute power(x) for each element, where x is an integer.
    /// Each pvar element is the base which is raised to the x power.
    pvar& power( T p ) 
    {
      //make a new pvar.
      pvar* plv = new pvar(this->getShape());
      T* tmpV{plv -> getLinearizedData()};

      if constexpr (std::is_floating_point<T>::value)
		     {

#pragma acc kernels loop independent
		       for (int i=0; i!=length; ++i)
			 tmpV[i]= (T) std::pow(_v[i], p);
		     }
      else
	std::cerr <<"\n\nERROR: pvar.pow() is only implemented for floating point pvars.\n\n";
    
      return(*plv);
    }


    /// Compute mod(x) for each element. 
    ///
    /// Compute mod(x) for each element.
    pvar<T>& mod(T q)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
      //T* tmpW{a.getLinearizedData()};

      if constexpr (std::is_integral<T>::value)
		     {
#pragma acc kernels loop independent
		       for (int i=0; i!=length; ++i)
			 tmpV[i]= _v[i] % q ;
		     }
      else
	std::cerr <<"\n\nERROR: pvar.mod() is only implemented for integer pvars.\n\n";

      return(*plv);
    }


    /*
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Binary Commands: pvar op pvar
      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    */

    
    /// Compute max(x,y) for each element. 
    ///
    /// Compute max(x,y) for each element.
    template <typename U>
    pvar<T>& max(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= (T) std::max( _v[i], (T)tmpW[i] );
    
      return(*plv);
    }


    /// Compute min(x) for each element. 
    ///
    /// Compute min(x) for each element.
    template <typename U>
    pvar<T>& min(pvar<U>& a)
    {
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
      U* tmpW{a.getLinearizedData()};

#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= (T) std::min( _v[i], (T)tmpW[i] );
    
      return(*plv);
    }


    //----------------------------------------------------------------------
    // Uniary communications on a pvar
    //----------------------------------------------------------------------


    /// ShiftWithWrap()- shifts the values up or down the dimension by N positions
    ///
    /// When a value goes off one end it come back on the other.
    /// Positive distance shifts the values up the axis.
    /// Negative shifts toward element 0.
    /// @param distance how many positions to shift the values
    /// @param dimensionName the dimension to use for the shift
    /// @return P a pvar with values shifted.
    pvar<T>& shiftWithWrap(int distance, dimensionName x){
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
	
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= _v[(i+length-distance)%length];
    
      return(*plv);     
    }


    /// ShiftWithInsert()- shifts the values up or down the dimension by N positions
    ///
    /// When a value is shifted off one end of the dimension, the insertValue is added
    /// to the other end.
    /// Positive distance shifts the values up the axis.
    /// Negative shifts toward element 0.
    /// @param distance how many positions to shift the values
    /// @param dimensionName the dimension to use for the shift
    /// @param insertVal the value to be inserted in the space opened after the shift. 
    /// @return P a pvar with values shifted.
    pvar<T>& shiftWithInsert(int distance, dimensionName x, T insertVal){
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
	
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= ((0<=(i-distance))&&((i-distance)<length))?
	  _v[i-distance]: insertVal;
    
      return(*plv);     
    }
    

    /// Copies the value from the index location to the rest
    /// of the dimension.
    ///
    /// @param fromElementNumber the zero-based index number with the source value
    /// @param dimensionName the name of the demensionto use for the scan
    /// @return P a pvar with all values set to the value of the source element
    pvar<T>& flood(int fromElementNumber, dimensionName x ){
      //make a new pvar.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};
      T  fromVal{_v[fromElementNumber]};
      uint* _active{getShape()-> activeSet()};
	
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= (_active[i])? fromVal: tmpV[i];
    
      return(*plv);     
    }
    

    /// ScanAdd(): Adds the value from the current location to the RunningTotal
    /// of all lower indices in the dimension, if the location is active.
    ///
    /// Uses the Parallel Prefix Scan() Algorithm (which, in openacc,
    /// works without knowing
    /// the number of processors).
    /// Internally scales up the pvar's data to the next power of two.
    /// Scans across the
    /// active elements. The internal data is then scaled back down.
    /// Inactive elements remain untouched.

    pvar<T>& scanAdd(dimensionName x ){
      
      //make a new pvar for return.
      pvar* plv = new pvar<T>(this->getShape());
      T* tmpV{plv-> getLinearizedData()};

      uint* _active{getShape()-> activeSet()};

      assert (length < (INT_MAX/2) );

      // * Scale up to next power of two
      int powerOfTwo{this-> getShape()-> smallestContainingPowerOfTwo()};
      int logOfLength{this-> getShape()-> log2Steps()};

      T* tmpX= new T [powerOfTwo];

      // * Copy data across to power of two array from active elements.
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpX[i]= (_active[i])? _v[i]: 0;

      // Initialize Amdahl Logarithmic variables.
      int subtreeSize=1;
      int subtreeSize_previous=0;
      
      // Amdahl Logarithmic Loop - (Serial loop, log2 number of iters)
      for (int amdahl_LoopIndex=1;
	   amdahl_LoopIndex < logOfLength;
	   ++amdahl_LoopIndex)
	{
	  subtreeSize_previous= subtreeSize;
	  subtreeSize= subtreeSize*2;

	  // Thorp Parallel Expression (Parallel loop, independant iters)
#pragma acc kernels loop independent
	  for (int i=0; i!=(powerOfTwo/2); ++i){
	    tmpX[_snarkleDestAddr(i, subtreeSize, subtreeSize_previous)]
	      += tmpX[_snarkleSrcAddr(i, subtreeSize, subtreeSize_previous)];
	  }
	}

      // Copy back data from power of 2 array to returned varible.
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpV[i]= (_active[i])? tmpX[i]: tmpV[i];
      
      return(*plv);     
    }


    /// EnumerateActive(): Gives each active location a unique number
    ///
    /// Uses the Parallel scanAdd() Algo (which, in openacc,
    /// works without knowing
    /// the number of processors).
    /// Initializes each active element with 1.
    /// Internally scales up the pvar's data to the next power of two.
    /// Scans across the
    /// active elements. The internal data is then scaled back down.
    /// Inactive elements are set to 0.

    void enumerateActiveElements(){
      
      uint* _active{getShape()-> activeSet()};
      assert (length < (INT_MAX/2) );
      
      // * Scale up to next power of two
      int powerOfTwo{this-> getShape()-> smallestContainingPowerOfTwo()};
      int logOfLength{this-> getShape()-> log2Steps()};

      T* tmpX= new T [powerOfTwo];

      // * Copy data across to power of two array from active elements.
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	tmpX[i]= (_active[i])? 1: 0;

      // Initialize Amdahl Logarithmic variables.
      int subtreeSize=1;
      int subtreeSize_previous=0;
      
      // Amdahl Logarithmic Loop - (Serial loop, log2 number of iters)
      for (int amdahl_LoopIndex=1;
	   amdahl_LoopIndex < logOfLength;
	   ++amdahl_LoopIndex)
	{
	  subtreeSize_previous= subtreeSize;
	  subtreeSize= subtreeSize*2;

	  // Thorp Parallel Expression (Parallel loop, independant iters)
#pragma acc kernels loop independent
	  for (int i=0; i!=(powerOfTwo/2); ++i){
	    tmpX[_snarkleDestAddr(i, subtreeSize, subtreeSize_previous)]
	      += tmpX[_snarkleSrcAddr(i, subtreeSize, subtreeSize_previous)];
	  }
	}

      // Copy back data from power of 2 array to returned varible.
#pragma acc kernels loop independent
      for (int i=0; i!=length; ++i)
	_v[i]= (_active[i])? tmpX[i]-1: _v[i];
    }
    
    // End of pvar class
  };



  /*
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Shape Command Definitions
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */


  /*
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Shape: Active Set Operations
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

  
  /* Kept for reference. Each of the test(pvar&, pvar&) member
     functions looked like this before refactor into policy objects.
     Ref: stroustrup, a tour of c++, p. 86

  template<typename T>
  void shape::eq(pvar<T>& a, pvar<T>& b){
    T* tmpA{a.getLinearizedData()};
    T* tmpB{b.getLinearizedData()};

#pragma acc kernels loop independent
    for (int i=0; i!=countElements(); ++i)
      _activeSet[i] &= (tmpA[i] == tmpB[i]);
  }
  */
  /// Test for equality to pvar
  template<typename T>
  void shape::eq(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_eq<T>() );
  }


  
  /* Kept for reference. Each of the test(pvar&, scalarTest) member
     functions looked like this before refactor into policy objects.
     Ref: stroustrup, a tour of c++, p. 86

  template<typename T>
  void shape::eq(pvar<T>& a, T bScalar ){
    T* tmpA{a.getLinearizedData()};

#pragma acc kernels loop independent
    for (int i=0; i!=countElements(); ++i)
      _activeSet[i] &= (tmpA[i] == bScalar);
    }
  */



  /// Test for equality to scalar
  template<typename T>
  void shape::eq(pvar<T>& a, T bScalar ){
    GenericPvarTest_argPS( a, internal_argPS_eq<T>{bScalar} );
  }


  /* Kept for reference. Each of the test(pvar&, scalarTest) member
     functions looked like this before refactor into policy objects.
     Ref: stroustrup, a tour of c++, p. 86

<typename T>
  void shape::eq(pvar<T>& a){
    T* tmpA{a.getLinearizedData()};

#pragma acc kernels loop independent
    for (int i=0; i!=countElements(); ++i)
      _activeSet[i] &= ( tmpA[i] == (T) 0 );
  }
  */

  
  /// Test for equality to 0
  template<typename T>
  void shape::eq(pvar<T>& a){
    eq(a,(T)0);
  }

  
  /// Test for inequality to pvar
  template<typename T>
  void shape::ne(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_ne<T>() );
  }

  
  /// Test for inequality to scalar
  template<typename T>
  void shape::ne(pvar<T>& a, T bScalar ){
        GenericPvarTest_argPS( a, internal_argPS_ne<T>{bScalar} );
  }


  /// Test for inequality to 0
  template<typename T>
  void shape::ne(pvar<T>& a){
    ne(a,(T)0);
  }


  /// Test for x>y to pvar
  template<typename T>
  void shape::gt(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_gt<T>() );
  }


  /// Test for x > scalar: for all x
  template<typename T>
  void shape::gt(pvar<T>& a, T bScalar ){
        GenericPvarTest_argPS( a, internal_argPS_gt<T>{bScalar} );
  }


  /// Test for x>0
  template<typename T>
  void shape::gt(pvar<T>& a){
    gt(a,(T)0);
  }


  /// Test for x>=y to pvar
  template<typename T>
  void shape::ge(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_ge<T>() );
  }


  /// Test for x >= scalar: for all x
  template<typename T>
  void shape::ge(pvar<T>& a, T bScalar ){
        GenericPvarTest_argPS( a, internal_argPS_ge<T>{bScalar} );
  }


  /// Test for x>=0
  template<typename T>
  void shape::ge(pvar<T>& a){
    ge(a,(T)0);
  }


  /// Test for x<y to pvar
  template<typename T>
  void shape::lt(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_lt<T>() );
  }


  /// Test for x < scalar: for all x
  template<typename T>
  void shape::lt(pvar<T>& a, T bScalar ){
    GenericPvarTest_argPS( a, internal_argPS_lt<T>{bScalar} );
  }


  /// Test for x<0
  template<typename T>
  void shape::lt(pvar<T>& a){
    lt(a,(T)0);
  }


  /// Test for x<=y to pvar
  template<typename T>
  void shape::le(pvar<T>& a, pvar<T>& b){
    GenericPvarTest_argPP( a, b, internal_argPP_le<T>() );
  }


  /// Test for x <= scalar: for all x
  template<typename T>
  void shape::le(pvar<T>& a, T bScalar ){
        GenericPvarTest_argPS( a, internal_argPS_le<T>{bScalar} );
  }


  /// Test for x<=0
  template<typename T>
  void shape::le(pvar<T>& a){
    le(a,(T)0);
  }

  

  /// Perform a reduction using only elements in the active set

  // multAdd reduction using the active set.
  template<typename T>
  T shape::subsetReduce_MultAdd(pvar<T>& a, pvar<T>& b)
  {
    T result{0.0};

    T* tmpA{a.getLinearizedData()};
    T* tmpB{b.getLinearizedData()};
    
    uint* active= _activeSet;

#pragma acc parallel loop reduction(+:result)
    for (int i=0; i!=countElements(); ++i){
      result += (tmpA[i] * tmpB[i] * (T)active[i]);
    }
    
    return(result);
  }


  // Sum reduction using the active set.
  template<typename T>
  T shape::subsetReduce_Sum(pvar<T>& a)
  {
    T result{(T)0};
    
    T* tmpA{a.getLinearizedData()};
    uint* active= activeSet();

#pragma acc parallel loop reduction(+:result)
    for (int i=0; i<countElements(); ++i){
      result += (tmpA[i] * (T)active[i]);
    }
    
    return(result);
  }


    /// Finds the min() the pvar elements in the active set.
  template<typename T>
  T shape::subsetReduce_Min(pvar<T>& a)
  {
    T* tmpA{a.getLinearizedData()};
    uint* active= _activeSet;
    T result{std::numeric_limits<T>::max()};

#pragma acc parallel loop reduction(min:result)
    for (int i=0; i!=countElements(); ++i){
      result = (active[i])? min(result, tmpA[i]): result;
    }
    
    return(result);
  }

  
  /// Finds the min() of the pvar elements in the active set.
  template<typename T>
  T shape::subsetReduce_Max(pvar<T>& a)
  {
    T* tmpA{a.getLinearizedData()};
    uint* active= _activeSet;
    T result{std::numeric_limits<T>::lowest()};

#pragma acc parallel loop reduction(max:result)
    for (int i=0; i!=countElements(); ++i){
      result = (active[i])? max(result, tmpA[i]): result;
    }
    
    return(result);
  }


  /// Finds the min() of the pvar elements in the active set.
  template<typename T>
  T shape::subsetReduce_Or(pvar<T>& a)
  {
    T* tmpA{a.getLinearizedData()};
    uint* active= _activeSet;
    T result{0};

#pragma acc parallel loop reduction(|:result)
    for (int i=0; i!=countElements(); ++i){
      result = (active[i])? (result | tmpA[i]): result;
    }
    
    return(result);
  }
  
  
  /// This macro is the second half of defining the uniary pvar functions.
#define UNIARY_POSTAMBLE(FUNC)\
template<typename T>\
pvar<T>& pvar<T>::FUNC(){\
return GenericPvarFunc_Uniary<T>( internal_uniary_##FUNC<T>() );}

  UNIARY_POSTAMBLE(sin);
  UNIARY_POSTAMBLE(cos);
  UNIARY_POSTAMBLE(tan);
  UNIARY_POSTAMBLE(asin);
  UNIARY_POSTAMBLE(acos);
  UNIARY_POSTAMBLE(atan);

  UNIARY_POSTAMBLE(sinh);
  UNIARY_POSTAMBLE(cosh);
  UNIARY_POSTAMBLE(tanh);
  UNIARY_POSTAMBLE(asinh);
  UNIARY_POSTAMBLE(acosh);
  UNIARY_POSTAMBLE(atanh);

  UNIARY_POSTAMBLE(sigmoid);
  UNIARY_POSTAMBLE(signum);
  UNIARY_POSTAMBLE(ln);
  UNIARY_POSTAMBLE(log2);
  UNIARY_POSTAMBLE(log10);
  UNIARY_POSTAMBLE(exp);
  UNIARY_POSTAMBLE(sqrt);
  UNIARY_POSTAMBLE(cbrt);
  UNIARY_POSTAMBLE(ceil);
  UNIARY_POSTAMBLE(floor);
  UNIARY_POSTAMBLE(trunc);
  UNIARY_POSTAMBLE(round);
  //  UNIARY_POSTAMBLE(abs);
  UNIARY_POSTAMBLE(erf);
  UNIARY_POSTAMBLE(erfc);
  
  
  // End of Namespace  
}
