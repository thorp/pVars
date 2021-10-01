/*******************************************************************************

pVars - The Parallel Variables Library

File: pvars.cpp

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

#include <new>
#include "pvars.h"

using namespace pvarLib;

//********************
//********************
// Shape
//********************
//********************

shape::shape(shape* s)
  : _activeSet{new uint[s-> countElementsInDimension(dimensionName::x)]},
    _form{s-> getShapeName()},
    _nDimensions{s-> countDimensions()}
{
  assert(countDimensions() == 1);
  _dimensionSizes= new int[1]{ s-> countElementsInDimension(dimensionName::x) };

  setPowerOfTwoVars();

#pragma acc kernels loop independent
  for (int i=0; i!= countElementsInDimension(dimensionName::x); ++i){
    _activeSet[i]= 1;
  }

}

shape::shape(int x_size)
  : _activeSet{new uint[x_size]},
    _form{shapeName::array_1d},
    _nDimensions{1}
{
  assert(countDimensions() == 1);
  _dimensionSizes= new int[1]{ x_size };

  setPowerOfTwoVars();
  
#pragma acc kernels loop independent
  for (int i=0; i!= countElementsInDimension(dimensionName::x); ++i){
    _activeSet[i]= 1;
  }
}


/// Fill a pvar filled with random numbers
/// 
/// For floating point they are between 0 and 1.
/// For int they are 0 to RAND_MAX/2
///
/// Unfortunately we either have to dive into non-portable cuda
/// or generate randoms on the host and move them over.
/// For now, we have opted to generate on the host and move them. Slow :-(
void shape::initPseudoRandom(){
  // Generate random numbers on host then move over to device.
  // As recommended by Mat Colgrave-
  // 2021 ref: https://forums.developer.nvidia.com/t/
  // ... simple-way-of-generating-normally-distributed
  // ... -random-numbers-with-curand-in-parallel-loop/165603

  // Generate random ints on host.
  _vectorOfRandomInts= new int[this-> countElements()];
  std::srand(std::time(nullptr));

  for (int i=0; i != countElements(); ++i)
    _vectorOfRandomInts[i]= std::rand()/2;
}


/// Gets pointer to array of RandomInts 
int* shape::getRandomInts(){
  return &_vectorOfRandomInts[0];
}


/// Updates the array with ( r(i) + (0.001)Rand_Max ) % (0.5 * rand_max)
void shape::pseudorandomize(){
  
#pragma acc kernels loop independent
  for (int i=0; i!= countElements() ; ++i){
    _vectorOfRandomInts[i]+= (_halfMaxRandomInt/19);
    // A prime palendrome, 75557!
    _vectorOfRandomInts[i]= _vectorOfRandomInts[i] % _halfMaxRandomInt;
  }
}
