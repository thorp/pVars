/*******************************************************************************

pVars - The Parallel Variables Library

File: testable.cpp

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
#include <string>

#include "testable.h"

const int _lineLength{80};

testable::testable(std::string testSuiteNameArg)
{
  testSuiteName = testSuiteNameArg;
  std::cout <<"\n<< " <<testSuiteName <<" Test Suite >>\n\n";
}


testable::~testable()
{
  if (!currentTest.empty())
    {
      std::cout <<"\n[Passed " <<verificationsPassed
		<<"/" <<verifications <<" verifications]\n\n";

      if ( verifications == verificationsPassed)
	testsPassed++;

      std::cout <<"Summary [" <<testSuiteName <<"]: " <<testsPassed
		<<" out of " <<testsTotal <<" tests passed.\n";
    }
}



void testable::test(std::string testNameArg)
{
  if (!currentTest.empty())
    {
      std::cout <<"\n[Passed " <<verificationsPassed <<"/"
		<<verifications <<" verifications]\n\n";
      if ( verifications == verificationsPassed)
	testsPassed++;
    }
  
  currentTest = testNameArg;
  
  std::cout <<"===== " <<testNameArg <<" =====\nTesting: ";
  starCount= 9; // Set to the number of characters in "Testing"
  verifications=0;
  verificationsPassed=0;
  testsTotal++;
}


bool testable::verify( bool x, std::string msg)
{
  verifications++;
  
  if (!x)
    { 
      if (verifications > 1)
	std::cout <<"\n";
  
      std::cout <<">>> Verification " <<verifications
		<<" FAILED: " <<msg <<"\n";
    }
  else{
    if (starCount > _lineLength)
      {
	std::cout <<"\n*";
	starCount= 1;
      }
    else
      {
	std::cout <<"*";
	starCount++;
      }
    verificationsPassed++;
  }

  return !x;
}
