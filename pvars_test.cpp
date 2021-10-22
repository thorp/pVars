/*******************************************************************************

pVars - The Parallel Variables Library

File: pvars_test.cpp

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
#include <cmath>
#include <random>
//#include <assert.h>

#include "pvars.h"
#include "testable.h"


using namespace std;

int main()
{
  using namespace pvarLib;

  // Init testing framework
  testable testSuite("pVars");

  
  // Test creating shapes.
  testSuite.test("Shape Creation");

  shape v100(100);
  shape* vector100Shape; //= new shape(100);
  vector100Shape= &v100;
  testSuite.verify( vector100Shape, "shape vector100Shape not created");

  shape* vectorShape= new shape(65536*1024);
  testSuite.verify( vectorShape, "shape vectorShape not created");

  
  
  // Test Initialization
  testSuite.test("pVar Initialization Tests");

  {
    pvar<double> addScalar(vectorShape);
    testSuite.verify( addScalar[5]==0.0,
		      "pVar addScalar not initialized correctly.");


    pvar<double> constAssignment(vectorShape, 8.0);
    testSuite.verify( constAssignment[5]==8.0,
		      "pVar constAssignment not initialized correctly.");

  
    // Initialize to dimensional index
    pvar<double> p_id(vectorShape, dimensionName::x);
    testSuite.verify( p_id[5]==5.0,
		      "pVar p_id not initialized correctly.");
  } 


  
  // Index Assignment
  testSuite.test("Index Assignment");

  {
    pvar<double> testIA_var1(vectorShape, 10.0);
    testSuite.verify( testIA_var1[5]==10.0,
		      "pVar testIA_var1 not initialized correctly.");

    testIA_var1.addIndex(dimensionName::x);
    testSuite.verify( testIA_var1[6]==16.0,
		      "pVar testIA_var1, bad value.");
  }


  // pvar ops
  testSuite.test("Basic pVar Operations");

  {
    // - Constant assignment
    pvar<double> testBPA_var1(vectorShape, 10.0);
    testSuite.verify( testBPA_var1[5]==10.0,
		      "pVar testBPA_var1 not initialized correctly.");
    testBPA_var1= 1.0;
    testSuite.verify( testBPA_var1[5]==1.0,
		      "testBPA_var1 was not set to constant correctly.");


    // - Add constant

    pvar<double> p_id(vectorShape, dimensionName::x);
    testSuite.verify( p_id[5]==5.0,
		      "pVar p_id not initialized correctly.");

    pvar<double> testBPA_var2(vectorShape);
    testSuite.verify( testBPA_var2[5]==0.0,
		      "pVar testBPA_var2 not initialized correctly.");

    testBPA_var2= p_id +1.0;
    testSuite.verify( testBPA_var2[5]==6.0,
		      "pVar testBPA_var2 incorrect result from add const.");


    // - Add pVar
    pvar<double> p_add(vectorShape);
    testSuite.verify( p_add[5]==0.0,
		      "pVar p_add not initialized correctly.");

    p_add= p_id + p_id;
    testSuite.verify( p_add[5]==10.0,
		      "pVar p_add incorrect result from pVar addition.");


    // - Subtract pVar
    pvar<double> p_subtr(vectorShape);
    testSuite.verify( p_subtr[5]==0.0,
		      "pVar p_subtr not initialized correctly.");

    p_subtr= p_id - p_id;
    testSuite.verify( p_subtr[5]==0.0,
		      "pVar p_subtr incorrect result from pVar subtraction.");

    
    // - Multiply pVar
    pvar<double> p_mult(vectorShape);
    testSuite.verify( p_mult[5]==0.0,
		      "pVar p_mult not initialized correctly.");

    p_mult= p_id * p_id;
    testSuite.verify( p_mult[5]==25.0,
		      "pVar p_mult incorrect result from pVar multiplcation.");
  }


  // Test Random
  testSuite.test("Random Function");

  {
    pvar<float> p_maxTest(vectorShape);  
    testSuite.verify( p_maxTest[5]==0.0,
		      "pVar p_maxTest not initialized correctly.");

    vectorShape-> initPseudoRandom();
    p_maxTest.setToRandom();

    float testRF_var1= p_maxTest[5];

    if ( testSuite.verify ( (p_maxTest[5]>=0.0)&(p_maxTest[5]<=1.0),
			    "Invalid random number generated.") )
      {
	std::cout <<"Invalid random number: " <<testRF_var1 <<std::endl;
      }
  }

  
  // Test reductions
  testSuite.test("Reductions");


  // Parallel reduction mult-add test
  {
    pvar<double> p_id(vectorShape, dimensionName::x); 
    double fmadd_result{0};
    fmadd_result= p_id.reduce_MultAdd(p_id);

    // Serial reduction mult-add test
    double sum=0.0;
    double if2=0.0;
    for (int i=0; i< vectorShape-> countElements(); ++i){
      if2=(double)i;
      sum+= if2*if2;
    }

    if (
	testSuite.verify ( abs((fmadd_result/sum)-1) < 1.0e-10,
				"Invalid parallel reduce_multAdd() result.") )
      {
	cout << "[Parallel result: " <<fmadd_result <<"]" <<endl;
	cout << "[  Serial result: " << sum         <<"]" <<endl <<endl;
      }
  }

  
  // Test reduce_MaxElement()
  { 
    pvar<float> testRED_var1(vectorShape, 0.2);  
    testSuite.verify( abs(testRED_var1[5]-0.2)< 1.0e-6,
		      "pVar testRED_var1 not initialized correctly.");

    testRED_var1[200]= 0.8;
    
    float maxResult{0.0};
    maxResult= testRED_var1.reduce_MaxElement();

    testSuite.verify( abs(maxResult-0.8) < 1.0e-6,
		      "pVar.reduce_MaxElement() generated incorrect result.");
  }

    
  // Test reduce_Sum()
  {
    pvar<int> p_sumTest(vectorShape, 1);
    testSuite.verify( p_sumTest[5]==1.0,
		      "pVar p_sumTest not initialized correctly.");

    int sum_result{0};
    sum_result= p_sumTest.reduce_Sum();
    
    testSuite.verify( sum_result==(vectorShape->countElements()),
		      "pVar.reduce_Sum() generated incorrect result.");
  }



  // Test function type mismatch error.
  testSuite.test("Invalid Condition");
  {
    pvar<int> p_sumTest(vectorShape, 1);
    testSuite.verify( p_sumTest[5]==1,
		      "pVar p_sumTest not initialized correctly.");

    pvar<int> p_intTest(vectorShape, 1);
    testSuite.verify( p_intTest[5]==1,
		      "pVar p_sumTest not initialized correctly.");

     cout <<"* Expect error since there is no INT implementation.\n";
     p_intTest= p_sumTest.ln();
  }


  
  // Test some of the uniary ops.
  testSuite.test("Uniary Operations");
  {
    // Create test data for uniary test  
    pvar<float> p_testData(vectorShape);
    p_testData.setToRandom();
    p_testData= (p_testData*2.0)-1.0;

#define MY_SIGNUM(X) (((0.0) < X) - (X < (0.0)))
    // Test uniary signum().
    pvar<float> p_signumTest(vectorShape);
    p_signumTest= p_testData.signum();
    testSuite.verify( p_signumTest[5] == (MY_SIGNUM(p_testData[5])),
		      "Serial vs parallel compare, signum failed.");

 
    // Test of abs() to see is abs(forInt) or fabs(forFloat) is called
    pvar<float> floatAbsTest = pvar<double>(vectorShape);
    pvar<float>floatTestData(vectorShape, dimensionName::x);
    floatTestData= floatTestData - 1000.001;
    floatAbsTest= floatTestData.abs();
    testSuite.verify( floatAbsTest[5] == fabs(floatTestData[5]),
		      "float abs test: incorrect result.");


    pvar<int>iTestData(vectorShape, dimensionName::x);
    pvar<int>intAbsTest(vectorShape);
    iTestData= iTestData - 1000;
    intAbsTest= iTestData.abs();
    testSuite.verify( intAbsTest[5] == abs(iTestData[5]),
		      "int abs test: incorrect result.");


    /// Test power(x)
    pvar<double> p_uniaryTest2(vectorShape, dimensionName::x);
    p_uniaryTest2= p_testData.power(2);
    testSuite.verify(abs(p_uniaryTest2[5]-std::pow(p_testData[5], 2)) < 1.0e-6,
		      "power() test: incorrect result.");

  
    // Test a uniary function for several values
    // Define the uniary function to be tested in _UNIFUNC
    // Note test of double to float pvar
    pvar<float> p_uniaryTest(vectorShape);

    //----------------------------------------------------------------------
    // Just change the uniary function here and it gets tested with and
    // without the active set.
    //----------------------------------------------------------------------
#define _UNIFUNC acos
#define USB(y) #y
#define US(x) USB(x)

    p_uniaryTest= p_testData._UNIFUNC();
    for (auto j=10; j<16;++j){
      testSuite.verify( abs(p_uniaryTest[j]-_UNIFUNC(p_testData[j])) < 1.0e6,
			"uniary function multi test failed.");
    }

  }

  
  // Test a float and an int the binary operation.
  testSuite.test("Binary Operations");
  {
    pvar<float> binaryTest_a(vectorShape);
    pvar<float> binaryTest_b(vectorShape);
    pvar<float> binaryTest_result(vectorShape);

    binaryTest_a.setToRandom();
    binaryTest_b.setToRandom();
  
    binaryTest_result= binaryTest_a.min(binaryTest_b);
    testSuite.verify( abs(binaryTest_result[14]-
			  (min(binaryTest_a[14],binaryTest_b[14])))
		      < 1.0e-6,
		      "Test of elementwise min failed.");


    
    pvar<int> binaryIntTest(vectorShape, dimensionName::x);
    pvar<int> binaryIntResult(vectorShape);  

    binaryIntResult = binaryIntTest.mod(11);
    testSuite.verify( binaryIntResult[14] == (binaryIntTest[14]%11),
		      "Elementwise binary mod() produced wrong result.");
  }


  // Test some of the binary ops.
  testSuite.test("Active Set");
  {
    vectorShape-> resetActiveSet();

    pvar<int> binaryIntResult(vectorShape);
    pvar<int> binaryIntTest(vectorShape, dimensionName::x);
    
    vectorShape-> ge(binaryIntTest, 2);
    testSuite.verify((binaryIntTest[14]>2)==(vectorShape-> activityStatus(14)),
		     "Shape activityStatus set incorrectly.");

  
    // Test Uniary operators with active and inactive elements of the
    // active set.
    //----------------------------------------------------------------------
    pvar<float> p_testData(vectorShape);
    p_testData.setToRandom();
    p_testData= (p_testData*2.0)-1.0;

    pvar<float> p_uniaryTest(vectorShape, -5);

    // Only works for elements >= 2
    p_uniaryTest= p_testData._UNIFUNC();

    // Test inactive elements
    testSuite.verify( abs(p_uniaryTest[1]-(-5)) < 1.0e6,
		      "inactive element uniary function multi test failed.");

    // Test active elements
    for (auto j=10; j<16;++j){
      testSuite.verify( abs(p_uniaryTest[j]-_UNIFUNC(p_testData[j])) < 1.0e6,
			"active set uniary function multi test failed.");
    }


    
    // Test shape activityStatus() fcn
    //----------------------------------------------------------------------
    vectorShape-> resetActiveSet();
    
    vectorShape-> eq(binaryIntTest, 14);

    testSuite.verify(vectorShape-> activityStatus(14),
		     "Activity status not true in eq()"); 

    testSuite.verify(!vectorShape-> activityStatus(110),
		     "Activity status not false in eq()");     



    /*
      pvar<int> binaryIntResult(vectorShape);  //remove
      pvar<int> binaryIntTest(vectorShape, dimensionName::x);//remove
      pvar<double> p_id(vectorShape, dimensionName::x); 
    */

    // Test areAllElementsActive() fcn 
    //----------------------------------------------------------------------
    vectorShape-> resetActiveSet();
    pvar<double> p_id(vectorShape, dimensionName::x); 

    testSuite.verify(vectorShape-> areAllElementsActive(),
		     "Test1: areAllElementsActive() should have returned TRUE.");     
    
    vectorShape-> gt(p_id, -10.0); // set all elements with pid > (-10)
    testSuite.verify(vectorShape-> areAllElementsActive(),
		     "Test2: areAllElementsActive() should have returned TRUE.");     


    // Test turning off some elements.
    //----------------------------------------------------------------------
    vectorShape-> resetActiveSet();
    vectorShape-> le(p_id, 10.0);

    pvar<double> constAssignment(vectorShape, 8.0);
    constAssignment= 1.0;

    float fResult{0};
    fResult= vectorShape-> subsetReduce_Sum(constAssignment);

    testSuite.verify( fResult == 11,
		     "Wrong number of items in active set.");     

    vectorShape-> gt(p_id, -10.0);
    testSuite.verify(!vectorShape-> areAllElementsActive(),
		     "areAllElementsActive() should have returned False.");     

  
    // TEST active set min() here.
    //----------------------------------------------------------------------
    fResult=0;

    // Continue with the prior active set
    //vectorShape-> resetActiveSet();
    vectorShape-> gt(p_id, 4.0);
    
    fResult= vectorShape-> subsetReduce_Min(p_id);

    testSuite.verify( fResult == 5,
		      "min() should return 5 for activeSet with 4<id<=10");
  
    testSuite.verify( vectorShape-> countActiveElements() == 6,
		      "Number of active elements should be 6.");

    testSuite.verify(!vectorShape-> areAllElementsInactive(),
		     "Some elements should be active.");
    
    testSuite.verify(!vectorShape-> areAllElementsActive(),
		     "Some elements should be inactive.");
  }

  { 
    // Test subsetReduction MultAdd on even elements test
    // Sum the squares of all odd number < 100.
    //----------------------------------------------------------------------
    vector100Shape->resetActiveSet();
    
    int imadd_result{0};
    pvar<int> testMA_var1(vector100Shape);
    pvar<int> testMA_selfAddress(vector100Shape, dimensionName::x);

    testMA_var1= testMA_selfAddress.mod(2);
    vector100Shape-> eq(testMA_var1, 1);
    
    imadd_result= vector100Shape-> subsetReduce_MultAdd(testMA_selfAddress,
							testMA_selfAddress);
    
    int sum{0};
    for (int i=1; i < 101; i=i+2) sum += (i*i);
      
    testSuite.verify( imadd_result == sum,
		      "Result of subsetReduce_MultAdd on vector100Shape wrong.");
  }


  
  testSuite.test("Spread Communication");
  {
    // --TEST Flood()  
    pvar<int> spreadTest1dData(vectorShape);
    pvar<int> spreadTest1d(vectorShape);
    int fromElementNumber{5};
    string msg{""};

    vectorShape-> resetActiveSet();
    spreadTest1dData.setToRandom();
    spreadTest1d= 
      spreadTest1dData.flood(fromElementNumber, dimensionName::x );

    for (int i=210; i<281; i=i+7)
      {
	msg= "Flooded data is incorrect in element:"
	  +std::to_string(i) +" , (" +std::to_string(spreadTest1d[i]) +")";
	testSuite.verify( spreadTest1d[fromElementNumber] == spreadTest1d[i],
			  msg);
      }
  }

  
  // --TEST scanAdd() and shiftWithWrap()  
#define _SHORT_VECTOR_LEN 11
  shape* shortVectorShape= new shape(_SHORT_VECTOR_LEN);
  
  {
    pvar<int> scanTest1dData(shortVectorShape);
    pvar<int> scanTest1d(shortVectorShape);

    scanTest1dData= 1;
    scanTest1d= 
      scanTest1dData.scanAdd(dimensionName::x);

    scanTest1d= scanTest1d.shiftWithWrap( 1, dimensionName::x );

    string msg{""};
    for (auto i=0; i<(_SHORT_VECTOR_LEN); i++){
      msg= "Scanned data is incorrect in element:"
	+std::to_string(i) +" , (" +std::to_string(scanTest1d[i]) +")";
      testSuite.verify( scanTest1d[i] == ((i==0)?11:i), msg);
    }
  }
  

  // Enumeration test
  {
    // -- TEST Enumerate even elements.
    pvar<int> enumTest(shortVectorShape);
    pvar<int> evenElements(shortVectorShape);

    // Turn on only the even elements
    evenElements.setToIndex(dimensionName::x);
    evenElements= evenElements.mod(2);
    shortVectorShape-> eq(evenElements, 0);
  
    enumTest.enumerateActiveElements();
    
    string msg{""};
    for (auto i=0; i < ( shortVectorShape-> countElements() ); ++i){
      msg= "Enumeration " +std::to_string(i) +"is incorrect. Value: "
	+std::to_string(enumTest[i]);

      if (shortVectorShape-> activityStatus(i))
	{
	  testSuite.verify( enumTest[i] == i/2, msg);
	}
      else
	{
	  testSuite.verify( enumTest[i] == 0, msg);
	}
    }
  }
  

  return 0;
}
