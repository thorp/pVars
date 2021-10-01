// parlate_test.cpp

#include <iostream> 
#include <new>
#include <cmath>
#include <random>
#include "pvars.h"

using namespace std;

int main()
{
  using namespace pvarLib;

  // Declare variables.
  shape* vectorShape= new shape(65536*1024);


  // Initialization
  cout <<"\n===== Initialization Tests =====\n";
  pvar<double> addScalar(vectorShape);
  pvar<double> constAssignment(vectorShape, 8.0);

  //Print results
  cout <<"* Default init (check 0): " <<addScalar[5]<< endl;
  cout <<"* Constant init (check 8.0): " <<constAssignment[5]<< endl;


  // Index Assignment
  cout <<"\n===== Index Assignment Tests =====\n";
  pvar<double> p_id(vectorShape, dimensionName::x);
  constAssignment.addIndex(dimensionName::x);

  //Print results
  cout <<"* Dimension Index init (check 5.0): " <<p_id[5]<< endl;
  cout <<"* Dimension Index assignment (check 14.0): "
       <<constAssignment[6]<< endl;

  
  // pvar ops
  cout <<"\n===== Basic pvar ops Tests =====\n";
  pvar<double> p_add(vectorShape);
  pvar<double> p_subtr(vectorShape);
  pvar<double> p_mult(vectorShape);
  {
    constAssignment= 1.0;
    addScalar= p_id + 1.0;
    p_add= p_id + p_id;
    p_subtr= p_id - p_id;
    p_mult= p_id * p_id;
  }
  // pvar ops results
  cout <<"* Dimension Index assignment (check 1.0): "
       <<constAssignment[6]<< endl;
  cout <<"* ID test (check 5): " <<p_id[5]<< endl;
  cout <<"* Scalar Addition test (check 6): " <<addScalar[5]<< endl;
  cout <<"* Addition test (check 10): " <<p_add[5]<< endl;
  cout <<"* Subtraction test (check 0): " <<p_subtr[5]<< endl;
  cout <<"* Multiplication test (check 25): " <<p_mult[5]<< endl;


  // Test Random
  cout <<"\n===== Random() Tests =====\n";
  pvar<float> p_maxTest(vectorShape);  

  vectorShape-> initPseudoRandom();
  p_maxTest.setToRandom();

  //Print out results
  cout <<"* Random test (check 0.0 to 1.0): " <<p_maxTest[5] <<endl;

  
  // Test reductions
  cout <<"\n===== Reduction Tests =====\n";  

  // reduction *+ test
  double fmadd_result{0};
  fmadd_result= p_id.reduce_MultAdd(p_id);
  cout <<"* Reduce_MultAdd test: " <<fmadd_result <<endl;

  double sum=0.0;
  double if2=0.0;
  for (int i=0; i< vectorShape-> countElements(); ++i){
    if2=(double)i;
    sum+= if2*if2;
  }
  cout << " [Serial check: " << sum <<"]" <<endl;

  
  // Test reduce_MaxElement()
  // maxTest defined in random() tests above
  float max_result{0.0};
  max_result= p_maxTest.reduce_MaxElement();
  cout <<"* Max test (Range 0-1): " <<max_result <<endl;


  // Test reduce_Sum()
  pvar<int> p_sumTest(vectorShape, 1);
  int sum_result{0};
  sum_result= p_sumTest.reduce_Sum();
  cout <<"* Sum test (check " <<vectorShape->countElements()  <<" ): "
       <<sum_result <<endl;

  
  // Test function type mismatch error.
  cout <<"\n===== Invalid Condition Tests =====\n";  
  pvar<int> p_inttest(vectorShape, 1);
  cout <<"* Expect error for no FP implementation\n";
  p_inttest= p_sumTest.ln();


  // Test some of the uniary ops.
  cout <<"\n===== Uniary Op Tests =====\n";  

  // Create test data for uniary test  
  pvar<float> p_testData(vectorShape);
  p_testData= p_maxTest * 2.0 - 1.0;


  // Test uniary signum().
  pvar<float> p_signumTest(vectorShape);
  // Note test of double to float pvar
  p_signumTest= p_testData.signum();
  cout <<"* Parallel signum() test: " <<p_testData[5] <<" "
       <<p_signumTest[5] <<endl;
  
  cout <<"* Parallel signum() test: " <<p_testData[10] <<" "
       <<p_signumTest[10] <<endl;
  
  cout <<"* Parallel signum() test: " <<p_testData[15] <<" "
       <<p_signumTest[15] <<endl;

 
  // Test of abs() vs fabs()
  pvar<float> floatAbsTest = pvar<double>(vectorShape);
  pvar<float>floatTestData(vectorShape, dimensionName::x);
  floatTestData= floatTestData - 1000.001;
  floatAbsTest= floatTestData.abs();
  cout <<"* Parallel fabs test: " <<floatTestData[15] <<" "
       <<floatAbsTest[15] <<endl;

  pvar<int>iTestData(vectorShape, dimensionName::x);
  pvar<int>intAbsTest(vectorShape);
  iTestData= iTestData - 1000;
  intAbsTest= iTestData.abs();
  cout <<"* Parallel abs test: " <<iTestData[15] <<" "
       <<intAbsTest[15] <<endl;


  /// Test power(x)
  pvar<double> p_uniaryTest2(vectorShape, dimensionName::x);
  p_uniaryTest2= p_testData.power(2);
  cout <<"* Parallel Uniary test (power): " <<p_testData[4] <<" "
       <<p_uniaryTest2[4] <<endl;

  // Place uniary function to be test here.
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
    cout <<"* Parallel Uniary test (" << US(_UNIFUNC)
	 <<"): " <<p_testData[j] <<" "
	 <<p_uniaryTest[j] <<endl;
  }
  
  // Test some of the binary ops.
  cout <<"\n===== Binary Op Tests =====\n";  

  pvar<float> binaryTest_a(vectorShape);
  pvar<float> binaryTest_b(vectorShape);
  pvar<float> binaryTest_result(vectorShape);

  binaryTest_a.setToRandom();
  binaryTest_b.setToRandom();
  
  binaryTest_result= binaryTest_a.min(binaryTest_b);
  cout <<"* Parallel binary test (min): "
       <<binaryTest_a[14] <<", "
       <<binaryTest_b[14] <<"= "
       <<binaryTest_result[14] <<endl;

  pvar<int> binaryIntTest(vectorShape, dimensionName::x);
  pvar<int> binaryIntResult(vectorShape);  
  binaryIntResult = binaryIntTest.mod(11);

  cout <<"* Parallel binary test (mod): "
       <<binaryIntTest[14] <<"%11 = "
       <<binaryIntResult[14] <<endl;


  // Test some of the binary ops.
  cout <<"\n===== ActiveSet Tests =====\n";

  vectorShape-> resetActiveSet();

  vectorShape-> ge(binaryIntResult, 2);
  cout <<"* Parallel activeSet test (ge 2): "
       <<binaryIntResult[14] <<" check(1) "
       <<vectorShape-> activityStatus(14) <<endl;

  // Test Uniary operators with active set.
  //----------------------------------------------------------------------
  p_uniaryTest= 0;
  p_uniaryTest= p_testData._UNIFUNC();
  for (auto j=10; j<16;++j){
    cout <<"* Parallel Uniary test (" << US(_UNIFUNC)
	 <<"): " <<p_testData[j] <<" "
	 <<p_uniaryTest[j] <<endl;
  }
   

  vectorShape-> resetActiveSet();

  vectorShape-> eq(binaryIntResult, binaryIntResult );
  cout <<"* Parallel activityStatus test ([14] eq scalar): "
       <<binaryIntResult[14] <<" check(1) "
       <<vectorShape-> activityStatus(14) <<endl;
  cout <<"* Parallel activityStatus test ([12] eq scalar): "
       <<binaryIntResult[12] <<" check(0) "
       <<vectorShape-> activityStatus(12) <<endl;


  // subsetReduction *+ test
  fmadd_result=0;

  // Create an active set of even elements
  vectorShape-> resetActiveSet();
  binaryIntTest= binaryIntTest.mod(2);
  vectorShape-> eq(binaryIntTest, 1);
  fmadd_result= vectorShape-> subsetReduce_MultAdd(p_id, p_id);
  cout <<"* subsetReduce_MultAdd test: " <<fmadd_result <<endl;

  sum=0.0;
  if2=0.0;
  for (int i=0; i< vectorShape-> countElements(); ++i){
    if2=(double)i;
    sum+= if2*if2;
  }
  cout << " [Serial check: " << sum/2 <<"]" <<endl;



  // Test areAllElementsActive();
  vectorShape-> resetActiveSet();
  cout << "* Test reset(): All elements are active [1]: "
       <<vectorShape-> areAllElementsActive() <<endl;
  vectorShape-> gt(p_id, -10.0); // set all elements with pid > (-10)
  cout << "* Test Set all true: All elements are active [1]: "
       <<vectorShape-> areAllElementsActive() <<endl;

  vectorShape-> le(p_id, 10.0);
  float fResult{0};
  constAssignment= 1.0;
  fResult= vectorShape-> subsetReduce_Sum(constAssignment);

  cout << "* Sum of elements with id<=10, [sum=11], actual:"
       <<fResult <<endl;

  vectorShape-> gt(p_id, -10.0);
  cout << "* Test Set all true:All elements are active [0]: "
       <<vectorShape-> areAllElementsActive() <<endl;

  
  // TEST active set min here.
  fResult=0;
  vectorShape-> gt(p_id, 4.0);
  fResult= vectorShape-> subsetReduce_Min(p_id);
  cout << "* Min elements with 4<id<=10, [min=5], actual:"
       <<fResult <<endl;
  
  cout << "* Number of active elements [6]: "
       <<vectorShape-> countActiveElements() <<endl;
  cout << "* All elements are inactive [false,0]: "
       <<vectorShape-> areAllElementsInactive() <<endl;
  cout << "* All elements are active [false,0]: "
       <<vectorShape-> areAllElementsActive() <<endl;

  cout <<"\n===== Spread Communication Tests =====\n";  

  cout <<"\n--TEST Flood()\n";  
  pvar<int> spreadTest1dData(vectorShape);
  pvar<int> spreadTest1d(vectorShape);
  int fromElementNumber{5};

  vectorShape-> resetActiveSet();
  spreadTest1dData.setToRandom();
  spreadTest1d= 
    spreadTest1dData.flood(fromElementNumber, dimensionName::x );
  cout <<"* Element number " <<fromElementNumber 
       <<" *spread* on dimension x, 1d, from element 5 to all. " <<endl;
  cout <<"    Element 5: "<< spreadTest1d[fromElementNumber]
       <<" Element 45: " << spreadTest1d[45] <<endl;


  cout <<"\n--TEST scanAdd() and shiftWithWrap() \n";  
  #define _SHORT_VECTOR_LEN 11
  shape* shortVectorShape= new shape(_SHORT_VECTOR_LEN);
  //shortVectorShape-> resetActiveSet();
  
  pvar<int> scanTest1dData(shortVectorShape);
  pvar<int> scanTest1d(shortVectorShape);

  scanTest1dData= 1;
  scanTest1d= 
    scanTest1dData.scanAdd(dimensionName::x);

  scanTest1d= scanTest1d.shiftWithWrap( 1, dimensionName::x );

  for (auto i=0; i<(_SHORT_VECTOR_LEN); i++){
    cout <<"* Element number " <<i 
	 <<" shifted accumulated value: " << scanTest1d[i] <<endl;
  }
  cout<< "Check values: 11 [1-10]\n";


  // Enumeration test
  cout <<"\n-- TEST Enumerate even elements.\n";  
  pvar<int> enumTest(shortVectorShape);
  pvar<int> evenElements(shortVectorShape);

  // turn on only the even elements
  evenElements.setToIndex(dimensionName::x);
  evenElements= evenElements.mod(2);
  shortVectorShape-> eq(evenElements, 0);
  
  enumTest.enumerateActiveElements();
  for (auto i=0; i < ( shortVectorShape-> countElements() ); i++){
    if (shortVectorShape-> activityStatus(i)){
      cout <<"* Element number " <<i 
	   <<" enumeration value: " << enumTest[i] <<endl;
    }
  }

  

  return 0;
}
