#include <iostream>
#include <cassert>
#include <cstdlib>

#include <lionBase64.h>

using namespace lion;

void runEncodeTests () {

  //test1, encode3({'a','a','a'}) == "YWFh"
  char* testStr3 = (char*)malloc(3*sizeof(char));
  testStr3[0] = 'a';
  testStr3[1] = 'a';
  testStr3[2] = 'a';
  assert(base64Encode3Bytes(testStr3) == "YWFh");

  //test2, encode3({'x','y','z'}) == "eHl6"
  testStr3[0] = 'x';
  testStr3[1] = 'y';
  testStr3[2] = 'z';
  assert(base64Encode3Bytes(testStr3) == "eHl6");

  //test3, encode2({'a','a'}) == "YWE="
  char* testStr2 = (char*)malloc(2*sizeof(char));
  testStr2[0] = 'a';
  testStr2[1] = 'a';
  assert(base64Encode2Bytes(testStr2) == "YWE=");

  //test4, encode2({'x','y'}) == "eHk=""
  testStr2[0] = 'x';
  testStr2[1] = 'y';
  assert(base64Encode2Bytes(testStr2) == "eHk=");

  //test5, encode1('a') == "YQ=="
  char testStr;
  testStr = 'a';
  assert(base64Encode1Byte(testStr) == "YQ==");

  //test6, encode1('z') == "eg=="
  testStr = 'z';
  assert(base64Encode1Byte(testStr) == "eg==");

  //test7, encodeStr({'a','b','c','d','e','f'}) == "YWJjZGVm"
  char* testStr6 = (char*)malloc(6*sizeof(char));
  testStr6[0] = 'a';
  testStr6[1] = 'b';
  testStr6[2] = 'c';
  testStr6[3] = 'd';
  testStr6[4] = 'e';
  testStr6[5] = 'f';
  assert(base64Encode(testStr6, 6) == "YWJjZGVm");

  //test8, encodeStr({'u','v','w','x','y','z'}) == "dXZ3eHl6"
  testStr6[0] = 'u';
  testStr6[1] = 'v';
  testStr6[2] = 'w';
  testStr6[3] = 'x';
  testStr6[4] = 'y';
  testStr6[5] = 'z';
  assert(base64Encode(testStr6, 6) == "dXZ3eHl6");

  //test9, encodeStr({'a','b','c','d','e','f','g'}) == "YWJjZGVmZw=="
  char* testStr7 = (char*)malloc(7*sizeof(char));
  testStr7[0] = 'a';
  testStr7[1] = 'b';
  testStr7[2] = 'c';
  testStr7[3] = 'd';
  testStr7[4] = 'e';
  testStr7[5] = 'f';
  testStr7[6] = 'g';
  assert(base64Encode(testStr7, 7) == "YWJjZGVmZw==");

  //test10, encodeStr({'t','u','v','w','x','y','z'}) == "dHV2d3h5eg=="
  testStr7[0] = 't';
  testStr7[1] = 'u';
  testStr7[2] = 'v';
  testStr7[3] = 'w';
  testStr7[4] = 'x';
  testStr7[5] = 'y';
  testStr7[6] = 'z';
  assert(base64Encode(testStr7, 7) == "dHV2d3h5eg==");

  //test11, encodeStr({'a','b','c','d','e','f','g','h'}) == "YWJjZGVmZ2g="
  char* testStr8 = (char*)malloc(8*sizeof(char));
  testStr8[0] = 'a';
  testStr8[1] = 'b';
  testStr8[2] = 'c';
  testStr8[3] = 'd';
  testStr8[4] = 'e';
  testStr8[5] = 'f';
  testStr8[6] = 'g';
  testStr8[7] = 'h';
  assert(base64Encode(testStr8, 8) == "YWJjZGVmZ2g=");

  //test12, encodeStr({'s','t','u','v','w','x','y','z'}) == "c3R1dnd4eXo="
  testStr8[0] = 's';
  testStr8[1] = 't';
  testStr8[2] = 'u';
  testStr8[3] = 'v';
  testStr8[4] = 'w';
  testStr8[5] = 'x';
  testStr8[6] = 'y';
  testStr8[7] = 'z';
  assert(base64Encode(testStr8, 8) == "c3R1dnd4eXo=");

  free(testStr3);
  free(testStr2);
  free(testStr6);
  free(testStr7);
  free(testStr8);

  std::cout << "Encode tests pass!" << std::endl;
}

void runDecodeTests () {

  //test1, test decode table
  std::string base64CharsNoEquals = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "abcdefghijklmnopqrstuvwxyz"
                                    "0123456789+/";
  for ( unsigned int i = 0; i < base64CharsNoEquals.length(); i++ ) {
    unsigned int base64Code = getDecodedBase64Char(base64CharsNoEquals[i]);
    assert(i == base64Code);
  }

  //test2, decode4({'Y','W','F','h'}) == "aaa"
  char* testStr4 = (char*)malloc(4*sizeof(char));
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = 'F';
  testStr4[3] = 'h';
  assert(base64Decode4Bytes(testStr4) == "aaa");

  //test3, decode4({'Y','W','F','='}) == "aa"
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = 'F';
  testStr4[3] = '=';
  assert(base64Decode4Bytes(testStr4) == "aa");

  //test4, decode4({'Y','W','=','='}) == "a"
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = '=';
  testStr4[3] = '=';
  assert(base64Decode4Bytes(testStr4) == "a");

  //test5, decode({'e','H','l','6'}) == "xyz"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = 'l';
  testStr4[3] = '6';
  assert(base64Decode4Bytes(testStr4) == "xyz");

  //test6, decode({'e','H','l','='}) == "xy"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = 'l';
  testStr4[3] = '=';
  assert(base64Decode4Bytes(testStr4) == "xy");

  //test7, decode({'e','H','=','='}) == "x"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = '=';
  testStr4[3] = '=';
  assert(base64Decode4Bytes(testStr4) == "x");

  //test8, decode({'Y','W','J','j','Z','G','V','m'}) == "abcdef"
  // *** FUNCTION NOT YET IMPLEMENTED ***
  char* testStr8 = (char*)malloc(8*sizeof(char));
  testStr8[0] = 'Y';
  testStr8[1] = 'W';
  testStr8[2] = 'J';
  testStr8[3] = 'j';
  testStr8[4] = 'Z';
  testStr8[5] = 'G';
  testStr8[6] = 'V';
  testStr8[7] = 'm';
  // assert(base64Decode(testStr8) == "abcdef");

  //test9, decode({'d','X','Z','3','e','H','l','6'}) = "uvwxyz"
  // *** FUNCTION NOT YET IMPLEMENTED ***
  testStr8[0] = 'd';
  testStr8[1] = 'X';
  testStr8[2] = 'Z';
  testStr8[3] = '3';
  testStr8[4] = 'e';
  testStr8[5] = 'H';
  testStr8[6] = 'l';
  testStr8[7] = '6';
  // assert(base64Decode(testStr8) == "uvwxyz");

  //test10, decode({'d','n','d','4','e','X','o','='}) == "vwxyz"
  // *** FUNCTION NOT YET IMPLEMENTED ***
  testStr8[0] = 'd';
  testStr8[1] = 'n';
  testStr8[2] = 'd';
  testStr8[3] = '4';
  testStr8[4] = 'e';
  testStr8[5] = 'X';
  testStr8[6] = 'o';
  testStr8[7] = '=';
  // assert(base64Decode(testStr8) == "vwxyz");

  //test11, decode({'d','3','h','5','e','g','=','='}) = "wxyz"
  // *** FUNCTION NOT YET IMPLEMENTED ***
  testStr8[0] = 'd';
  testStr8[1] = '3';
  testStr8[2] = 'h';
  testStr8[3] = '5';
  testStr8[4] = 'e';
  testStr8[5] = 'g';
  testStr8[6] = '=';
  testStr8[7] = '=';
  // assert(base64Decode(testStr8) == "wxyz");


  free(testStr8);
  free(testStr4);

  std::cout << "Decode tests pass!" << std::endl;
}

int main (int argc, char** argv) {
    
  runEncodeTests();
  runDecodeTests();

  return 0;
}
