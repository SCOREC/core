#include <iostream>
#include <pcu_util.h>

#include <lionBase64.h>

void runEncodeTests ()
{
  //test1, encode3({'a','a','a'}) == "YWFh"
  char testStr3[3];
  testStr3[0] = 'a';
  testStr3[1] = 'a';
  testStr3[2] = 'a';
  PCU_ALWAYS_ASSERT(lion::base64Encode3Bytes(testStr3) == "YWFh");

  //test2, encode3({'x','y','z'}) == "eHl6"
  testStr3[0] = 'x';
  testStr3[1] = 'y';
  testStr3[2] = 'z';
  PCU_ALWAYS_ASSERT(lion::base64Encode3Bytes(testStr3) == "eHl6");

  //test3, encode2({'a','a'}) == "YWE="
  char testStr2[2];
  testStr2[0] = 'a';
  testStr2[1] = 'a';
  PCU_ALWAYS_ASSERT(lion::base64Encode2Bytes(testStr2) == "YWE=");
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr2, 2) == "YWE=");

  //test4, encode2({'x','y'}) == "eHk=""
  testStr2[0] = 'x';
  testStr2[1] = 'y';
  PCU_ALWAYS_ASSERT(lion::base64Encode2Bytes(testStr2) == "eHk=");
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr2, 2) == "eHk=");

  //test5, encode1('a') == "YQ=="
  char testStr1[1];
  testStr1[0] = 'a';
  PCU_ALWAYS_ASSERT(lion::base64Encode1Byte(testStr1[0]) == "YQ==");
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr1, 1) == "YQ==");

  //test6, encode1('z') == "eg=="
  testStr1[0] = 'z';
  PCU_ALWAYS_ASSERT(lion::base64Encode1Byte(testStr1[0]) == "eg==");
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr1, 1) == "eg==");

  //test7, encodeStr({'a','b','c','d','e','f'}) == "YWJjZGVm"
  char testStr6[6];
  testStr6[0] = 'a';
  testStr6[1] = 'b';
  testStr6[2] = 'c';
  testStr6[3] = 'd';
  testStr6[4] = 'e';
  testStr6[5] = 'f';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr6, 6) == "YWJjZGVm");

  //test8, encodeStr({'u','v','w','x','y','z'}) == "dXZ3eHl6"
  testStr6[0] = 'u';
  testStr6[1] = 'v';
  testStr6[2] = 'w';
  testStr6[3] = 'x';
  testStr6[4] = 'y';
  testStr6[5] = 'z';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr6, 6) == "dXZ3eHl6");

  //test9, encodeStr({'a','b','c','d','e','f','g'}) == "YWJjZGVmZw=="
  char testStr7[7];
  testStr7[0] = 'a';
  testStr7[1] = 'b';
  testStr7[2] = 'c';
  testStr7[3] = 'd';
  testStr7[4] = 'e';
  testStr7[5] = 'f';
  testStr7[6] = 'g';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr7, 7) == "YWJjZGVmZw==");

  //test10, encodeStr({'t','u','v','w','x','y','z'}) == "dHV2d3h5eg=="
  testStr7[0] = 't';
  testStr7[1] = 'u';
  testStr7[2] = 'v';
  testStr7[3] = 'w';
  testStr7[4] = 'x';
  testStr7[5] = 'y';
  testStr7[6] = 'z';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr7, 7) == "dHV2d3h5eg==");

  //test11, encodeStr({'a','b','c','d','e','f','g','h'}) == "YWJjZGVmZ2g="
  char testStr8[8];
  testStr8[0] = 'a';
  testStr8[1] = 'b';
  testStr8[2] = 'c';
  testStr8[3] = 'd';
  testStr8[4] = 'e';
  testStr8[5] = 'f';
  testStr8[6] = 'g';
  testStr8[7] = 'h';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr8, 8) == "YWJjZGVmZ2g=");

  //test12, encodeStr({'s','t','u','v','w','x','y','z'}) == "c3R1dnd4eXo="
  testStr8[0] = 's';
  testStr8[1] = 't';
  testStr8[2] = 'u';
  testStr8[3] = 'v';
  testStr8[4] = 'w';
  testStr8[5] = 'x';
  testStr8[6] = 'y';
  testStr8[7] = 'z';
  PCU_ALWAYS_ASSERT(lion::base64Encode(testStr8, 8) == "c3R1dnd4eXo=");

  std::cout << "Encode tests pass!" << std::endl;
}

void runDecodeTests ()
{
  //test1, test decode table
  std::string base64CharsNoEquals = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "abcdefghijklmnopqrstuvwxyz"
                                    "0123456789+/";
  for ( unsigned int i = 0; i < base64CharsNoEquals.length(); i++ )
  {
    unsigned int base64Code =
        lion::getDecodedBase64Char(base64CharsNoEquals[i]);
    PCU_ALWAYS_ASSERT(i == base64Code);
  }

  //test2, decode4({'Y','W','F','h'}) == "aaa"
  char testStr4[4];
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = 'F';
  testStr4[3] = 'h';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "aaa");
  PCU_ALWAYS_ASSERT(lion::base64Decode("YWFh") == "aaa");

  //test3, decode4({'Y','W','F','='}) == "aa"
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = 'F';
  testStr4[3] = '=';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "aa");
  PCU_ALWAYS_ASSERT(lion::base64Decode("YWF=") == "aa");

  //test4, decode4({'Y','W','=','='}) == "a"
  testStr4[0] = 'Y';
  testStr4[1] = 'W';
  testStr4[2] = '=';
  testStr4[3] = '=';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "a");
  PCU_ALWAYS_ASSERT(lion::base64Decode("YW==") == "a");

  //test5, decode({'e','H','l','6'}) == "xyz"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = 'l';
  testStr4[3] = '6';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "xyz");
  PCU_ALWAYS_ASSERT(lion::base64Decode("eHl6") == "xyz");

  //test6, decode({'e','H','l','='}) == "xy"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = 'l';
  testStr4[3] = '=';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "xy");
  PCU_ALWAYS_ASSERT(lion::base64Decode("eHl=") == "xy");

  //test7, decode({'e','H','=','='}) == "x"
  testStr4[0] = 'e';
  testStr4[1] = 'H';
  testStr4[2] = '=';
  testStr4[3] = '=';
  PCU_ALWAYS_ASSERT(lion::base64Decode4Bytes(testStr4) == "x");
  PCU_ALWAYS_ASSERT(lion::base64Decode("eH==") == "x");

  //test8, decode("YWJjZGVm") == "abcdef"
  PCU_ALWAYS_ASSERT(lion::base64Decode("YWJjZGVm") == "abcdef");

  //test9, decode("dXZ3eHl6") = "uvwxyz"
  PCU_ALWAYS_ASSERT(lion::base64Decode("dXZ3eHl6") == "uvwxyz");

  //test10, decode("dnd4eXo=") == "vwxyz"
  PCU_ALWAYS_ASSERT(lion::base64Decode("dnd4eXo=") == "vwxyz");

  //test11, decode("d3h5eg==") = "wxyz"
  PCU_ALWAYS_ASSERT(lion::base64Decode("d3h5eg==") == "wxyz");

  std::cout << "Decode tests pass!" << std::endl;
}

void runCombinedTests()
{
  std::string resultStr;

  //test1, encode and decode the alphabet in lover case
  std::string alphabet0 = "abcdefghijklmnopqrstuvwxyz";
  resultStr =
      lion::base64Decode(
          lion::base64Encode( alphabet0.c_str(), alphabet0.length() )
      );
  PCU_ALWAYS_ASSERT(resultStr  == alphabet0);

  //test2, encode and decode more characters
  std::string alphabet1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                          "abcdefghijklmnopqrstuvwxyz"
                          "0123456789"
                          "!@#$^&*()-_=+[]{}:;<>,.?/|~`";
  resultStr =
      lion::base64Decode(
          lion::base64Encode(alphabet1.c_str(), alphabet1.length() )
      );

  PCU_ALWAYS_ASSERT(resultStr == alphabet1);

  std::cout << "Combined tests pass!" << std::endl;
}

int main ( )
{
  runEncodeTests();
  runDecodeTests();
  runCombinedTests();

  return 0;
}
