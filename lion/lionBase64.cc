#include "lionBase64.h"

#include <string>
#include <stdlib.h>

namespace lion {

// ===========================================================================
// =========================================================================== 
//  Implementation of base64.h. documentation of functions found in base64.h 
// ===========================================================================
// ===========================================================================

static const std::string base64EncodeTable =  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                              "abcdefghijklmnopqrstuvwxyz"
                                              "0123456789+/=";

// ===========================================================================

static const unsigned char base64DecodeTable[256] =
{ 
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  // ====================================
  0xFF,0xFF,0xFF,0x3E,0xFF,0xFF,0xFF,0x3F,
  0x34,0x35,0x36,0x37,0x38,0x39,0x3A,0x3B,
  0x3C,0x3D,0xFF,0xFF,0xFF,0x00,0xFF,0xFF,
  0xFF,0x00,0x01,0x02,0x03,0x04,0x05,0x06,
  0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,
  0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,
  0x17,0x18,0x19,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x20,
  0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x28,
  0x29,0x2A,0x2B,0x2C,0x2D,0x2E,0x2F,0x30,
  0x31,0x32,0x33,0xFF,0xFF,0xFF,0xFF,0xFF,
  // ====================================
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
  0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF
};

// ===========================================================================

char getBase64Char (int index) {

  return base64EncodeTable[index];
}

// ===========================================================================

unsigned char getDecodedBase64Char (unsigned char c) {

  return base64DecodeTable[c];
}

// ===========================================================================

std::string base64Encode3Bytes (char* bytes) {

  std::string encoded;

  //first 6 bits of byte1
  encoded += getBase64Char((bytes[0] >> 2) & 0x3F); 

  //last 2 bits of byte1 and the first 4 of byte2
  encoded += getBase64Char( ((bytes[0] << 4) & 0x30) | ((bytes[1] >> 4) & 0x0F));

  //last 4 bits of byte2 and the first 2 of byte3
  encoded += getBase64Char( ((bytes[1] << 2) & 0x3C) | ((bytes[2] >> 6) & 0x03));

  //last 6 of byte3
  encoded += getBase64Char( bytes[2] & 0x3F ); 

  return encoded;
}

// ===========================================================================

std::string base64Encode2Bytes (char* bytes) {

  std::string encoded;

  //first 6 bits of byte1
  encoded += getBase64Char((bytes[0] >> 2) & 0x3F); 

  //last 2 bits of byte1 combined with the the first 4 of byte2
  encoded += getBase64Char( ((bytes[0] << 4) & 0x30) | ((bytes[1] >> 4) & 0x0F));

  //last 4 bits of byte2 combined with the first 2 of byte3
  encoded += getBase64Char( ((bytes[1] << 2) & 0x3C) );

  //add a padding char since there is no byte 3
  encoded += '=';

  return encoded;
}

// ===========================================================================

std::string base64Encode1Byte (char byte) {

  std::string encoded;

  //first 6 bits of the byte
  encoded += getBase64Char((byte >> 2) & 0x3F); 

  //last 2 bits of the byte
  encoded += getBase64Char((byte << 4) & 0x30);

  //add the padding chars since there is no byte 2 or 3
  encoded += "==";

  return encoded;
}

// ===========================================================================

std::string base64Encode (const char* input,
                          const unsigned long len ) {

  std::string encoded;
  char* inputChars = (char*)malloc(3*sizeof(char));
  unsigned int index = 0;

  //encode all the string in 3 byte sections, this loop won't encode the last
  // 1 or 2 bytes if the len is not divisible by 3
  for ( ; index <= len - 3; index += 3 ) {
    inputChars[0] = input[index];
    inputChars[1] = input[index+1];
    inputChars[2] = input[index+2];
    encoded += base64Encode3Bytes(inputChars);
  }

  //encode the last 2 bytes if there are 2 bytes left at the end of the string
  // i.e. if (len % 3 == 2)
  if ( len - index == 2 ) {
    inputChars[0] = input[index];
    inputChars[1] = input[index+1];
    encoded += base64Encode2Bytes(inputChars);
    index += 2;
  }

  //encode the last 1 byte if there is one byte left at the end of the string
  // i.e. if (len % 3 == 1)
  if ( len - index == 1 ) {
    char inputChar = input[index];
    encoded += base64Encode1Byte(inputChar);
    index++;
  }

  free(inputChars);

  return encoded;
}

// ===========================================================================

std::string base64Decode4Bytes (char* bytes){
  
  std::string decoded;

  //get the Base64 codes for the input chars
  unsigned char* base64Chars = (unsigned char*)malloc(4*sizeof(char));
  base64Chars[0] = getDecodedBase64Char(bytes[0]);
  base64Chars[1] = getDecodedBase64Char(bytes[1]);
  base64Chars[2] = getDecodedBase64Char(bytes[2]);
  base64Chars[3] = getDecodedBase64Char(bytes[3]);

  //check that the input was all Base64 chars, getDecodedBase64Char returns
  // 0xFF if the char isn't a Base64 char
  if ( base64Chars[0] == 0xFF || base64Chars[1] == 0xFF
      || base64Chars[2] == 0xFF || base64Chars[3] == 0xFF ) {
    return ""; //input had non-Base64 chars
  }

  //convert the 6 bits of char1 and the first 2 bits of char2 to plaintext
  decoded += (char)(((base64Chars[0] << 2) & 0xFC) | ((base64Chars[1] >> 4) & 0x03));
  
  //check if there's padding instead of a Base64 char at the 2 index, decode if
  // there isn't
  if ( bytes[2] != '=' ) {
    decoded += (char)(((base64Chars[1] << 4) & 0xF0) | ((base64Chars[2] >> 2) & 0x0F));
  }

  //check if there's padding instead of a Base64 char at the 2 and 3 indicies, 
  // decode if there isn't
  if ( bytes[2] != '=' && bytes[3] != '=') {
    decoded += (char)(((base64Chars[2] << 6) & 0xC0) | (base64Chars[3] & 0x3F));
  }

  free(base64Chars);

  return decoded;
}

// ===========================================================================

std::string base64Decode (std::string encoded) {

  if ( encoded.length() % 4 != 0 ) {
    return "";
  }

  std::string decoded;
  char* charsToDecode = (char*)malloc(4*sizeof(char));

  for ( unsigned int index = 0; index < encoded.length(); index += 4 ) {
    charsToDecode[0] = encoded[index];
    charsToDecode[1] = encoded[index+1];
    charsToDecode[2] = encoded[index+2];
    charsToDecode[3] = encoded[index+3];
    decoded += base64Decode4Bytes(charsToDecode);
  }

  free(charsToDecode);
  
  return decoded;
}

}
