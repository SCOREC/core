#ifndef LION_BASE64_H
#define LION_BASE64_H

#include <string>

namespace lion {

/*
Function getBase64Char:
  gets the char at the index given to it from the encode table
*/
char getBase64Char (int index);

// ===========================================================================

/*
Function getDecodedBase64Char:
  gets the Base64 value of an ASCII char sent to it. Uses the  
  base64DecodeTable. We have to do this because our Base64 string is
  represented as ASCII chars which have different codes than Base64.
*/
unsigned char getDecodedBase64Char (unsigned char c);

// ===========================================================================

/*
Function base64Encode3Bytes:
  Encodes 3 bytes sent to it into Base64

Arguments:
  char* bytes - char of array of bytes to be encoded - the first three positions
                will be accessed so len(bytes) must be >= 3

Returns:
  std::string - string of length 4 of the encoded bytes
*/
std::string base64Encode3Bytes (char* bytes);

// ===========================================================================

/*
Function base64Encode2Bytes:
  Encodes 2 bytes sent to it into Base64. Used a when the end of a string is 
  being encoded and the string is of a length not divisible by 3. Since Base64
  encodes in 3 byte sections a padding char, '=', must be added to take the 
  place of the missing byte.

Arguments:
  char* bytes - char array of bytes to be encoded, the first two positions
                will be accessed so len(bytes) must be >= 2

Returns:  
  std::string - string of length 4 of the encoded bytes including padding
*/
std::string base64Encode2Bytes (char* bytes);

// ===========================================================================

/*
Function base64Encode1Byte:
  Encodes a single byte sent to it into Base64. Used when the end of a string is
  being encode and the tsring is of a length not divisible by 3. Since Base64
  encodes in 3 byte sections, two padding chars, "==", must be added to take the
  place of the missing bits.

Arguments:
  char byte - char of the byte to be encoded

Returns:
  std::string - string of length 4 of the encoded bytes including padding 
*/
std::string base64Encode1Byte (char byte);

// ===========================================================================

/*
Function base64Encode:
  Encodes a series of bytes specified by the arguments into a string of Base64
  chars that will represent them

Arguments:
  char* input - pointer to start of byte string to be encoded
  long len - number of bytes to be encoded

Returns:
  std::string - Base64 encoded string of the input bytes
*/
std::string base64Encode (const char* input,
                          const unsigned long len );

// ===========================================================================

/*
Function base64Decode4Bytes:
  Decodes 4 bytes send to it from Base64 to plaintext,

Arguments:
  char* bytes - char array of bytes to be decoded into plaintext, the first
                four positions will be accessed so len(bytes) must be >= 4

Returns:
  std::string - string of length 1, 2, or 3 based on the amount of padding in
                the chars sent to the function. len(string) will be 1 if the
                input has two padding chars ('='*2). len(string) will be 2 if
                the input has one padding char ('='). len(string) will be 3
                if the input has no padding chars.
*/
std::string base64Decode4Bytes (char* bytes);

// ===========================================================================

/*
TODO: PR4: Implement

Function base64Decode:
  Decodes a string of Base64 chars to the bytes they represent

Arguments:
  std::string encoded - Base64 encoded string

Returns:
  std::string - decoded string of bytes from the encoded data
*/
std::string base64Decode (std::string encoded);

}

#endif
