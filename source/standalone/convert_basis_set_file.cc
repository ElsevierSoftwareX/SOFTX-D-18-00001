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

/** @file convert_basis_set_file.cc

    \brief Program that can be used to convert a file downloaded from
    the EMSL Basis Set Library (in SuperMolecule format) to the format
    expected by the Ergo program.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <string>
#include <cassert>

static int getChargeForAtomName(const std::string & atomName) {
  if(atomName == "HYDROGEN") return 1;
  if(atomName == "HELIUM") return 2;
  if(atomName == "LITHIUM") return 3;
  if(atomName == "BERYLLIUM") return 4;
  if(atomName == "BORON") return 5;
  if(atomName == "CARBON") return 6;
  if(atomName == "NITROGEN") return 7;
  if(atomName == "OXYGEN") return 8;
  if(atomName == "FLUORINE") return 9;
  if(atomName == "NEON") return 10;
  if(atomName == "SODIUM") return 11;
  if(atomName == "MAGNESIUM") return 12;
  if(atomName == "ALUMINUM") return 13;
  if(atomName == "SILICON") return 14;
  if(atomName == "PHOSPHOROUS") return 15;
  if(atomName == "SULFUR") return 16;
  if(atomName == "CHLORINE") return 17;
  if(atomName == "ARGON") return 18;
  if(atomName == "POTASSIUM") return 19;
  if(atomName == "CALCIUM") return 20;
  if(atomName == "SCANDIUM") return 21;
  if(atomName == "TITANIUM") return 22;
  if(atomName == "VANADIUM") return 23;
  if(atomName == "CHROMIUM") return 24;
  if(atomName == "MANGANESE") return 25;
  if(atomName == "IRON") return 26;
  if(atomName == "COBALT") return 27;
  if(atomName == "NICKEL") return 28;
  if(atomName == "COPPER") return 29;
  if(atomName == "ZINC") return 30;
  if(atomName == "GALLIUM") return 31;
  if(atomName == "GERMANIUM") return 32;
  if(atomName == "ARSENIC") return 33;
  if(atomName == "SELENIUM") return 34;
  if(atomName == "BROMINE") return 35;
  if(atomName == "KRYPTON") return 36;
  return -1;
}

static std::string getSecondWordFromLine(const std::string & str) {
  int len = str.length();
  char s[len+1];
  strcpy(s, str.c_str());
  assert(len > 3);
  assert(s[0] == '$');
  assert(s[1] == ' ');
  int idx = 2;
  assert(s[idx] != ' ');
  while(idx < len) {
    if(s[idx] == ' ')
      break;
    idx++;
  }
  int nChars = idx - 2;
  char ss[nChars+1];
  memcpy(ss, &s[2], nChars);
  ss[nChars] = '\0';
  std::string resultStr = ss;
  return resultStr;
}

static bool is_digit(char c) {
  if(c >= '0' && c <= '9')
    return true;
  return false;
}

static bool checkIfLineHasThreeNumbers(const std::string & str) {
  int len = str.length();
  char s[len+1];
  strcpy(s, str.c_str());
  int nDigitsFound = 0;
  int idx = 0;
  while(idx < len) {
    if(s[idx] == ' ') {
      idx++;
      continue;
    }
    if(is_digit(s[idx])) {
      // Digit found. Check how many digits follow.
      int nDigits = 1;
      for(int k = 1; k < len; k++) {
	if(is_digit(s[idx+k]))
	  nDigits++;
	else
	  break;
      }
      idx += nDigits;
      nDigitsFound++;
    }
    else
      return false;
  }
  if(nDigitsFound == 3)
    return true;
  return false;
}

int main(int argc, char* argv[])
{
  printf("convert_basis_set_file 1.0\n");
  printf("Written by Elias Rudberg\n");
  printf("Source modified on Mon  9 Nov 13:22:39 CET 2015\n");
  if(argc != 3) {
    printf("usage: convert_basis_set_file infile outfile\n");
    return -1;
  }

  const char* inFileName = argv[1];
  const char* outFileName = argv[2];

  printf("inFileName = '%s', outFileName = '%s'\n", inFileName, outFileName);

  FILE* inFile = fopen(inFileName, "rb");
  if(!inFile) {
    printf("Error opening inFile '%s' for reading.\n", inFileName);
    return -1;
  }
  FILE* outFile = fopen(outFileName, "wb");
  if(!outFile) {
    printf("Error opening outFile '%s' for writing.\n", outFileName);
    return -1;
  }

  const int MAXFILESIZE = 8888888;
  std::vector<char> buf(MAXFILESIZE);
  memset(&buf[0], 0x00, MAXFILESIZE);
  
  if(fread(&buf[0], 1, MAXFILESIZE, inFile) <= 0) {
    printf("Error reading inFile\n");
    return -1;
  }
  if(buf[MAXFILESIZE-1] != '\0') {
    printf("Error: zero not found at end of buffer. File too large?\n");
    return -1;
  }

  // Count number of lines in file
  int nLines = 0;
  for(int i = 0; i < MAXFILESIZE; i++) {
    if(buf[i] == '\n')
      nLines++;
  }
  nLines++;

  std::vector<std::string> lines(nLines);

  const char* p = &buf[0];
  int lineCount = 0;
  while(*p != '\0') {
    // Find end of line
    const char* q = p;
    while(*q != '\n' && *q != '\0')
      q++;
    int nChars = q - p;
    char lineStr[nChars+1];
    memcpy(lineStr, p, nChars);
    lineStr[nChars] = '\0';
    lines[lineCount] = lineStr;
    lineCount++;
    p = q;
    if(*q == '\n')
      p++;
  }
  printf("lineCount = %d\n", lineCount);

  const std::string str_s = "$ S-TYPE FUNCTIONS";
  const std::string str_p = "$ P-TYPE FUNCTIONS";
  const std::string str_d = "$ D-TYPE FUNCTIONS";
  const std::string str_f = "$ F-TYPE FUNCTIONS";
  const std::string str_g = "$ G-TYPE FUNCTIONS";
  const std::string str_h = "$ H-TYPE FUNCTIONS";
  const std::string str_i = "$ I-TYPE FUNCTIONS";

  // Check how many atom types there are
  int nAtomTypes = 0;
  for(int i = 0; i < lineCount; i++) {
    std::string & currLine = lines[i];
    if(currLine == str_s)
      nAtomTypes++;
  }
  printf("nAtomTypes = %d\n", nAtomTypes);
  assert(nAtomTypes >= 1);
  
  std::vector<std::string> linesToInsert(nAtomTypes);
  int linesToInsertCount = 0;

  // OK, now we have extracted the lines.
  // Look for lines containing three integer numbers.
  int lineIdx = 0;
  while(lineIdx < lineCount) {
    std::string & currLine = lines[lineIdx];
    if(checkIfLineHasThreeNumbers(currLine)) {
      // Now previous line must be one of the following strings:
      assert(lineIdx > 5);
      std::string & prevLine = lines[lineIdx-1];
      if(prevLine != str_s &&
	 prevLine != str_p &&
	 prevLine != str_d &&
	 prevLine != str_f &&
	 prevLine != str_g &&
	 prevLine != str_h &&
	 prevLine != str_i) {
	printf("ERROR: string like 'X-TYPE FUNCTIONS' not found where expected.\n");
	return -1;
      }
      if(prevLine == str_s) {
	// Now we found a place where info about nuclear charge should be inserted.
	std::string & prevLine2 = lines[lineIdx-2];
	std::string atomName = getSecondWordFromLine(prevLine2);
	int charge = getChargeForAtomName(atomName);
	if(charge <= 0) {
	  printf("ERROR: getChargeForAtomName failed for atomName = '%s'\n", atomName.c_str());
	  return -1;
	}
	printf("atomName = '%s', charge = %d\n", atomName.c_str(), charge);
	char s[88];
	sprintf(s, "a %d", charge);
	std::string lineToInsert = s;
	linesToInsert[linesToInsertCount] = lineToInsert;
	linesToInsertCount++;
      }
    }
    lineIdx++;
  }
  assert(linesToInsertCount == nAtomTypes);

  int nLinesFinal = nLines + nAtomTypes;
  std::vector<std::string> linesFinal(nLinesFinal);
  // No go through all lines again, creating linesFinal.
  lineIdx = 0;
  int lineIdx2 = 0;
  int atomTypeCounter = 0;
  while(lineIdx < lineCount) {
    // Find next str_s line
    int foundIdx = -1;
    for(int idxTmp = lineIdx; idxTmp < lineCount; idxTmp++) {
      if(lines[idxTmp] == str_s) {
	foundIdx = idxTmp;
	break;
      }
    }
    assert(foundIdx >= 0);
    while(lineIdx < foundIdx)
      linesFinal[lineIdx2++] = lines[lineIdx++];
    linesFinal[lineIdx2++] = linesToInsert[atomTypeCounter++];
    linesFinal[lineIdx2++] = lines[lineIdx++];
    if(atomTypeCounter == nAtomTypes) {
      while(lineIdx < lineCount)
	linesFinal[lineIdx2++] = lines[lineIdx++];
    }
  }

  for(int i = 0; i < nLinesFinal; i++)
    fprintf(outFile, "%s\n", linesFinal[i].c_str());
  fclose(inFile);
  fclose(outFile);

  printf("Done, file '%s' created OK, nLinesFinal = %d.\n", outFileName, nLinesFinal);  
  
  return 0;
}
